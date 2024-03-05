#!/usr/bin/env python3
"""Convert GreenGenes style FASTA+TSV into SINTAX style annotated FASTA file.

As of v0.2.0 of the script, the input files can optionally be provided inside
Qiime QZA files (ZIP files with a single in a .../data/... subdirectory).
"""

import argparse
import gzip
import io
import sys
import zipfile

from Bio.SeqIO.FastaIO import SimpleFastaParser

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.2.0")
    sys.exit(0)

# Parse Command Line
usage = """\
The input files should be a FASTA file without species annotation, and a simple
two-column tab separated plain text TSV file mapping the FASTA identifiers to
taxonomic information (GreenGenes database style with entries like semicolon
g underscore underscore genus). These can be provided inside Qiime QZA/QZV files
which are expected to contain a single data file each, here we are expecting
'.../data/dna-sequences.fasta' and '.../data/taxonomy.tsv' respectively. The
output is a FASTA file with SINTAX style taxonomic annotation (like comma g
colon genus).
"""

# Example FASTA:
#
#    >GBA28166-15
#    AACATTATATTTAATTTTCGGTGTCTGGGCAGGCCTGATTGGCACATCCTTAAGAATTTT...
#
# Example TSV:
#
#   Feature ID (tab) Taxon
#   GBA28166-15 (tab) p__Arthropoda;c__Insecta;o__Archaeognatha;f__Machilidae;
#   sf__Petrobiinae;g__Pedetontinus;s__Pedetontinus luanchuanensis
#
# Example of desired output FASTA using SINTAX style taxonomic annotation:
#
#    >GBA28166-15;tax=p:Arthropoda,c:Insecta,o:Archaeognatha,f:Machilidae,
#    g:Pedetontinus,s:Pedetontinus_luanchuanensis
#    AACATTATATTTAATTTTCGGTGTCTGGGCAGGCCTGATTGGCACATCCTTAAGAATTTT...
#
# Note SINTAX (and the NCBI taxonomy) does not have a sub-family rank.

parser = argparse.ArgumentParser(
    prog="gg_to_sintax.py",
    description="Combine GreenGeens style FASTA & TSV in SINTAX style FASTA.",
    epilog=usage,
)
parser.add_argument(
    "-i",
    "--input",
    metavar="FILE",
    required=True,
    help="Input FASTA filename, or QZA/QZV file containing "
    "'.../data/dna-sequences.fasta'.",
)
parser.add_argument(
    "-t",
    "--taxonomy",
    metavar="FILE",
    required=True,
    help="Input TSV taxonomy, or QZA/QZV file containing "
    "'.../data/taxonomy.tsv', with column one matching the FASTA file IDs.",
)
parser.add_argument(
    "-o",
    "--output",
    dest="output",
    default="/dev/stdout",
    metavar="FASTA",
    help="SINTAX style FASTA output filename, default stdout.",
)

if len(sys.argv) == 1:
    sys.exit("ERROR: Invalid command line, try -h or --help.")
options = parser.parse_args()


def convert_taxonomy(gg_style):
    """Convert GeneGenes entry with double under score etc into SINTAX style."""
    return (
        gg_style.replace("__", ":")
        .replace("; ", ",")
        .replace(";", ",")
        .replace(" ", "_")
    )


def qza_open(qza_filename, mode="rt", suffix=None):
    """Open the only data file within a QZA/QZV style ZIP file for reading."""
    if mode not in ("r", "rb", "rt"):
        raise ValueError(f"Unsupported mode {mode!r}, expected r, rb or rt only.")
    zip = zipfile.ZipFile(qza_filename)
    data = [_ for _ in zip.namelist() if "/data/" in _]
    if suffix:
        data = [_ for _ in data if _.endswith(suffix)]
    if len(data) != 1:
        if suffix:
            sys.exit(
                "ERROR: Expected one .../data/... file "
                "ending %r in Quiime file, found %i:\n%s"
                % (suffix, len(data), "\n".join(data))
            )
        else:
            sys.exit(
                "ERROR: Expected one .../data/... file "
                "in Quiime file, found %i:\n%s" % (len(data), "\n".join(data))
            )
    name = data[0]
    sys.stderr.write(f"Opening {name} from Qiime file.\n")
    if mode == "rb":
        return zip.open(data[0])
    else:
        return io.TextIOWrapper(zip.open(data[0]), encoding="utf-8")


def merge(fasta_filename, tsv_filename, output_fasta):
    """Extract FASTA file of unknown sequences."""
    taxonomy = {}

    with gzip.open(tsv_filename, "rt") if tsv_filename.endswith(".gz") else qza_open(
        tsv_filename, "rt", ".tsv"
    ) if tsv_filename.endswith((".qza", ".qzv")) else open(tsv_filename) as handle:
        for line in handle:
            idn, tax = line.rstrip().split("\t")
            taxonomy[idn] = convert_taxonomy(tax)

    with gzip.open(fasta_filename, "rt") if fasta_filename.endswith(
        ".gz"
    ) else qza_open(fasta_filename, "rt", ".fasta") if fasta_filename.endswith(
        (".qza", ".qzv")
    ) else open(fasta_filename) as handle:
        with open(output_fasta, "w") as output:
            for title, seq in SimpleFastaParser(handle):
                idn = title.split(None, 1)[0]
                try:
                    output.write(f">{idn};tax={taxonomy[idn]}\n{seq}\n")
                except KeyError:
                    sys.exit(f"ERROR: Missing taxonomy for: {title}")


merge(options.input, options.taxonomy, options.output)
