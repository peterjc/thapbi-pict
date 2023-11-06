#!/usr/bin/env python3
"""Convert GeneGenes style FASTA+TSV into SINTAX style annotated FASTA file."""
import argparse
import gzip
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)

# Parse Command Line
usage = """\
The input file should be a FASTA file without species annotation, and a simple
two-column tab separated plain text TSV file mapping the FASTA identifiers to
taxonomic information (Greengenes database style with entries like semicolon
g underscore underscore genus). The output is a FASTA file with SINTAX style
taxonomic annotation (like comma g colon genus).
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

# TODO: Optionally support extracting inputs from Qiime2 archive files (*.qza)?

parser = argparse.ArgumentParser(
    prog="gg_to_sintax.py",
    description="Combine GreenGeens style FASTA & TSV in SINTAX style FASTA.",
    epilog=usage,
)
parser.add_argument(
    "-i",
    "--input",
    metavar="FASTA",
    required=True,
    help="Input FASTA filename.",
)
parser.add_argument(
    "-t",
    "--taxonomy",
    metavar="TSV",
    required=True,
    help="Input TSV taxonomy, with column one matching the FASTA file IDs.",
)
parser.add_argument(
    "-o",
    "--output",
    dest="output",
    default="/dev/stdout",
    metavar="FASTA",
    help="SINTAX style FASTA output filename, defaults stdout.",
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


def merge(fasta_filename, tsv_filename, output_fasta):
    """Extract FASTA file of unknown sequences."""
    taxonomy = {}

    with gzip.open(tsv_filename, "rt") if tsv_filename.endswith(".gz") else open(
        tsv_filename
    ) as handle:
        for line in handle:
            idn, tax = line.rstrip().split("\t")
            taxonomy[idn] = convert_taxonomy(tax)

    with gzip.open(fasta_filename, "rt") if fasta_filename.endswith(".gz") else open(
        fasta_filename
    ) as handle:
        with open(output_fasta, "w") as output:
            for title, seq in SimpleFastaParser(handle):
                idn = title.split(None, 1)[0]
                try:
                    output.write(f">{idn};tax={taxonomy[idn]}\n{seq}\n")
                except KeyError:
                    sys.exit(f"ERROR: Missing taxonomy for: {title}")


merge(options.input, options.taxonomy, options.output)
