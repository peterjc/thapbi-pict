#!/usr/bin/env python3
"""Extract NCBI references matching unknown FASTA entries."""

import argparse
import sys
from collections import defaultdict

from Bio.SeqIO.FastaIO import SimpleFastaParser

from thapbi_pict.db_import import parse_ncbi_fasta_entry

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.4")
    sys.exit(0)

# Parse Command Line
usage = """\
The input file should be a FASTA file of observed amplicon sequences (e.g. the
all reads output, or your unknowns), with a reference FASTA file of published
untrimmed sequences expected to contain the amplicon, and optionally another
FASTA file which is a subset of those already trimmed.

Complete example as part of DB update:

$ ./Oomycota_ITS1_search.sh

That updates Oomycota_ITS1_search.fasta (all search results) and also
Oomycota_ITS1_w32.fasta (those with expected 32bp leader) which will be
included in the database.

Enare the observed NCBI entries file is empty, build a TEMP database:

$ rm -rf Oomycota_ITS1_obs.fasta; touch Oomycota_ITS1_obs.fasta
$ ./build_ITS1_DB.sh  # TEMP DB!

Then rerun the pipeline using this temp DB to get a full list of unknowns
without any previously obsered entries from NCBI:

$ ### run pipeline here with identity classifier! ###
$ ../scripts/unknowns.py -i thapbi_pict.ITS1.identity.tsv -o unknowns.fasta

$ ../scripts/missed_refs.py -i unknowns.fasta \
                   -f Oomycota_ITS1_search.fasta \
                   -x Oomycota_ITS1_w32.fasta \
                      Phytophthora_ITS1_curated.fasta \
                      Nothophytophthora_ITS1_curated.fasta \
                      Peronosporales_ITS1_curated.fasta \
                   -o Oomycota_ITS1_obs.fasta

$ ./build_ITS1_DB.sh  # Finished DB

i.e. Generate unknowns.fasta, then use it with this script to make
Oomycota_ITS1_obs.fasta for input to the database.
"""

parser = argparse.ArgumentParser(
    prog="missed_refs.py",
    description="Extract NCBI references matching unknown FASTA entries.",
    epilog=usage,
)
parser.add_argument(
    "-i",
    "--input",
    type=str,
    default="/dev/stdin",
    metavar="FILE",
    help="Input amplicons FASTA file. Required, default stdin.",
)
parser.add_argument(
    "-f",
    "--fasta",
    type=str,
    required=True,
    metavar="FILE",
    help="Input untrimmed reference FASTA file. Required.",
)
parser.add_argument(
    "-x",
    "--exclude",
    type=str,
    metavar="FILE",
    nargs="+",
    help="FASTA file of trimmed references to exclude. Optional.",
)
parser.add_argument(
    "-o",
    "--output",
    dest="output",
    default="/dev/stdout",
    type=str,
    metavar="FASTA",
    help="Output FASTA filename, defaults to stdout.",
)
parser.add_argument(
    "-l",
    "--left",
    default="TTTCCGTAGGTGAACCTGCGGAAGGATCATTA",
    metavar="SEQ",
    help="String to assume at the start of the amplicon, will be added to "
    "partial matches. Default TTTCCGTAGGTGAACCTGCGGAAGGATCATTA for ITS1.",
)

if len(sys.argv) == 1:
    sys.exit("ERROR: Invalid command line, try -h or --help.")
options = parser.parse_args()


def reject_title(title):
    """Check if a species name be rejected."""
    return title.split(None, 2)[1].lower() in (
        "uncultured",
        "unidentified",
        "environmental",
    )


def species_heuristics(text):
    """Apply NCBI species name heuristics."""
    # Ignore the (zero) taxid this returns
    return parse_ncbi_fasta_entry(text, [])[1]


def generate_references(
    input_fasta, reference_fasta, output_fasta, exclude_fastas, left, sep=";"
):
    """Extract FASTA file of reference sequences."""
    references = defaultdict(set)
    with open(reference_fasta) as handle:
        for title, seq in SimpleFastaParser(handle):
            if reject_title(title):
                continue
            # Take accession.X genus species (hopefully)
            references[seq.upper()].add(
                title.split(".", 1)[0] + ".X " + species_heuristics(title)
            )
    sys.stderr.write(f"Loaded {len(references)} references\n")
    exclude = set()
    if exclude_fastas:
        for exclude_fasta in exclude_fastas:
            sys.stderr.write(f"Will exclude entries from {exclude_fasta}\n")
            with open(exclude_fasta) as handle:
                for title, seq in SimpleFastaParser(handle):
                    if reject_title(title):
                        continue
                    exclude.add(seq.upper())
        sys.stderr.write(
            f"Will exclude {len(exclude)} refs from {len(exclude_fastas)} files\n"
        )

    drop_via_exclude = 0
    seq_dict = defaultdict(set)
    with open(input_fasta) as handle:
        for _, seq in SimpleFastaParser(handle):
            seq = seq.upper()
            if seq == (
                "TGAACCTGCGGAAGGATCATTACCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAA"
                "CCCTTTTATTGGGGGCTTCTGTCTGGTCTGGCTTCGGCTGGATTGGGTGGCGGCTCTATC"
                "ATGGCGACCGCTCTGAGCTTCGGCCTGGAGCTAGTAGCCCACTTTTTAAACCCATTCTTA"
                "ATTACTGAACAAACT"
            ):
                # Seems to have been truncated...
                continue
            if seq in exclude:
                drop_via_exclude += 1
                continue
            # sys.stderr.write(f"{title}\n")
            target = seq[len(left) :] if seq.startswith(left) else seq
            # Will try looking for matches which can be extended
            for ref in references:
                if target not in ref:
                    continue
                common = seq
                while common not in ref and len(left) + len(common) > len(seq):
                    common = common[1:]
                assert common in ref
                assert target in common
                assert len(left) + len(common) >= len(seq)
                masked = left[: len(seq) - len(common)].lower() + common
                assert (
                    masked.upper() == seq
                ), f"{left[:len(seq) - len(common)]} + {common} != {seq} from {ref}"
                seq_dict[masked].update(references[ref])
                # sys.stderr.write(f"{title} {references[ref]}\n")

    with open(output_fasta, "w") as out_handle:
        for seq in sorted(seq_dict, key=lambda s: (s.upper(), s)):
            out_handle.write(
                ">%s\n%s\n"
                % (
                    sep.join(sorted(seq_dict[seq])),
                    seq,
                )
            )

    if drop_via_exclude:
        sys.stderr.write(
            f"Dropped {drop_via_exclude} via the {len(exclude_fastas)} exclude files\n"
        )


generate_references(
    options.input, options.fasta, options.output, options.exclude, options.left
)
