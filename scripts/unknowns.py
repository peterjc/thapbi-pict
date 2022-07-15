#!/usr/bin/env python3
"""Extract FASTA of unknowns from THAPBI PICT all-reads output."""
import argparse
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser

from thapbi_pict.utils import abundance_from_read_name

# from collections import Counter

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)

# Parse Command Line
usage = """\
The input file should be a THAPBI PICT 'intermediate' TSV file
with at least three columns:

1. Sequence name in format <MD5>_<count>
2. List of classified NCBI taxids
3. List of classified species names

and a matching FASTA file using the same sequence names.
"""

parser = argparse.ArgumentParser(
    prog="unknowns.py",
    description="Extract FASTA file of sequences with no classifier output",
    epilog=usage,
)
parser.add_argument(
    "-i",
    "--input",
    default="/dev/stdin",
    metavar="FILE",
    help="Input TSV filename, default stdin.",
)
parser.add_argument(
    "-f",
    "--fasta",
    required=True,
    metavar="FILE",
    help="Input FASTA file. Required.",
)
parser.add_argument(
    "-o",
    "--output",
    dest="output",
    default="/dev/stdout",
    metavar="FASTA",
    help="Output FASTA filename (defaults to stdout)",
)
parser.add_argument(
    "-a",
    "--abundance",
    type=int,
    default=1000,
    help=(
        "Minimum abundance to require before importing a sequence, "
        "over-and-above whatever was used to prepare the FASTA file. "
        "Default here is 100, ten times the default used for the "
        "classification pipeline - be cautious what goes in your marker DB."
    ),
)

if len(sys.argv) == 1:
    sys.exit("ERROR: Invalid command line, try -h or --help.")
options = parser.parse_args()


def filter_unclassifed(input_filename, input_fasta, output_fasta, abundance):
    """Extract FASTA file of unknown sequences."""
    wanted = set()
    with open(input_filename) as handle:
        line = handle.readline().rstrip("\n")
        if not line.startswith("#") or "\t" not in line:
            sys.exit("ERROR: Invalid TSV input file.")
        header = [_.rstrip() for _ in line.split("\t")]
        # e.g.
        # #ITS1/sequence-name (tab) taxid (tab) genus-species
        if (
            not header[0].endswith("/sequence-name")
            or not header[1] == "taxid"
            or not header[2] == "genus-species"
        ):
            sys.exit("ERROR: Header does not match THAPBI PICT classifier TSV output.")

        for line in handle:
            parts = [_.rstrip() for _ in line.rstrip("\n").split("\t")]
            if len(parts) != len(header):
                sys.exit("ERROR: Inconsistent field counts")
            if abundance_from_read_name(parts[0]) < abundance:
                continue
            if not parts[2]:
                wanted.add(parts[0])
                # sys.stderr.write(f"DEBUG: {parts[0]} is unknown\n")

    sys.stderr.write(
        f"Found {len(wanted)} entries without a species from the classifier\n"
    )

    with open(input_fasta) as handle:
        with open(output_fasta, "w") as output:
            for title, seq in SimpleFastaParser(handle):
                if title.split(None, 1)[0] in wanted:
                    output.write(f">{title}\n{seq}\n")


filter_unclassifed(options.input, options.fasta, options.output, options.abundance)
