#!/usr/bin/env python3
"""Extract FASTA of unknowns from THAPBI PICT all-reads output."""
import argparse
import sys

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)

# Parse Command Line
usage = """\
The input file should be a THAPBI PICT 'read' report TSV file with columns
including MD5, <method>-predictions, Marker-sequence, Sample-count, and
Total-abundance.
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
    metavar="TSV",
    help="Input TSV read report filename, default stdin.",
)
parser.add_argument(
    "-o",
    "--output",
    dest="output",
    default="/dev/stdout",
    metavar="FASTA",
    help="Output FASTA filename, defaults stdout.",
)
parser.add_argument(
    "-s",
    "--samples",
    type=int,
    metavar="INTEGER",
    default=5,
    help=(
        "Minimum number of samples to require before using an unknown "
        "sequence. Default 5."
    ),
)
parser.add_argument(
    "-a",
    "--abundance",
    type=int,
    metavar="INTEGER",
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


def filter_unclassifed(input_filename, output_fasta, abundance, samples):
    """Extract FASTA file of unknown sequences."""
    with open(input_filename) as handle:
        with open(output_fasta, "w") as output:
            line = handle.readline().rstrip("\n")
            while line.startswith("#\t"):
                line = handle.readline().rstrip("\n")
            if not line.startswith("#Marker\tMD5\t"):
                sys.exit("ERROR: Invalid TSV input file.")
            header = [_.rstrip() for _ in line.split("\t")]
            md5_col = 1
            predictions_col = 2
            seq_col = header.index("Marker-sequence")
            samples_col = header.index("Sample-count")
            abundance_col = header.index("Total-abundance")
            for line in handle:
                parts = [_.rstrip() for _ in line.rstrip("\n").split("\t")]
                if len(parts) != len(header):
                    sys.exit("ERROR: Inconsistent field counts")
                if int(parts[samples_col]) < samples:
                    continue
                if int(parts[abundance_col]) < abundance:
                    continue
                if parts[seq_col] in ("", "-"):
                    continue
                if parts[predictions_col] not in ("", "-"):
                    continue
                output.write(
                    f">{parts[md5_col]}_{parts[abundance_col]} "
                    f"in {parts[samples_col]} samples\n{parts[seq_col]}\n"
                )


filter_unclassifed(options.input, options.output, options.abundance, options.samples)
