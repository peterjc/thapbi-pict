#!/usr/bin/env python3
# Copyright 2022-2024 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Extract FASTA of unknowns from THAPBI PICT TSV output."""
import argparse
import sys

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.2")
    sys.exit(0)

# Parse Command Line
usage = """\
The input file should be either a THAPBI PICT classifier TSV file (with
columns including Marker/MD5_abundance, Sequence, and genus-species), or
'read' report TSV file (with columns including MD5, <method>-predictions,
Marker-sequence, Sample-count, and Total-abundance).
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
        "Default here is 1000, ten times the default used for the "
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
            while line.startswith("#") and not line.startswith(
                ("#Marker/MD5_abundance\t", "#Marker\tMD5\t")
            ):
                line = handle.readline().rstrip("\n")
            if line.startswith("#Marker/MD5_abundance\t"):
                # Tally file
                header = [_.rstrip() for _ in line.split("\t")]
                try:
                    seq_col = header.index("Sequence")
                except ValueError:
                    sys.exit("ERROR: Invalid TSV input (missing 'Sequence' column).")
                try:
                    predictions_col = header.index("genus-species")
                except ValueError:
                    sys.exit(
                        "ERROR: Invalid TSV input (missing 'genus-species' column)."
                    )
                assert seq_col < predictions_col
                # Assumes all columns between "Marker/MD5_abundance"
                # and "Sequence" are sample counts
                for line in handle:
                    parts = [_.rstrip() for _ in line.rstrip("\n").split("\t")]
                    if len(parts) != len(header):
                        sys.exit("ERROR: Inconsistent field counts")
                    marker_md5, total = parts[0].rsplit("_")
                    total = int(total)
                    marker, md5 = marker_md5.split("/")
                    counts = [int(_) for _ in parts[1:seq_col]]
                    assert sum(counts) == total, f"{parts[0]} vs sum {sum(counts)}"
                    counts = [_ for _ in counts if _]
                    if len(counts) < samples:
                        continue
                    if total < abundance:
                        continue
                    if parts[predictions_col] not in ("", "-"):
                        continue
                    output.write(
                        f">{parts[0]} in {len(counts)} samples\n{parts[seq_col]}\n"
                    )
            elif line.startswith("#Marker\tMD5\t"):
                # Read file
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
            else:
                hint = "\t".join(line.rsplit("\t")[0:2])
                sys.exit(
                    f"ERROR: Invalid TSV input. Didn't expect line starting: {hint}"
                )


filter_unclassifed(options.input, options.output, options.abundance, options.samples)
