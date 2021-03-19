#!/usr/bin/env python3
"""Pool THAPBI PICT sample report using metadata."""
from __future__ import print_function

import sys
from optparse import OptionParser

import numpy as np

# from collections import Counter

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)

# Parse Command Line
usage = """Example usage:

$ python pooling.py -i input.tsv -o pooled -c 1,2,5

The input file should be a THAPBI PICT sample summary report as
a plain text tab separated variable (TSV) file. Give a filename
stem as the output argument (will make TSV and Excel files). The
column argument is which metadata to retain and group by.
"""

parser = OptionParser(usage=usage)
parser.add_option(
    "-i",
    "--input",
    dest="input",
    default="/dev/stdin",
    metavar="FILE",
    help="Input TSV filename.",
)
parser.add_option(
    "-c",
    "--columns",
    type=str,
    default="1",
    metavar="COLUMNS",
    help="Comma separated list (e.g, '1,3,5') of columns from the input "
    "table. Only these metadata columns will be retained, and will be "
    "used for pooling.",
)
parser.add_option(
    "-p",
    "--pending",
    type=str,
    default="0",
    metavar="COLUMN",
    help="Columns from the input table containing Y, Yes, True (in any case) "
    "where further sequencing is pending. Such entries will get a row of ??? "
    "(in addition to a row of values if any were non-zero).",
)
parser.add_option(
    "-o",
    "--output",
    dest="output",
    default="pooled",
    metavar="STEM",
    help="Output filename stem (defaults to 'pooled')",
)

options, args = parser.parse_args()
if args:
    sys.exit("ERROR: Invalid command line, try -h or --help.")


def pool(input_filename, output_stem, columns_str, column_pending):
    """Pool samples to make a more consise report."""
    try:
        value_cols = [int(_) - 1 for _ in columns_str.split(",")]
    except ValueError:
        sys.exit(
            "ERROR: Columns should be a comma separated list"
            f" of positive integers, not {columns_str!r}."
        )
    if min(value_cols) < 0:
        sys.exit("ERROR: Invalid column, should all be positive.")
    try:
        column_pending = int(column_pending) - 1
    except ValueError:
        sys.exit(
            "ERROR: Pending column should be a positive integer (or zero for none)."
        )
    if column_pending < 0:
        column_pending = None

    meta_samples = {}
    meta_species = {}
    meta_pending = {}

    with open(input_filename) as handle:
        line = handle.readline().rstrip("\n")
        if not line.startswith("#") or "\t" not in line:
            sys.exit("ERROR: Invalid TSV input file.")
        header = line.split("\t")
        try:
            sample_col = header.index("Sequencing sample")
            count_col = header.index("Seq-count")
            unk_col = header.index("Unknown")
        except IndexError:
            sys.exit("ERROR: Header does not match THAPBI PICT sample report.")
        if max(value_cols) >= min(sample_col, count_col, unk_col):
            sys.exit(
                f"ERROR: Requested column {max(value_cols)+1} not in metadata range."
            )
        if column_pending is not None and column_pending >= min(
            sample_col, count_col, unk_col
        ):
            sys.exit("ERROR: Pending column not in metadata range.")
        sp_headers = header[unk_col:]
        meta_headers = [header[_] for _ in value_cols]

        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                sys.exit("ERROR: Inconsistent field counts")
            meta = tuple(parts[_] for _ in value_cols)
            if column_pending is None:
                meta_pending[meta] = False
            else:
                meta_pending[meta] = parts[column_pending].upper().strip() in (
                    "Y",
                    "YES",
                    "TRUE",
                )
            samples = [_ for _ in parts[sample_col].split(";") if _ != "-"]
            sp_counts = parts[unk_col:]
            # If this has not been sequenced, treat as zero for summation
            sp_counts = np.array([0 if _ == "-" else int(_) for _ in sp_counts], np.int)
            if meta in meta_samples:
                meta_samples[meta].update(samples)
                meta_species[meta] += sp_counts
            else:
                meta_samples[meta] = set(samples)
                meta_species[meta] = sp_counts

    with open(output_stem + ".tsv", "w") as handle:
        handle.write(
            "\t".join(meta_headers)
            + "\tSamples-sequenced\t"
            + "\t".join(sp_headers)
            + "\n"
        )
        for meta, sp_counts in meta_species.items():
            if any(sp_counts) or not meta_pending[meta]:
                # If all zero AND pending, don't print this line - do ??? only
                handle.write(
                    "\t".join(meta)
                    + "\t"
                    + str(len(meta_samples[meta]))
                    + "\t"
                    + "\t".join(str(_) for _ in sp_counts)
                    + "\n"
                )
            if meta_pending[meta]:
                handle.write(
                    "\t".join(meta)
                    + "\t"
                    + str(len(meta_samples[meta]))
                    + "\t"
                    + "\t".join("?" for _ in sp_counts)
                    + "\n"
                )


pool(options.input, options.output, options.columns, options.pending)
