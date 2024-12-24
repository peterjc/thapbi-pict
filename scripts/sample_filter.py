#!/usr/bin/env python3
# Copyright 2023 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Pool THAPBI PICT sample report using metadata."""

import argparse
import re
import sys
from collections import Counter

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.1.0")
    sys.exit(0)

# Parse Command Line
usage = """\
The input file should be THAPBI PICT sample-tally TSV intermediate file(s),
optionally all with classifier output included (taxid and genus-species columns).
The output is a subset filtering the sample names according to the regex. e.g.

$ python sample_filter.py -i input.onebp.tsv -o subet.onebp.tsv -r "^N[00-99]-"

When used without a regex this script simply merges the input files.
"""

# TODO - Offer -a for minimum per seq per sample abundance threshold?

parser = argparse.ArgumentParser(
    prog="sample_filter.py",
    description="Select samples in a THAPBI PICT sample-tally TSV file with a regex",
    epilog=usage,
)
parser.add_argument(
    "-i",
    "--input",
    dest="input",
    metavar="FILE",
    required=True,
    nargs="+",
    help="Input TSV filename(s), at least one filename required.",
)
parser.add_argument(
    "-o",
    "--output",
    dest="output",
    default="-",
    metavar="STEM",
    help="Output filename stem (defaults to '-' meaning TSV only to stdout)",
)
parser.add_argument(
    "-r",
    "--regex",
    type=str,
    default="",
    metavar="REGEX",
    help="Regular expression applied to sample names.",
)
if len(sys.argv) == 1:
    sys.exit("ERROR: Invalid command line, try -h or --help.")
options = parser.parse_args()

if not options.regex and len(options.input) == 1:
    sys.stderr.write(
        "WARNING: With one input and no regex, you can just copy the file\n"
    )


def apply_filter(seqs, seq_meta, sample_headers, counts, compiled_regex):
    """Apply regex filter."""
    sample_headers = {
        sample: value
        for sample, value in sample_headers.items()
        if compiled_regex.match(sample)
    }
    # Don't actually need to drop the redundant entries in counts dict, but can:
    counts = {
        (idn, marker, sample): a
        for (idn, marker, sample), a in counts.items()
        if sample in sample_headers
    }
    # Can now potentially drop some sequences:
    totals = Counter()
    for (idn, marker, _sample), a in counts.items():
        if a:
            totals[idn, marker] += a
    seqs = {
        (idn, marker): seq for (idn, marker), seq in seqs.items() if totals[idn, marker]
    }
    seq_meta = {
        (idn, marker): value
        for (idn, marker), value in seq_meta.items()
        if totals[idn, marker]
    }
    sys.stderr.write(
        f"Regex reduced to {len(sample_headers)} samples and {len(seqs)} sequences\n"
    )
    return seqs, seq_meta, sample_headers, counts


def sample_filter(input_filenames, output_filename, regex):
    """Filter samples using a regular expression (regex)."""
    pattern = re.compile(regex) if regex else None

    try:
        from thapbi_pict.utils import parse_sample_tsv
    except ImportError:
        sys.exit(
            "ERROR: Couldn't import Python function thapbi_pict.utils.parse_sample_tsv"
        )

    all_seqs = {}
    all_seq_meta = {}
    all_sample_headers = {}
    all_counts = {}

    for input_filename in input_filenames:
        seqs, seq_meta, sample_headers, counts = parse_sample_tsv(input_filename)
        sys.stderr.write(
            f"Loaded {len(sample_headers)} samples and {len(seqs)} sequences"
            f" from {input_filename}\n"
        )
        if regex:
            seqs, seq_meta, sample_headers, counts = apply_filter(
                seqs, seq_meta, sample_headers, counts, pattern
            )
        all_seqs.update(seqs)
        all_seq_meta.update(seq_meta)
        all_sample_headers.update(sample_headers)
        all_counts.update(counts)

    try:
        from thapbi_pict.utils import export_sample_tsv
    except ImportError:
        sys.exit(
            "ERROR: Couldn't import Python function thapbi_pict.utils.export_sample_tsv"
        )

    export_sample_tsv(
        output_filename, all_seqs, all_seq_meta, all_sample_headers, all_counts
    )
    sys.stderr.write(
        f"Wrote {len(all_sample_headers)} samples and {len(all_seqs)} sequences"
        f" to {output_filename}\n"
    )


sample_filter(
    options.input,
    options.output,
    options.regex,
)
