# Copyright 2023 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Apply UNOISE read-correction to denoise FASTA file(s).

This implements the ``thapbi_pict denoise ...`` command, which is a simplified
version of the ``thapbi_pict sample-tally ...`` command intended to be easier
to use outside the THAPBI PICT pipeline.
"""
import gzip
import os
import sys
from collections import Counter
from collections import defaultdict

from Bio.SeqIO.FastaIO import SimpleFastaParser

from .sample_tally import unoise
from .utils import abundance_from_read_name
from .utils import md5seq


def main(
    inputs,
    output,
    total_min_abundance=0,
    min_length=0,
    max_length=sys.maxsize,
    gzipped=False,  # output
    debug=False,
):
    """Implement the ``thapbi_pict denoise`` command.

    This is a simplified version of the ``thapbi_pict sample-tally`` command
    which pools one or more FASTA input files before running the UNOISE read
    correction algorithm to denoise the dataset.

    Arguments min_length and max_length are applied while loading the input
    FASTA file(s).

    Arguments total_min_abundance is applied after read correction.

    Results sorted by decreasing abundance, then alphabetically by sequence.
    """
    if isinstance(inputs, str):
        inputs = [inputs]
    assert isinstance(inputs, list)
    assert inputs

    assert "-" not in inputs

    if os.path.isdir(output):
        sys.exit("ERROR: Output directory given, want a FASTA filename.")

    totals = Counter()
    for filename in inputs:
        if debug:
            sys.stderr.write(f"DEBUG: Parsing {filename}\n")
        with open(filename) as handle:
            for _, seq in SimpleFastaParser(handle):
                seq = seq.upper()
                if min_length <= len(seq) <= max_length:
                    a = abundance_from_read_name(_.split(None, 1)[0])
                    totals[seq] += a

    if totals:
        sys.stderr.write(
            f"Loaded {len(totals)} unique sequences from {sum(totals.values())}"
            f" in total within length range, max abundance {max(totals.values())}\n"
        )
    else:
        sys.stderr.write("WARNING: Loaded zero sequences within length range\n")

    if debug:
        sys.stderr.write("DEBUG: Starting UNOISE algorithm...\n")
    unoise_gamma = 4
    corrections = unoise(
        {seq: a for seq, a in totals.items() if a >= unoise_gamma},
        debug=False,
    )
    new_totals = defaultdict(int)
    for seq, a in totals.items():
        if totals[seq] < unoise_gamma:
            # Ignoring as per UNOISE algorithm
            assert seq not in corrections
            continue
        seq = corrections[seq]
        new_totals[seq] += a
    sys.stderr.write(
        f"UNOISE reduced unique ASVs from {len(totals)} to {len(new_totals)}, "
        f"max abundance now {max(new_totals.values(), default=0)}\n"
    )
    totals = new_totals
    del new_totals, corrections

    values = sorted(
        ((count, seq) for seq, count in totals.items() if count >= total_min_abundance),
        # (sort by marker), then put the highest abundance entries first:
        key=lambda x: (-x[0], x[1]),
    )
    del totals

    if output == "-":
        if gzipped:
            raise ValueError("Does not support gzipped output to stdout.")
        out_handle = sys.stdout
    elif gzipped:
        out_handle = gzip.open(output, "wt")
    else:
        out_handle = open(output, "w")

    for count, seq in values:
        md5 = md5seq(seq)
        out_handle.write(f">{md5}_{count}\n{seq}\n")
    if output != "-":
        out_handle.close()
