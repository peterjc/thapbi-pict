# Copyright 2018-2020 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Prepare a non-redundant FASTA file using MD5 naming.

This implements the ``thapbi_pict fasta-nr ...`` command, and does part
of the work of the ``thapbi_pict prepare-reads`` command.
"""
import os
import sys
from collections import Counter

from Bio.Seq import reverse_complement
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .prepare import save_nr_fasta
from .utils import abundance_from_read_name


def main(
    inputs,
    revcomp,
    output,
    min_abundance=0,
    min_length=0,
    max_length=sys.maxsize,
    debug=False,
):
    """Implement the ``thapbi_pict fasta-nr`` command."""
    if not inputs:
        inputs = []
    if not revcomp:
        revcomp = []
    if isinstance(inputs, str):
        inputs = [inputs]
    if isinstance(revcomp, str):
        revcomp = [revcomp]
    assert isinstance(inputs, list)
    assert isinstance(revcomp, list)
    assert inputs or revcomp

    if os.path.isdir(output):
        if len(inputs + revcomp) == 1:
            output = os.path.join(output, os.path.basename((inputs + revcomp)[0]))
            if debug:
                sys.stderr.write(
                    "DEBUG: Single input with directory as ouput,"
                    f" writing to {output}\n"
                )
        else:
            sys.exit(
                "ERROR: Output directory can only be used with a single input file"
            )

    counts = Counter()
    for filename in inputs:
        # Assuming FASTA for now
        if debug:
            sys.stderr.write(f"DEBUG: Parsing {filename}\n")
        with open(filename) as handle:
            for _, seq in SimpleFastaParser(handle):
                if min_length <= len(seq) <= max_length:
                    a = abundance_from_read_name(_.split(None, 1)[0])
                    counts[seq.upper()] += a
    for filename in revcomp:
        if debug:
            sys.stderr.write(f"DEBUG: Parsing {filename} (will reverse complement)\n")
        with open(filename) as handle:
            for _, seq in SimpleFastaParser(handle):
                if min_length <= len(seq) <= max_length:
                    a = abundance_from_read_name(_.split(None, 1)[0])
                    counts[reverse_complement(seq.upper())] += a

    if counts:
        sys.stderr.write(
            f"Loaded {len(counts)} unique sequences from {sum(counts.values())}"
            f" in total within length range, max abundance {max(counts.values())}\n"
        )
    else:
        sys.stderr.write("WARNING: Loaded zero sequences within length range\n")

    accepted_total, accepted_count = save_nr_fasta(
        counts, output, min_abundance=min_abundance
    )
    sys.stderr.write(f"Saved {accepted_count} unique sequences\n")
