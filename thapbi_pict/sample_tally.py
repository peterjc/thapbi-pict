# Copyright 2022 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Prepare a non-redundant TSV file using MD5 naming.

This implements the ``thapbi_pict sample-tally ...`` command.
"""
import gzip
import os
import sys
from collections import Counter
from collections import defaultdict

from Bio.SeqIO.FastaIO import SimpleFastaParser

from .prepare import load_fasta_header
from .utils import abundance_from_read_name
from .utils import file_to_sample_name
from .utils import md5seq


def main(
    inputs,
    output,
    min_abundance=0,
    min_length=0,
    max_length=sys.maxsize,
    gzipped=False,  # output
    debug=False,
):
    """Implement the ``thapbi_pict sample-tally`` command.

    Arguments min_length and max_length are applied while loading the input
    per-sample FASTA files.

    Argument min_abundance is applied to the output using the total per ASV
    (not to the per-sample abundance in individual input files).

    """
    if isinstance(inputs, str):
        inputs = [inputs]
    assert isinstance(inputs, list)
    assert inputs

    if os.path.isdir(output):
        sys.exit("ERROR: Output directory given, want a filename.")

    totals = Counter()
    counts = defaultdict(int)
    samples = set()
    sample_headers = {}
    for filename in inputs:
        # Assuming FASTA for now
        if debug:
            sys.stderr.write(f"DEBUG: Parsing {filename}\n")
        sample = file_to_sample_name(filename)
        assert sample not in samples, f"ERROR: Duplicate stem from {filename}"
        samples.add(sample)
        sample_headers[sample] = load_fasta_header(filename)
        assert "raw_fastq" in sample_headers[sample], sample_headers[sample]
        with open(filename) as handle:
            for _, seq in SimpleFastaParser(handle):
                seq = seq.upper()
                if min_length <= len(seq) <= max_length:
                    a = abundance_from_read_name(_.split(None, 1)[0])
                    totals[seq] += a
                    counts[seq, sample] += a

    if totals:
        sys.stderr.write(
            f"Loaded {len(totals)} unique sequences from {sum(totals.values())}"
            f" in total within length range, max abundance {max(totals.values())}\n"
        )
    else:
        sys.stderr.write("WARNING: Loaded zero sequences within length range\n")

    samples = sorted(samples)
    values = sorted(
        ((count, seq) for seq, count in totals.items() if count >= min_abundance),
        # put the highest abundance entries first:
        key=lambda x: (-x[0], x[1:]),
    )
    del totals

    # TODO - avoid double definition here and in summary code:
    stats_fields = (
        "Raw FASTQ",
        "Flash",
        "Marker",
        "Cutadapt",
        "Threshold pool",
        "Threshold",
        "Max non-spike",
        "Max spike-in",
        "Singletons",
    )

    if output == "-":
        if gzipped:
            raise ValueError("Does not support gzipped output to stdout.")
        out_handle = sys.stdout
    elif gzipped:
        out_handle = gzip.open(output, "wt")
    else:
        out_handle = open(output, "w")
    missing = [None] * len(samples)
    for stat in stats_fields:
        stat_values = [
            sample_headers[sample].get(stat.lower().replace(" ", "_"), None)
            for sample in samples
        ]
        if stat_values == missing:
            # Expected for "Max non-spike" and "Max spike-in" if no synthetic controls:
            # sys.stderr.write(f"WARNING: Missing all {stat} values in FASTA headers\n")
            continue
        if None in stat_values:
            sys.stderr.write(f"WARNING: Missing some {stat} values in FASTA headers\n")
        # Using "-" as missing value to match default in summary reports
        out_handle.write(
            "\t".join(
                ["#" + stat]
                + ["-" if _ is None else str(_) for _ in stat_values]
                + ["\n"]
            )
        )
        del stat_values
    out_handle.write("\t".join(["#MD5_abundance"] + samples + ["Sequence\n"]))
    for count, seq in values:
        data = "\t".join(str(counts[seq, sample]) for sample in samples)
        out_handle.write(f"{md5seq(seq)}_{count}\t{data}\t{seq}\n")
    if output != "-":
        out_handle.close()
