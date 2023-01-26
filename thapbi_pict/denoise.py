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
from math import floor
from math import log2

from Bio.SeqIO.FastaIO import SimpleFastaParser
from rapidfuzz.distance import Levenshtein
from rapidfuzz.process import extract
from rapidfuzz.process import extract_iter

from .utils import abundance_from_read_name
from .utils import md5seq


def unoise(counts, unoise_alpha=2.0, unoise_gamma=4, abundance_based=True, debug=False):
    """Apply UNOISE2 algorithm.

    Argument counts is an (unsorted) dict of sequences (for the same amplicon
    marker) as keys, with their total abundance counts as values.
    """
    debug = False  # too noisy otherwise
    if not counts:
        return {}
    if abundance_based:
        # size ordered abundance-based greedy clustering (AGC),
        # where choices are sorted by decreasing abundance.
        # Don't need to calculate all the distances:
        search_function = extract_iter
    else:
        # distance-based greedy clustering (DGC), must compute all distances
        # in order to sort on them. Can't use ``extractOne`` as the single
        # entry it picks may not pass the dynamic threshold.
        search_function = extract

    top_a = max(counts.values())  # will become first centroid
    last_a = None
    cutoff = 0
    centroids = defaultdict(set)
    high_abundance_centroids = None
    # Start by sorting sequences by abundance, largest first
    for a, query in sorted(
        ((a, seq) for (seq, a) in counts.items() if a >= unoise_gamma),
        key=lambda x: (-x[0], x[1]),
    ):
        # All the larger abundances have been processed already
        # Towards the end have lots of entries with the same abundance
        # so cache these dynamic thresholds use to narrow search pool
        if a != last_a:
            last_a = a
            # Fist centroid is largest, so gives upper bound on d
            cutoff = min(
                floor((log2(top_a / a) - 1) / unoise_alpha) if a * 2 < top_a else 0, 10
            )
            high_abundance_centroids = [
                seq for seq in centroids if counts[seq] >= a * 2 ** (unoise_alpha + 1)
            ]
        for choice, dist, _index in search_function(
            query,
            high_abundance_centroids,
            scorer=Levenshtein.distance,
            score_cutoff=cutoff,
        ):
            if a * 2 ** (unoise_alpha * dist + 1) <= counts[choice]:
                centroids[choice].add(query)
                if debug:
                    sys.stderr.write(
                        f"DEBUG: unoise C:{md5seq(choice)}_{counts[choice]} "
                        f"<-- Q:{md5seq(query)}_{a} dist {dist}\n"
                    )
                break
        else:
            # New centroid
            centroids[query].add(query)
            # print(query, "new")

    corrections = {}
    for seq, choices in centroids.items():
        for _ in choices:
            corrections[_] = seq
        assert corrections[seq] == seq, "Centroid missing"
    return corrections


def main(
    inputs,
    output,
    total_min_abundance=0,
    min_length=0,
    max_length=sys.maxsize,
    unoise_alpha=2.0,
    unoise_gamma=4,
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
    corrections = unoise(
        totals,
        unoise_alpha=unoise_alpha,
        unoise_gamma=unoise_gamma,
        debug=debug,
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
        # put the highest abundance entries first:
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
