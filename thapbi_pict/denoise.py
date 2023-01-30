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
import tempfile
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
from .utils import run
from .versions import check_tools


def unoise(
    counts, unoise_alpha=2.0, unoise_gamma=4, abundance_based=False, debug=False
):
    """Apply UNOISE2 algorithm.

    Argument counts is an (unsorted) dict of sequences (for the same amplicon
    marker) as keys, with their total abundance counts as values.
    """
    debug = False  # too noisy otherwise
    if not counts:
        return {}

    top_a = max(counts.values())  # will become first centroid
    last_a = None
    cutoff = 0
    centroids = defaultdict(set)
    high_abundance_centroids = None
    if abundance_based:
        sys.stderr.write("Starting UNOISE abundance-based greedy clustering (AGC)\n")
    else:
        sys.stderr.write("Starting UNOISE distance-based greedy clustering (DGC)\n")
    # Start by sorting sequences by abundance, largest first
    for a, query in sorted(
        ((a, seq) for (seq, a) in counts.items() if a >= unoise_gamma),
        key=lambda x: (-x[0], x[1]),
    ):
        if debug:
            sys.stderr.write(f"DEBUG: UNOISE on {query}_{a}\n")
        # All the larger abundances have been processed already
        # Towards the end have lots of entries with the same abundance
        # so cache these dynamic thresholds use to narrow search pool
        if a != last_a:
            last_a = a
            # Fist centroid is largest, so gives upper bound on d
            cutoff = min(
                floor((log2(top_a / a) - 1) / unoise_alpha) if a * 2 < top_a else 0, 10
            )
            # dist = 1 is a lower bound
            high_abundance_centroids = [
                seq for seq in centroids if counts[seq] >= a * 2 ** (unoise_alpha + 1)
            ]
        if abundance_based:
            # size ordered abundance-based greedy clustering (AGC),
            # where choices are sorted by decreasing abundance.
            # Don't need to calculate all the distances:
            for choice, dist, _index in extract_iter(
                query,
                high_abundance_centroids,
                scorer=Levenshtein.distance,
                score_cutoff=cutoff,
            ):
                if a * 2 ** (unoise_alpha * dist + 1) <= counts[choice]:
                    centroids[choice].add(query)
                    if debug:
                        sys.stderr.write(
                            f"DEBUG: unoise-agc C:{md5seq(choice)}_{counts[choice]} "
                            f"<-- Q:{md5seq(query)}_{a} dist {dist}\n"
                        )
                    break
            else:
                # New centroid
                centroids[query].add(query)
                # print(query, "new")
        else:
            # distance-based greedy clustering (DGC), must compute all distances
            # in order to sort on them. Can't use ``extractOne`` as the single
            # entry it picks may not pass the dynamic threshold.
            candidates = sorted(
                [
                    (dist, choice)
                    for (choice, dist, _index) in extract(
                        query,
                        high_abundance_centroids,
                        scorer=Levenshtein.distance,
                        score_cutoff=cutoff,
                    )
                    if a * 2 ** (unoise_alpha * dist + 1) <= counts[choice]
                ],
                key=lambda x: (x[0], counts[x[1]], x[1]),
            )
            if candidates:
                if debug:
                    sys.stderr.write(
                        f"DEBUG: unoise-dgc {md5seq(query)} has "
                        f"{len(candidates)} candidates\n"
                    )
                    if len(candidates) > 1:
                        for i, (dist, choice) in enumerate(candidates):
                            sys.stderr.write(
                                f"DEBUG: #{i+1} {md5seq(choice)}_{counts[choice]} "
                                f"dist {dist}\n"
                            )
                # Take the closest one - tie break on abundance then sequence
                dist, choice = candidates[0]
                centroids[choice].add(query)
                if debug:
                    sys.stderr.write(
                        f"DEBUG: unoise-dgc C:{md5seq(choice)}_{counts[choice]} "
                        f"<-- Q:{md5seq(query)}_{a} dist {dist}\n"
                    )
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


def vsearch(
    counts,
    unoise_alpha=2.0,
    unoise_gamma=4,
    abundance_based=False,
    tmp_dir=None,
    debug=False,
    cpu=0,
):
    """Invoke VSEARCH to run its reimplementation of the UNOISE3 algorithm.

    Argument counts is an (unsorted) dict of sequences (for the same amplicon
    marker) as keys, with their total abundance counts as values.
    """
    check_tools(["vsearch"], debug)

    if tmp_dir:
        # Up to the user to remove the files
        tmp_obj = None
        shared_tmp = tmp_dir
    else:
        tmp_obj = tempfile.TemporaryDirectory()
        shared_tmp = tmp_obj.name

    if debug:
        sys.stderr.write(f"DEBUG: Shared temp folder {shared_tmp}\n")

    input_fasta = os.path.join(shared_tmp, "noisy.fasta")
    output_tsv = os.path.join(shared_tmp, "clustering.tsv")
    md5_to_seq = {}
    with open(input_fasta, "w") as handle:
        # Output using MD5 with USEARCH naming: md5;size=abundance
        prev_a = None
        for seq, a in counts.items():
            if a >= unoise_gamma:
                md5 = md5seq(seq)
                handle.write(f">{md5};size={a}\n{seq}\n")
                md5_to_seq[md5] = seq
            if prev_a is not None and prev_a < a:
                sys.exit("ERROR: Input to VSEARCH was not sorted by abundance")
            prev_a = a
        del prev_a

    cmd = [
        "vsearch",
        "--unoise_alpha",
        str(unoise_alpha),
        "--minsize",
        str(unoise_gamma),
        "--sizein",
        "--sizeout",
        "--cluster_unoise",
        input_fasta,
        # "--centroids",
        # output_fasta,
        "--uc",
        output_tsv,
    ]
    if abundance_based:
        cmd += ["--sizeorder", "--maxaccepts", "30"]
        # Note --sizeorder is documented to only takes effect if --maxaccepts
        # is higher than one, then it turns on abundance-based greedy
        # clustering (AGC), in contrast to the default distance-based greedy
        # clustering (DGC).
    if cpu:
        cmd += ["--threads", str(cpu)]
    run(cmd, debug=debug)
    corrections = {}
    with open(output_tsv) as handle:
        for line in handle:
            if not line:
                continue
            parts = line.split("\t")
            md5 = parts[8].split(";", 1)[0]  # strip the size information
            seq = md5_to_seq[md5]
            if len(parts) != 10:
                sys.exit(
                    f"ERROR: Found {len(parts)} fields in vsearch --uc output, not 10"
                )
            if parts[0] == "S":
                # centroid
                corrections[seq] = seq
            elif parts[0] == "H":
                # hit, mapped to a centroid
                centroid_md5 = parts[9].split(";", 1)[0]
                corrections[seq] = md5_to_seq[centroid_md5]
            elif parts[0] == "C":
                # cluster summary (should be one for each S line)
                pass
            else:
                sys.exit(
                    f"ERROR: vsearch --uc output line started {parts[0]}, not S, H or C"
                )
    del md5_to_seq
    return corrections


def read_correction(
    algorithm,
    counts,
    unoise_alpha=2.0,
    unoise_gamma=4,
    abundance_based=False,
    tmp_dir=None,
    debug=False,
    cpu=0,
):
    """Apply builtin UNOISE algorithm or invoke an external tool like VSEARCH.

    Argument counts is an (unsorted) dict of sequences (for the same amplicon
    marker) as keys, with their total abundance counts as values.
    """
    if algorithm == "unoise":
        # Does not need tmp_dir, cpu
        return unoise(counts, unoise_alpha, unoise_gamma, abundance_based, debug=False)
    elif algorithm == "vsearch":
        return vsearch(
            counts, unoise_alpha, unoise_gamma, abundance_based, tmp_dir, debug, cpu
        )
    elif algorithm == "-":
        sys.exit("ERROR: Internal function denoise_algorithm called without algorithm.")
    else:
        sys.exit(
            f"ERROR: denoise_algorithm called with {algorithm} (unknown algorithm)."
        )


def main(
    inputs,
    output,
    denoise_algorithm,
    total_min_abundance=0,
    min_length=0,
    max_length=sys.maxsize,
    unoise_alpha=2.0,
    unoise_gamma=4,
    gzipped=False,  # output
    tmp_dir=None,
    debug=False,
    cpu=0,
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
    assert denoise_algorithm != "-"

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
        sys.stderr.write(
            f"DEBUG: Starting read-correction with {denoise_algorithm}...\n"
        )
    corrections = read_correction(
        denoise_algorithm,
        totals,
        unoise_alpha=unoise_alpha,
        unoise_gamma=unoise_gamma,
        tmp_dir=tmp_dir,
        debug=debug,
        cpu=cpu,
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
