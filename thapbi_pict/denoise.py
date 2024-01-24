# Copyright 2023-2024 by Peter Cock, The James Hutton Institute.
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
from time import time

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

    If not specified (i.e. set to zero or None), unoise_alpha defaults to 2.0
    and unoise_gamma to 4.

    Returns a dict mapping input sequences to centroid sequences, and an empty
    dict (no chimera detection performed).
    """
    debug = False  # too noisy otherwise
    if not counts:
        return {}, {}
    if not unoise_alpha:
        unoise_alpha = 2.0
    if not unoise_gamma:
        unoise_gamma = 4

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
                # UNOISE merges query into centroid if skew(M,C) <= beta(d),
                #   (query abundance) / (centroid abundance) <= 1 / 2**(alpha*dist +1)
                #   (query abundance) * 2**(alpha*dist + 1) <= centroid abundance
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
                        limit=None,
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
    return corrections, {}


def usearch(
    counts,
    unoise_alpha=None,
    unoise_gamma=None,
    abundance_based=False,
    tmp_dir=None,
    debug=False,
    cpu=0,
):
    """Invoke USEARCH to run its implementation of the UNOISE3 algorithm.

    Assumes v10 or v11 (or later if the command line API is the same).
    Parses the four columns tabbed output.

    Returns a dict mapping input sequences to centroid sequences, and a dict
    of MD5 checksums of any sequences flagged as chimeras.
    """
    check_tools(["usearch"], debug)

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
        # Sorting sequences by abundance, largest first
        # If gamma known, can pre-filter to drop trace level entries
        for a, seq in sorted(
            (
                (a, seq)
                for (seq, a) in counts.items()
                if not unoise_gamma or a >= unoise_gamma
            ),
            key=lambda x: (-x[0], x[1]),
        ):
            md5 = md5seq(seq)
            handle.write(f">{md5};size={a}\n{seq}\n")
            md5_to_seq[md5] = seq

    cmd = [
        "usearch",
        "-unoise3",
        input_fasta,
        # "-zotus",
        # output_fasta,
        "-tabbedout",
        output_tsv,
    ]
    if unoise_alpha:
        cmd += ["-unoise_alpha", str(unoise_alpha)]
    if unoise_gamma:
        cmd += ["-minsize", str(unoise_gamma)]
    if abundance_based:
        pass
    if cpu:
        # Currently triggers:
        # WARNING: Option -threads not used
        cmd += ["--threads", str(cpu)]
    run(cmd, debug=debug)
    corrections = {}
    chimeras = {}
    amp_md5 = {}  # maps USEARCH amplicon (centroid) names to MD5
    with open(output_tsv) as handle:
        for line in handle:
            if not line:
                continue
            parts = line.split("\t")
            md5 = parts[0].split(";", 1)[0]  # strip the size information
            seq = md5_to_seq[md5]
            if len(parts) not in (3, 4):
                sys.exit(
                    f"ERROR: Found {len(parts)} fields in "
                    "usearch -tabbedout output, not 3 or 4"
                )
            if parts[1] == "denoise":
                if parts[2].startswith("amp") and len(parts) == 3:
                    # centroid
                    corrections[seq] = seq
                    amp_md5[parts[2].strip("\n").lower()] = md5
                elif parts[2] in ("bad", "shifted") and ";top=" in parts[3]:
                    # hit, mapped to a centroid
                    centroid_md5 = parts[3].split(";top=", 1)[1].split(";", 1)[0]
                    corrections[seq] = md5_to_seq[centroid_md5]
                else:
                    sys.exit(f"ERROR: Unexpected usearch denoise line: {line!r}")
            elif parts[1] == "chfilter":
                if parts[2] == "zotu\n":
                    pass
                elif parts[2] == "chimera":
                    # Expect something like this in final field:
                    # dqm=0;dqt=33;div=18.2;top=(R);parentL=Amp4;parentR=Amp13;
                    # dqm=0;dqt=7;div=3.4;top=(L);parentL=Amp8;parentR=Amp4;
                    # dqm=0;dqt=9;div=4.5;top=Amp54;parentL=Amp25;parentR=Amp13;
                    # Note Amp with upper case A, verus lower case above.
                    # Note not clear to me what top field means...
                    fields = dict(
                        _.split("=", 1)
                        for _ in parts[3].strip("\n").strip(";").lower().split(";")
                    )
                    chimeras[md5] = (
                        amp_md5[fields["parentl"]] + "/" + amp_md5[fields["parentr"]]
                    )
                else:
                    sys.exit(f"ERROR: Unexpected usearch chfilter line: {line!r}")
            else:
                sys.exit(f"ERROR: Unexpected usearch tabbedout line: {line!r}")
    del md5_to_seq
    return corrections, chimeras


def vsearch(
    counts,
    unoise_alpha=None,
    unoise_gamma=None,
    abundance_based=False,
    tmp_dir=None,
    debug=False,
    cpu=0,
):
    """Invoke VSEARCH to run its reimplementation of the UNOISE3 algorithm.

    Argument counts is an (unsorted) dict of sequences (for the same amplicon
    marker) as keys, with their total abundance counts as values.

    Returns a dict mapping input sequences to centroid sequences, and a dict
    of MD5 checksums of any sequences flagged as chimeras.
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
    denoised_fasta = os.path.join(shared_tmp, "denoised.fasta")
    chimera_tsv = os.path.join(shared_tmp, "chimeras.tsv")
    md5_to_seq = {}
    with open(input_fasta, "w") as handle:
        # Output using MD5 with USEARCH naming: md5;size=abundance
        # Sorting sequences by abundance, largest first
        # If we can, pre-filter using gamma
        for a, seq in sorted(
            (
                (a, seq)
                for (seq, a) in counts.items()
                if not unoise_gamma or a >= unoise_gamma
            ),
            key=lambda x: (-x[0], x[1]),
        ):
            md5 = md5seq(seq)
            handle.write(f">{md5};size={a}\n{seq}\n")
            md5_to_seq[md5] = seq

    cmd = [
        "vsearch",
        "--sizein",
        "--sizeout",
        "--cluster_unoise",
        input_fasta,
        # "--centroids",
        # output_fasta,
        "--uc",
        output_tsv,
    ]
    if unoise_alpha:
        cmd += ["-unoise_alpha", str(unoise_alpha)]
    if unoise_gamma:
        cmd += ["-minsize", str(unoise_gamma)]
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

    with open(denoised_fasta, "w") as handle:
        # TODO: Wasteful as will re-compute the post-correction counts
        post_counts = defaultdict(int)
        for seq, a in counts.items():
            if seq not in corrections:
                # Ignored as per UNOISE algorithm
                if unoise_gamma:
                    assert (
                        a < unoise_gamma
                    ), f"{md5seq(seq)} total {a} vs {unoise_gamma}"
                continue
            seq = corrections[seq]
            post_counts[seq] += a
        for seq, a in sorted(post_counts.items(), key=lambda x: (-x[1], x[0])):
            if unoise_gamma:
                assert a >= unoise_gamma
            assert corrections[seq] == seq
            handle.write(f">{md5seq(seq)};size={a}\n{seq}\n")
        del post_counts
    cmd = [
        "vsearch",
        "--sizein",
        "--sizeout",
        "--uchime3_denovo",
        denoised_fasta,
        "--uchimeout",
        chimera_tsv,
    ]
    if cpu:
        cmd += ["--threads", str(cpu)]
    run(cmd, debug=debug)
    chimeras = {}
    with open(chimera_tsv) as handle:
        for line in handle:
            if not line or line.endswith("N\n"):
                continue
            parts = line.split("\t", 5)
            md5 = parts[1].split(";", 1)[0]  # strip the size information
            assert len(md5) == 32, line
            chimeras[md5] = parts[2].split(";", 1)[0] + "/" + parts[3].split(";", 1)[0]

    return corrections, chimeras


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

    Argument algorithm is a string, "unoise-l" for our reimplementation of the
    UNOISE2 algorithm, or "usearch" or "vsearch" to invoke those tools at the
    command line.

    Argument counts is an (unsorted) dict of sequences (for the same amplicon
    marker) as keys, with their total abundance counts as values.

    Returns a dict mapping input sequences to centroid sequences, and dict of
    any chimeras detected (empty for some algorithms).
    """
    start = time()
    if algorithm == "unoise-l":
        # Does not need tmp_dir, cpu
        answer = unoise(
            counts, unoise_alpha, unoise_gamma, abundance_based, debug=False
        )
    elif algorithm == "usearch":
        # Does not need cpu?
        answer = usearch(
            counts, unoise_alpha, unoise_gamma, abundance_based, tmp_dir, debug, cpu
        )
    elif algorithm == "vsearch":
        answer = vsearch(
            counts, unoise_alpha, unoise_gamma, abundance_based, tmp_dir, debug, cpu
        )
    elif algorithm == "-":
        sys.exit("ERROR: Internal function denoise_algorithm called without algorithm.")
    else:
        sys.exit(
            f"ERROR: denoise_algorithm called with {algorithm} (unknown algorithm)."
        )
    time_corrections = time() - start
    sys.stderr.write(
        f"Spent {time_corrections:0.1f}s running {algorithm} for read-corrections\n"
    )
    return answer


def main(
    inputs,
    output,
    denoise_algorithm,
    total_min_abundance=0,
    min_length=0,
    max_length=sys.maxsize,
    unoise_alpha=None,  # e.g. 2.0,
    unoise_gamma=None,  # e.g. 4,
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
    corrections, chimeras = read_correction(
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
        if seq not in corrections:
            # Ignored as per UNOISE algorithm
            if unoise_gamma:
                assert (
                    totals[seq] < unoise_gamma
                ), f"{md5seq(seq)} total {totals[seq]} vs {unoise_gamma}"
            continue
        seq = corrections[seq]
        new_totals[seq] += a
    sys.stderr.write(
        f"{denoise_algorithm.upper()} reduced unique ASVs from "
        f"{len(totals)} to {len(new_totals)}, "
        f"max abundance now {max(new_totals.values(), default=0)}\n"
    )
    if chimeras:
        sys.stderr.write(
            f"{denoise_algorithm.upper()} flagged {len(chimeras)} as chimeras\n"
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
            msg = "Does not support gzipped output to stdout."
            raise ValueError(msg)
        out_handle = sys.stdout
    elif gzipped:
        out_handle = gzip.open(output, "wt")
    else:
        out_handle = open(output, "w")

    for count, seq in values:
        md5 = md5seq(seq)
        if md5 in chimeras:
            # Write any dict value as it is...
            out_handle.write(f">{md5}_{count} chimera {chimeras[md5]}\n{seq}\n")
        else:
            out_handle.write(f">{md5}_{count}\n{seq}\n")
    if output != "-":
        out_handle.close()
