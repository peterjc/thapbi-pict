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
from math import ceil

from Bio.SeqIO.FastaIO import SimpleFastaParser

from .prepare import load_fasta_header
from .prepare import load_marker_defs
from .utils import abundance_from_read_name
from .utils import file_to_sample_name
from .utils import is_spike_in
from .utils import md5seq


def main(
    inputs,
    synthetic_controls,
    negative_controls,
    output,
    session,
    spike_genus,
    fasta=None,
    min_abundance=100,
    min_abundance_fraction=0.001,
    total_min_abundance=0,
    min_length=0,
    max_length=sys.maxsize,
    gzipped=False,  # output
    debug=False,
):
    """Implement the ``thapbi_pict sample-tally`` command.

    Arguments min_length and max_length are applied while loading the input
    per-sample FASTA files.

    Arguments min_abundance and min_abundance_fraction are applied per-sample,
    increased by pool if negative or synthetic controls are given respectively.
    Comma separated string argument spike_genus is treated case insensitively.

    Results are sorted by marker, decreasing abundance then alphabetically by
    sequence.
    """
    if isinstance(inputs, str):
        inputs = [inputs]
    assert isinstance(inputs, list)
    assert inputs

    if synthetic_controls:
        # Possible in a pipeline setting may need to pass a null value,
        # e.g. -n "" or -n "-"
        synthetic_controls = [_ for _ in synthetic_controls if _ and _ != "-"]
    else:
        synthetic_controls = []
    if negative_controls:
        # Possible in a pipeline setting may need to pass a null value,
        # e.g. -n "" or -n "-"
        negative_controls = [_ for _ in negative_controls if _ and _ != "-"]
    else:
        negative_controls = []
    for filename in synthetic_controls + negative_controls:
        if filename not in inputs:
            inputs.append(filename)
    synthetic_controls = [file_to_sample_name(_) for _ in synthetic_controls]
    negative_controls = [file_to_sample_name(_) for _ in negative_controls]
    controls = set(negative_controls + synthetic_controls)
    if debug and not controls:
        sys.stderr.write("DEBUG: No control samples\n")

    if os.path.isdir(output):
        sys.exit("ERROR: Output directory given, want a filename.")

    marker_definitions = load_marker_defs(session, spike_genus)

    totals = Counter()
    counts = defaultdict(int)
    sample_counts = defaultdict(int)
    sample_marker_cutadapt = {}  # before any thresholds
    samples = set()
    sample_pool = {}
    sample_headers = {}
    for filename in inputs:
        # Assuming FASTA for now
        if debug:
            sys.stderr.write(f"DEBUG: Parsing {filename}\n")
        # TODO - Assuming just one marker for now:
        sample = file_to_sample_name(filename)
        assert sample not in samples, f"ERROR: Duplicate stem from {filename}"
        samples.add(sample)
        sample_headers[sample] = load_fasta_header(filename)
        marker = sample_headers[sample]["marker"]
        if marker not in marker_definitions:
            # Probably using the default DB by mistake!
            sys.exit(
                f"ERROR: Marker {marker} not in DB, should be one of "
                + ", ".join(marker_definitions)
            )
        assert "raw_fastq" in sample_headers[sample], sample_headers[sample]
        sample_marker_cutadapt[sample, marker] = int(sample_headers[sample]["cutadapt"])
        sample_pool[sample] = sample_headers[sample].get("threshold_pool", "default")
        with open(filename) as handle:
            for _, seq in SimpleFastaParser(handle):
                seq = seq.upper()
                if min_length <= len(seq) <= max_length:
                    a = abundance_from_read_name(_.split(None, 1)[0])
                    totals[marker, seq] += a
                    counts[marker, seq, sample] += a
                    sample_counts[marker, sample] += a

    if totals:
        sys.stderr.write(
            f"Loaded {len(totals)} unique sequences from {sum(totals.values())}"
            f" in total within length range, max abundance {max(totals.values())}\n"
        )
    else:
        sys.stderr.write("WARNING: Loaded zero sequences within length range\n")

    # TODO - support for multiple markers?
    pool_absolute_threshold = {}
    pool_fraction_threshold = {}
    max_spike_abundance = defaultdict(int)  # exporting in metadata
    max_non_spike_abundance = defaultdict(int)  # exporting in metadata
    for (marker, seq) in totals:
        spikes = marker_definitions[marker]["spike_kmers"]
        if is_spike_in(seq, spikes):
            for sample in samples:
                if counts[marker, seq, sample]:
                    max_spike_abundance[sample] = max(
                        max_spike_abundance[sample], counts[marker, seq, sample]
                    )
        else:
            for sample in samples:
                if counts[marker, seq, sample]:
                    max_non_spike_abundance[sample] = max(
                        max_non_spike_abundance[sample], counts[marker, seq, sample]
                    )
    if controls:
        if debug:
            sys.stderr.write("DEBUG: Applying dynamic abundance thresholds...\n")
        # Samples with zero non-spike-ins can't raise the threshold:
        for sample in max_non_spike_abundance:
            pool = sample_pool[sample]
            a = max_non_spike_abundance[sample]
            if sample in negative_controls and min_abundance < a:
                pool_absolute_threshold[pool] = max(
                    a, pool_absolute_threshold.get(pool, min_abundance)
                )
                if debug:
                    sys.stderr.write(
                        f"Negative control {sample} says increase {pool} "
                        f"absolute abundance threshold from {min_abundance} to {a}\n"
                    )
            elif debug and sample in negative_controls:
                sys.stderr.write(
                    f"Negative control {sample} suggests drop {pool} absolute "
                    f"abundance threshold from {min_abundance} to {a} (lower)\n"
                )
            if (
                sample in synthetic_controls
                and min_abundance_fraction * sample_marker_cutadapt[sample, marker] < a
            ):
                f = a / sample_marker_cutadapt[sample, marker]
                pool_fraction_threshold[pool] = max(
                    f, pool_fraction_threshold.get(pool, min_abundance_fraction)
                )
                if debug:
                    sys.stderr.write(
                        f"Synthetic control {sample} says increase {pool} "
                        "fractional abundance threshold from "
                        f"{min_abundance_fraction*100:.4f} to {f*100:.4f}%\n"
                    )
            elif debug and sample in synthetic_controls:
                f = a / sample_marker_cutadapt[sample, marker]
                sys.stderr.write(
                    f"Synthetic control {sample} suggests drop {pool} fractional "
                    f"abundance threshold from {min_abundance_fraction*100:.4f} "
                    f"to {f*100:.4f}% (lower)\n"
                )

    # Check for any dynamic abundance threshold increases from controls:
    sample_threshold = {}
    for sample in samples:
        if sample in controls:
            # Not dynamic (respect FASTA header which could be more)
            # Note complicated logic inherited from prepare-reads
            threshold = max(
                int(sample_headers[sample].get("threshold", 0)),
                # use half for a synthetic (fractional) control:
                ceil(min_abundance * 0.5)
                if sample in synthetic_controls
                else min_abundance,
                # use half for a negative (absolute) control:
                ceil(
                    sample_counts[marker, sample]
                    * min_abundance_fraction
                    * (0.5 if sample in negative_controls else 1.0)
                ),
            )
            if debug:
                sys.stderr.write(
                    f"DEBUG: Control {sample} in pool {pool} gets "
                    f"threshold {threshold}\n"
                )
        else:
            # Apply dynamic pool threshold from controls
            pool = sample_pool[sample]
            threshold = max(
                int(sample_headers[sample].get("threshold", 0)),
                pool_absolute_threshold.get(pool, min_abundance),
                ceil(
                    sample_counts[marker, sample]
                    * pool_fraction_threshold.get(pool, min_abundance_fraction)
                ),
            )
            if debug:
                sys.stderr.write(
                    f"DEBUG: {sample} in pool {pool} gets threshold {threshold}\n"
                )
        sample_threshold[sample] = threshold

    new_counts = defaultdict(int)
    new_totals = defaultdict(int)
    for (marker, seq, sample), a in counts.items():
        if a >= sample_threshold[sample]:
            new_counts[marker, seq, sample] = a
            new_totals[marker, seq] += a
    sys.stderr.write(
        f"Abundance thresholds reduced unique ASVs from {len(totals)} to "
        f"{len(new_totals)}.\n"
    )
    counts = new_counts
    totals = new_totals
    del new_totals, new_counts

    if total_min_abundance:
        new_counts = defaultdict(int)
        new_totals = defaultdict(int)
        for (marker, seq, sample), a in counts.items():
            if totals[marker, seq] >= total_min_abundance:
                new_counts[marker, seq, sample] = a
                new_totals[marker, seq] += a
        sys.stderr.write(
            f"Total abundance threshold {total_min_abundance} reduced unique "
            f"ASVs from {len(totals)} to {len(new_totals)}.\n"
        )
        counts = new_counts
        totals = new_totals
        del new_totals, new_counts

    samples = sorted(samples)
    values = sorted(
        ((marker, count, seq) for (marker, seq), count in totals.items()),
        # sort by marker, then put the highest abundance entries first:
        key=lambda x: (x[0], -x[1], x[2]),
    )
    del totals

    # TODO - avoid double definition here and in summary code:
    stats_fields = (
        "Raw FASTQ",
        "Flash",
        "Cutadapt",
        "Threshold pool",
        "Threshold",
        "Control",
        "Max non-spike",
        "Max spike-in",
        "Singletons",
    )

    if output == "-":
        if fasta == "-":
            sys.exit("ERROR: Don't use stdout for both TSV and FASTA output.")
        if gzipped:
            raise ValueError("Does not support gzipped output to stdout.")
        out_handle = sys.stdout
    elif gzipped:
        out_handle = gzip.open(output, "wt")
    else:
        out_handle = open(output, "w")
    missing = [None] * len(samples)
    for stat in stats_fields:
        if stat == "Max non-spike":
            if not spike_genus:
                continue
            stat_values = [max_non_spike_abundance[sample] for sample in samples]
        elif stat == "Max spike-in":
            if not spike_genus:
                continue
            stat_values = [max_spike_abundance[sample] for sample in samples]
        elif stat == "Control":
            terms = {
                (True, True): "Neg/Synth",
                (True, False): "Negative",
                (False, True): "Synthetic",
                (False, False): "Sample",
            }
            stat_values = [
                terms[sample in negative_controls, sample in synthetic_controls]
                for sample in samples
            ]
            del terms
        elif stat == "Threshold":
            stat_values = [sample_threshold[sample] for sample in samples]
        else:
            # Get from FASTA headers
            stat_values = [
                sample_headers[sample].get(stat.lower().replace(" ", "_"), None)
                for sample in samples
            ]
            if stat_values == missing:
                sys.stderr.write(
                    f"WARNING: Missing all {stat} values in FASTA headers\n"
                )
                continue
            if None in stat_values:
                sys.stderr.write(
                    f"WARNING: Missing some {stat} values in FASTA headers\n"
                )
        if stat == "Threshold pool":
            # Try to remove any common folder prefix like raw_data/
            # or C:/Users/... or /tmp/...
            common = os.path.commonpath(set(stat_values))
            if len(set(stat_values)) > 1 and common:
                if debug:
                    sys.stderr.write(
                        f"DEBUG: Dropping threshold pool common prefix {common}\n"
                    )
                stat_values = [_[len(common) + 1 :] for _ in stat_values]
        # Using "-" as missing value to match default in summary reports
        out_handle.write(
            "\t".join(
                ["#" + stat]
                + ["-" if _ is None else str(_) for _ in stat_values]
                + ["\n"]
            )
        )
        del stat_values
    out_handle.write("\t".join(["#Marker/MD5_abundance"] + samples + ["Sequence\n"]))

    if fasta == "-":
        if gzipped:
            raise ValueError("Does not support gzipped output to stdout.")
        fasta_handle = sys.stdout
    elif fasta and gzipped:
        fasta_handle = gzip.open(fasta, "wt")
    elif fasta:
        fasta_handle = open(fasta, "w")
    else:
        fasta_handle = None
    for marker, count, seq in values:
        data = "\t".join(str(counts[marker, seq, sample]) for sample in samples)
        md5 = md5seq(seq)
        out_handle.write(f"{marker}/{md5}_{count}\t{data}\t{seq}\n")
        if fasta_handle:
            # Does not export per-sample counts
            # TODO - Include the marker? Older fasta-nr command did not.
            fasta_handle.write(f">{md5}_{count}\n{seq}\n")
    if output != "-":
        out_handle.close()
    if fasta and fasta != "-":
        fasta_handle.close()
