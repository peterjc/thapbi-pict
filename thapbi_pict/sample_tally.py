# Copyright 2022-2024 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Prepare a non-redundant TSV file using MD5 naming.

This implements the ``thapbi_pict sample-tally ...`` command.
"""

from __future__ import annotations

import gzip
import os
import shutil
import sys
import tempfile
from collections import Counter
from collections import defaultdict
from math import ceil
from time import time
from typing import Optional
from typing import Union

from Bio.SeqIO.FastaIO import SimpleFastaParser

from .denoise import read_correction
from .prepare import load_marker_defs
from .utils import abundance_from_read_name
from .utils import export_sample_biom
from .utils import file_to_sample_name
from .utils import is_spike_in
from .utils import load_fasta_header
from .utils import md5seq


def main(
    inputs: Union[str, list[str]],
    synthetic_controls: list[str],
    negative_controls: list[str],
    output: str,
    session,
    marker: Optional[str] = None,
    spike_genus=None,
    fasta=None,
    min_abundance: int = 100,
    min_abundance_fraction: float = 0.001,
    total_min_abundance: int = 0,
    min_length: int = 0,
    max_length: int = sys.maxsize,
    denoise_algorithm: str = "-",
    unoise_alpha: Optional[float] = None,
    unoise_gamma: Optional[int] = None,
    gzipped: bool = False,  # output
    biom: Optional[str] = None,
    tmp_dir: Optional[str] = None,
    debug: bool = False,
    cpu: int = 0,
) -> None:
    """Implement the ``thapbi_pict sample-tally`` command.

    Arguments min_length and max_length are applied while loading the input
    per-sample FASTA files.

    Argument algorithm is a string, "-" for no read correction (denoising),
    "unoise-l" for our reimplementation of the UNOISE2 algorithm, or "usearch"
    or "vsearch" to invoke those tools at the command line.

    Arguments min_abundance and min_abundance_fraction are applied per-sample
    (after denoising if being used), increased by pool if negative or
    synthetic controls are given respectively. Comma separated string argument
    spike_genus is treated case insensitively.
    """
    if isinstance(inputs, str):
        inputs = [inputs]
    assert isinstance(inputs, list)
    assert inputs

    assert "-" not in inputs
    if synthetic_controls:
        # Ignore anything not also in the inputs
        synthetic_controls = [
            file_to_sample_name(_) for _ in synthetic_controls if _ in inputs
        ]
    else:
        synthetic_controls = []
    if negative_controls:
        # Ignore anything not also in the inputs
        negative_controls = [
            file_to_sample_name(_) for _ in negative_controls if _ in inputs
        ]
    else:
        negative_controls = []
    controls = set(negative_controls + synthetic_controls)
    if debug and not controls:
        sys.stderr.write("DEBUG: No control samples\n")

    if os.path.isdir(output):
        sys.exit("ERROR: Output directory given, want a filename.")

    marker_definitions = load_marker_defs(session, spike_genus)
    spikes: list[tuple[str, str, set[str]]] = []
    if marker:
        if not marker_definitions:
            sys.exit("ERROR: Marker given with -k / --marker not in DB.")
        marker_def = marker_definitions[marker]
    elif len(marker_definitions) > 1:
        sys.exit("ERROR: DB has multiple markers, please set -k / --marker.")
    else:
        # Implicit -k / --marker choice as only one in the DB
        marker, marker_def = marker_definitions.popitem()
    assert isinstance(marker_def["spike_kmers"], list)
    spikes = marker_def["spike_kmers"]
    del marker_def
    del marker_definitions
    assert marker

    totals: dict[str, int] = Counter()
    counts: dict[tuple[str, str], int] = Counter()
    sample_cutadapt: dict[str, int] = {}  # before any thresholds
    samples: Union[set[str], list[str]] = set()
    sample_pool: dict[str, str] = {}
    sample_headers: dict[str, dict] = {}
    for filename in inputs:
        # Assuming FASTA for now
        if debug:
            sys.stderr.write(f"DEBUG: Parsing {filename}\n")
        sample = file_to_sample_name(filename)
        assert sample not in samples, f"ERROR: Duplicate stem from {filename}"
        assert isinstance(samples, set)  # will later turn it into a list
        samples.add(sample)
        sample_headers[sample] = load_fasta_header(filename)
        if not sample_headers[sample]:
            if min_abundance_fraction or synthetic_controls:
                sys.exit(
                    f"ERROR: Missing FASTA header in {filename}, "
                    "required for fractional abundance threshold\n"
                )
            sys.stderr.write(f"WARNING: Missing FASTA header in {filename}\n")
            sample_cutadapt[sample] = 0
            # Will assume marker matches, and use default for threshold_pool
        else:
            if marker != sample_headers[sample]["marker"]:
                sys.exit(
                    f"ERROR: Expected marker {marker},"
                    f" not {sample_headers[sample]['marker']} from {filename}"
                )
            assert "raw_fastq" in sample_headers[sample], sample_headers[sample]
            sample_cutadapt[sample] = int(sample_headers[sample]["cutadapt"])
        sample_pool[sample] = sample_headers[sample].get("threshold_pool", "default")
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

    chimeras: dict[str, str] = {}
    if denoise_algorithm != "-":
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
        new_counts: dict[tuple[str, str], int] = defaultdict(int)
        new_totals: dict[str, int] = defaultdict(int)
        for (seq, sample), a in counts.items():
            if seq not in corrections:
                # Should have been ignored as per UNOISE algorithm
                if unoise_gamma:
                    assert (
                        totals[seq] < unoise_gamma
                    ), f"{md5seq(seq)} total {totals[seq]} vs {unoise_gamma}"
                continue
            seq = corrections[seq]
            new_counts[seq, sample] += a
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
        counts = new_counts
        totals = new_totals
        del new_totals, new_counts, corrections

    # Drop entries so low in abundance that they'll never be accepted, nor
    # increase the pool threshold. Drop before bothering to check if spikes.
    # Note factor of 0.5 due to heuristic below for threshold adjustment
    # (things slightly below the absolute threshold might still contribut to
    # raising the max non-spike in counts and this the automatic fractional
    # threshold)
    ultra_low = ceil(0.5 * min_abundance)
    if ultra_low > 1:
        if debug:
            sys.stderr.write("DEBUG: Dropping ultra low abundance entries\n")
        new_counts = defaultdict(int)
        new_totals = defaultdict(int)
        before = len(totals)
        del totals
        while counts:
            (seq, sample), a = counts.popitem()
            if a >= ultra_low:
                new_counts[seq, sample] = a
                new_totals[seq] += a
        if debug or len(new_totals) < before:
            sys.stderr.write(
                "Excluding ultra-low abundance entries reduced unique ASVs "
                f"from {before} to {len(new_totals)}.\n"
            )
        del before
        counts = new_counts
        totals = new_totals
        del new_totals, new_counts

    pool_absolute_threshold: dict[str, int] = {}
    pool_fraction_threshold: dict[str, float] = {}
    max_spike_abundance = {sample: 0 for sample in samples}  # exporting in metadata
    max_non_spike_abundance = {sample: 0 for sample in samples}  # exporting in metadata
    if spikes:
        if debug:
            sys.stderr.write(
                "DEBUG: About to identify spike-in synthetic sequences...\n"
            )
        start = time()
        for seq in totals:
            # Calling is_spike_in is relatively expensive, but will be of less
            # interest on the tail end low abundance samples.
            # Should we sort totals by count?
            assert totals[seq] >= ultra_low, seq
            if totals[seq] < min(max_spike_abundance.values()) and totals[seq] < min(
                max_non_spike_abundance.values()
            ):
                # Note counts[seq, sample] <= totals[seq] so will not be able to exceed
                continue
            if any(
                min(max_spike_abundance[sample], max_non_spike_abundance[sample])
                < counts[seq, sample]
                for sample in samples
            ):
                # This could raise the max (non-)spike-in abundance
                if is_spike_in(seq, spikes):
                    for sample in samples:
                        if max_spike_abundance[sample] < counts[seq, sample]:
                            max_spike_abundance[sample] = counts[seq, sample]
                else:
                    for sample in samples:
                        if max_non_spike_abundance[sample] < counts[seq, sample]:
                            max_non_spike_abundance[sample] = counts[seq, sample]
        time_spike_tagging = time() - start
        if debug:
            sys.stderr.write(
                f"DEBUG: Spent {time_spike_tagging:0.1f}s tagging spike-in sequences.\n"
            )
    else:
        for (_, sample), a in counts.items():
            max_non_spike_abundance[sample] = max(max_non_spike_abundance[sample], a)

    start = time()
    sample_threshold: dict[str, int] = {}
    if controls:
        if debug:
            sys.stderr.write("DEBUG: Applying dynamic abundance thresholds...\n")
        for sample in controls:
            pool = sample_pool[sample]
            a = max_non_spike_abundance[sample]
            # Ignore rest of pool. Respect FASTA header which could be more.
            # Note heuristic half thresholds inherited from prepare-reads
            sample_threshold[sample] = threshold = max(
                int(sample_headers[sample].get("threshold", 0)),
                # use half for a synthetic (fractional) control:
                (
                    ceil(min_abundance * 0.5)
                    if sample in synthetic_controls
                    else min_abundance
                ),
                # use half for a negative (absolute) control:
                ceil(
                    sample_cutadapt[sample]
                    * min_abundance_fraction
                    * (0.5 if sample in negative_controls else 1.0)
                ),
            )
            if a < threshold:
                if debug:
                    sys.stderr.write(
                        f"DEBUG: Control {sample} (in pool {pool}) does not "
                        f"exceed its threshold {threshold}\n"
                    )
                continue
            if debug:
                sys.stderr.write(
                    f"DEBUG: Control {sample} in (pool {pool}) gets "
                    f"threshold {threshold}, exceeded with {a}\n"
                )

            # Does this control raise the pool thresholds?
            if sample in negative_controls and min_abundance < a:
                pool_absolute_threshold[pool] = max(
                    a, pool_absolute_threshold.get(pool, min_abundance)
                )
                if debug:
                    sys.stderr.write(
                        f"DEBUG: Negative control {sample} says increase {pool} "
                        f"absolute abundance threshold from {min_abundance} to {a}\n"
                    )
            if (
                sample in synthetic_controls
                and min_abundance_fraction * sample_cutadapt[sample] < a
            ):
                f = a / sample_cutadapt[sample]
                pool_fraction_threshold[pool] = max(
                    f, pool_fraction_threshold.get(pool, min_abundance_fraction)
                )
                if f > 0.5:
                    sys.exit(
                        f"ERROR: Control {sample} suggests extremely high"
                        f" fractional abundance threshold {f*100:.1f}%\n"
                    )
                elif f > 0.05:
                    sys.stderr.write(
                        f"WARNING: Control {sample} suggests overly high"
                        f" fractional abundance threshold {f*100:.1f}%\n"
                    )
                elif debug:
                    sys.stderr.write(
                        f"DEBUG: Synthetic control {sample} says increase "
                        f"{pool} fractional abundance threshold from "
                        f"{min_abundance_fraction*100:.4f} to {f*100:.4f}%\n"
                    )

    for pool, value in pool_absolute_threshold.items():
        if value > min_abundance:
            sys.stderr.write(
                f"Negative controls for marker {marker} increased "
                f"pool {pool} absolute abundance threshold to {value}.\n"
            )
    for pool, fraction in pool_fraction_threshold.items():
        if fraction > min_abundance_fraction:
            sys.stderr.write(
                f"Negative controls for marker {marker} increased "
                f"pool {pool} fractional abundance threshold to "
                f"{fraction*100:0.4f}%\n."
            )

    # Apply any dynamic abundance threshold increases from controls:
    for sample in samples:
        if sample in controls:
            # Done above
            assert sample in sample_threshold, sample
            continue
        # Apply dynamic pool threshold from controls
        pool = sample_pool[sample]
        sample_threshold[sample] = threshold = max(
            int(sample_headers[sample].get("threshold", 0)),
            pool_absolute_threshold.get(pool, min_abundance),
            ceil(
                sample_cutadapt[sample]
                * pool_fraction_threshold.get(pool, min_abundance_fraction)
            ),
        )
        if debug and not controls:
            sys.stderr.write(
                f"DEBUG: {sample} in pool {pool} gets threshold {threshold}\n"
            )
    del pool_absolute_threshold
    del pool_fraction_threshold

    new_counts = defaultdict(int)
    new_totals = defaultdict(int)
    before = len(totals)
    del totals
    while counts:
        (seq, sample), a = counts.popitem()
        if a >= sample_threshold[sample]:
            new_counts[seq, sample] = a
            new_totals[seq] += a
    sys.stderr.write(
        f"Sample pool abundance thresholds reduced unique ASVs from {before} to "
        f"{len(new_totals)}.\n"
    )
    counts = new_counts
    totals = new_totals
    del new_totals, new_counts

    if total_min_abundance:
        new_counts = defaultdict(int)
        new_totals = defaultdict(int)
        while counts:
            (seq, sample), a = counts.popitem()
            if totals[seq] >= total_min_abundance:
                new_counts[seq, sample] = a
                new_totals[seq] += a
        if debug or len(new_totals) < len(totals):
            sys.stderr.write(
                f"Total abundance threshold {total_min_abundance} reduced "
                f"unique ASVs from {len(totals)} to {len(new_totals)}.\n"
            )
        counts = new_counts
        totals = new_totals
        del new_totals, new_counts

    if debug:
        time_thresholds = time() - start
        sys.stderr.write(
            f"DEBUG: Spent {time_thresholds:0.1f}s on abundance thresholds\n"
        )

    if chimeras:
        # Filter the chimeras list to drop low abundance entries
        # Would be better to call VSEARCH chimera detection here?
        chimeras = {
            md5seq(seq): chimeras[md5seq(seq)]
            for seq in totals
            if md5seq(seq) in chimeras
        }
        sys.stderr.write(
            f"Have {len(chimeras)} chimeras passing the abundance thresholds.\n"
        )
    del before

    # Cap max (non-)spike values by sample specific thresholds
    # (to match historic behaviour of pipeline via prepare-reads)
    for sample in samples:
        if max_spike_abundance[sample] < sample_threshold[sample]:
            max_spike_abundance[sample] = 0
        if max_non_spike_abundance[sample] < sample_threshold[sample]:
            max_non_spike_abundance[sample] = 0

    samples = sorted(samples)  # converting from set to list
    count_seq = sorted(
        ((count, seq) for seq, count in totals.items()),
        # (sort by marker), then put the highest abundance entries first:
        key=lambda x: (-x[0], x[1]),
    )
    del totals

    if debug:
        sys.stderr.write("DEBUG: Sorted, about to start output...\n")

    # TODO - avoid double definition here and in summary code:
    stats_fields: tuple[str, ...] = (
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
    if not spike_genus:
        stats_fields = tuple(
            stat
            for stat in stats_fields
            if stat not in ("Max non-spike", "Max spike-in")
        )
    sample_stats: dict[str, dict[str, Union[int, str, None]]] = {
        sample: {} for sample in samples
    }
    for stat in stats_fields:
        if stat == "Max non-spike":
            if not spike_genus:
                continue
            for sample in samples:
                sample_stats[sample][stat] = max_non_spike_abundance[sample]
        elif stat == "Max spike-in":
            if not spike_genus:
                continue
            for sample in samples:
                sample_stats[sample][stat] = max_spike_abundance[sample]
        elif stat == "Control":
            terms = {
                (True, True): "Neg/Synth",
                (True, False): "Negative",
                (False, True): "Synthetic",
                (False, False): "Sample",
            }
            for sample in samples:
                sample_stats[sample][stat] = terms[
                    sample in negative_controls, sample in synthetic_controls
                ]
            del terms
        elif stat == "Threshold":
            for sample in samples:
                sample_stats[sample][stat] = sample_threshold[sample]
        elif stat == "Threshold pool":
            # Get from FASTA headers
            stat_values = [
                sample_headers[sample].get("threshold_pool", None) for sample in samples
            ]
            # Try to remove any common folder prefix like raw_data/
            # or C:/Users/... or /tmp/... - note making str type explicit:
            common = os.path.commonpath([str(_) for _ in stat_values])
            if set(stat_values) == {None}:
                stat_values = [None] * len(stat_values)
            elif len(set(stat_values)) > 1 and common:
                if debug:
                    sys.stderr.write(
                        f"DEBUG: Dropping threshold pool common prefix {common}\n"
                    )
                stat_values = [_[len(common) + 1 :] for _ in stat_values]
            elif len(set(stat_values)) == 1:
                # discard possibly platform specific start
                common = common.rsplit(os.path.sep, 1)[-1]
                stat_values = [common] * len(stat_values)
            for sample, pool in zip(samples, stat_values):
                sample_stats[sample][stat] = pool
            del common, pool, stat_values
        else:
            # Get from FASTA headers
            for sample in samples:
                sample_stats[sample][stat] = sample_headers[sample].get(
                    stat.lower().replace(" ", "_"), None
                )

    start = time()
    # First, BIOM output for the AVS vs samples table
    if biom:
        if debug:
            sys.stderr.write(f"DEBUG: Starting BIOM output {biom}\n")
        if tmp_dir:
            # Using same file names, but in tmp folder:
            tmp_biom = os.path.join(tmp_dir, os.path.basename(biom))
        else:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".biom", delete=False
            ) as handle:
                handle.write("Put BIOM data here please\n")
                tmp_biom = handle.name
        if export_sample_biom(
            tmp_biom,
            # Created expected marker based sequence dict using md5:
            {(marker, md5seq(seq)): seq for total, seq in count_seq},
            {
                (marker, md5seq(seq)): {} for total, seq in count_seq
            },  # no sequence metadata
            sample_stats,  # internal sample metadata from FASTQ processing
            {
                (marker, md5seq(seq), sample): counts.get((seq, sample), 0)
                for sample in samples
                for (total, seq) in count_seq
            },
        ):
            # Move our temp file into position...
            shutil.move(tmp_biom, biom)
            if debug:
                sys.stderr.write(f"DEBUG: Wrote {biom}\n")
        else:
            sys.exit("ERROR: Missing optional Python library for BIOM output")
        if debug:
            time_output = time() - start
            sys.stderr.write(
                f"DEBUG: Spent {time_output:0.1f}s writing BIOM output file\n"
            )

    # Now the main TSV output for our pipeline...
    start = time()
    if output == "-":
        if fasta == "-":
            sys.exit("ERROR: Don't use stdout for both TSV and FASTA output.")
        if gzipped:
            msg = "Does not support gzipped output to stdout."
            raise ValueError(msg)
        out_handle = sys.stdout
    elif gzipped:
        out_handle = gzip.open(output, "wt")
    else:
        out_handle = open(output, "w")
    for stat in stats_fields:
        # Using "-" as missing value to match default in summary reports
        stat_values = [sample_stats[sample].get(stat, None) for sample in samples]
        missing = [_ is None for _ in stat_values]
        if all(missing):
            sys.stderr.write(f"WARNING: Missing all {stat} values in FASTA headers\n")
            continue  # don't export it
        elif any(missing):
            sys.stderr.write(f"WARNING: Missing some {stat} values in FASTA headers\n")
        del missing
        out_handle.write(
            "\t".join(
                ["#" + stat]
                + ["-" if _ is None else str(_) for _ in stat_values]
                + ([""] if chimeras else [])
                + ["\n"]
            )
        )
        del stat_values
    fields = ["Marker/MD5_abundance", *samples, "Sequence"]
    if chimeras:
        fields += ["Chimeras"]
    out_handle.write("#" + "\t".join(fields) + "\n")

    if fasta == "-":
        if gzipped:
            msg = "Does not support gzipped output to stdout."
            raise ValueError(msg)
        fasta_handle = sys.stdout
    elif fasta and gzipped:
        fasta_handle = gzip.open(fasta, "wt")
    elif fasta:
        fasta_handle = open(fasta, "w")
    else:
        fasta_handle = None
    for count, seq in count_seq:
        # This is the last time we'll use the counts[seq, *] values,
        # so use pop to gradually reduce the data held in memory.
        data = "\t".join(str(counts.pop((seq, sample), 0)) for sample in samples)
        md5 = md5seq(seq)
        if md5 in chimeras:
            out_handle.write(
                f"{marker}/{md5}_{count}\t{data}\t{seq}\t{chimeras[md5]}\n"
            )
        elif chimeras:
            # Empty field
            out_handle.write(f"{marker}/{md5}_{count}\t{data}\t{seq}\t\n")
        else:
            # Can omit the column
            out_handle.write(f"{marker}/{md5}_{count}\t{data}\t{seq}\n")
        if fasta_handle:
            # Does not export per-sample counts
            if md5 in chimeras:
                # Write any dict value as it is...
                fasta_handle.write(
                    f">{marker}/{md5}_{count} chimera {chimeras[md5]}\n{seq}\n"
                )
            else:
                fasta_handle.write(f">{marker}/{md5}_{count}\n{seq}\n")
    if output != "-":
        assert out_handle is not None
        out_handle.close()
    if fasta and fasta != "-":
        assert fasta_handle is not None
        fasta_handle.close()

    if debug:
        time_output = time() - start
        sys.stderr.write(f"DEBUG: Spent {time_output:0.1f}s writing TSV output file\n")
