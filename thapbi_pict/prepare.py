# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Prepare raw amplicon sequencing reads (trimming, merging, etc).

This implements the ``thapbi_pict prepare-reads ...`` command.
"""
import gzip
import os
import shutil
import sys
import tempfile
from collections import Counter
from math import ceil
from time import time

from Bio.Seq import reverse_complement
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from sqlalchemy.orm import aliased
from sqlalchemy.orm import contains_eager

from .db_orm import MarkerDef
from .db_orm import MarkerSeq
from .db_orm import SeqSource
from .db_orm import Taxonomy
from .utils import abundance_from_read_name
from .utils import abundance_values_in_fasta
from .utils import load_fasta_header
from .utils import md5seq
from .utils import primer_clean
from .utils import run
from .versions import check_tools

KMER_LENGTH = 31


def find_fastq_pairs(
    filenames_or_folders, ext=(".fastq", ".fastq.gz"), ignore_prefixes=None, debug=False
):
    """Interpret a list of filenames and/or foldernames.

    Returns a list of tuples (stem, left filename, right filename)
    where stem is intended for use in logging and output naming,
    and may include a directory name.

    The filenames will be normalised relative to the current directory
    (so that we can directly compare file lists which could have been
    defined inconsistently by the user).
    """
    if ignore_prefixes:
        assert isinstance(ignore_prefixes, tuple), ignore_prefixes
    if debug and ignore_prefixes:
        sys.stderr.write(
            "DEBUG: Finding FASTQ pairs ignoring prefixes: %s\n"
            % ", ".join(repr(_) for _ in ignore_prefixes)
        )

    answer = []
    for x in filenames_or_folders:
        if not os.path.isabs(x):
            x = os.path.relpath(x)
        x = os.path.normpath(x)
        if os.path.isdir(x):
            if debug:
                sys.stderr.write(f"Walking directory {x!r}\n")
            for root, _, files in os.walk(x, followlinks=True):
                for f in files:
                    if f.endswith(ext):
                        if ignore_prefixes and f.startswith(ignore_prefixes):
                            if debug:
                                sys.stderr.write(
                                    f"WARNING: Ignoring {os.path.join(root, f)}"
                                    " due to prefix\n"
                                )
                            continue
                        answer.append(os.path.join(root, f))
        elif os.path.isfile(x):
            if ignore_prefixes and os.path.split(x)[1].startswith(ignore_prefixes):
                sys.stderr.write(f"WARNING: Ignoring {x} due to prefix\n")
                continue
            answer.append(x)
        else:
            sys.exit(f"ERROR: {x!r} is not a file or a directory\n")
    # Warn if there were duplicates?
    answer = sorted(set(answer))
    if len(answer) % 2:
        sys.exit(f"ERROR: Found {len(answer)} FASTQ files, expected pairs\n")

    pairs = []
    while answer:
        left = answer.pop(0)
        right = answer.pop(0)
        if left.endswith(tuple("_I1_001" + _ for _ in ext)) and right.endswith(
            tuple("_I2_001" + _ for _ in ext)
        ):
            if debug:
                sys.stderr.write(
                    f"WARNING: Ignoring {left!r} and {right!r} due to suffix\n"
                )
            continue
        stem = None
        for suffix_left, suffix_right in zip(
            ("_R1_001", "_R1", "_1"), ("_R2_001", "_R2", "_2")
        ):
            if left.endswith(tuple(suffix_left + _ for _ in ext)) and right.endswith(
                tuple(suffix_right + _ for _ in ext)
            ):
                stem = left.rsplit(suffix_left, 1)[0]
                if stem != right.rsplit(suffix_right, 1)[0]:
                    sys.exit(
                        f"ERROR: Did not recognise {left!r} and {right!r} as a pair\n"
                    )
        if not os.path.isfile(left) and os.path.islink(left):
            sys.stderr.write(f"WARNING: Ignoring {left} as symlink broken\n")
            continue
        if not os.path.isfile(right) and os.path.islink(right):
            sys.stderr.write(f"WARNING: Ignoring {right} as symlink broken\n")
            continue
        for filename in (left, right):
            if not os.path.isfile(filename):
                if os.path.islink(filename):
                    sys.exit(f"ERROR: Broken symlink: {filename}")
                else:
                    # What might cause this - other than deletion after our dir listing?
                    sys.exit(f"ERROR: Missing file: {filename}")
            if not os.stat(filename).st_size:
                if filename.endswith(".gz"):
                    sys.exit(f"ERROR: Empty gzip file {filename}")
                else:
                    sys.stderr.write(f"WARNING: Empty file {filename}\n")
        if not stem:
            sys.exit(
                f"ERROR: Did not recognise pair naming for {left!r} and {right!r}\n"
            )
        pairs.append((stem, left, right))

    return pairs


def parse_cutadapt_stdout(stdout):
    r"""Extract FASTA count before and after cutadapt.

    >>> parse_cutadapt_stdout("...\nTotal reads processed: 5,869\n...\nReads written (passing filters): 5,861 (99.9%)\n...")
    (5869, 5861)
    """  # noqa: E501
    before = None
    after = None
    for line in stdout.strip().splitlines():
        words = line.strip().split()
        if words == ["No", "reads", "processed!"]:
            before = after = 0
        elif words[:-1] == ["Total", "reads", "processed:"]:
            before = int(words[3].replace(",", ""))
        elif words[:4] == [
            "Reads",
            "written",
            "(passing",
            "filters):",
        ]:
            after = int(words[4].replace(",", ""))
    if before is None or after is None:
        sys.exit(
            f"ERROR: Could not extract cutadapt before and after pair count:"
            f"\n\n{stdout}\n"
        )
    return before, after


# Example where not recording bad primers, want smaller number
assert (
    parse_cutadapt_stdout(
        """\
=== Summary ===

Total reads processed:                 106,089
Reads with adapters:                     1,511 (1.4%)
Reads that were too short:               3,271 (3.1%)
Reads that were too long:                    0 (0.0%)
Reads written (passing filters):         1,471 (1.4%)
""",
    )
    == (106089, 1471)
)


def run_cutadapt(
    long_in,
    out_template,
    marker_definitions,
    flip=False,
    debug=False,
    cpu=0,
):
    """Run cutadapt on a single file (i.e. after merging paired FASTQ).

    The input and/or output files may be compressed as long as they
    have an appropriate suffix (e.g. gzipped with ``.gz`` suffix).

    Returns FASTA count before and after cutadapt.
    """
    assert "{name}" in out_template
    # Currently at least, cannot set these at the adapter level...
    min_len = min(_["min_length"] for _ in marker_definitions.values())
    max_len = max(_["max_length"] for _ in marker_definitions.values())

    cmd = ["cutadapt", "--fasta", "--discard-untrimmed"]
    if cpu:
        cmd += ["-j", str(cpu)]
    if min_len:
        cmd += ["-m", str(min_len)]
    if max_len:
        cmd += ["-M", str(max_len)]
    if flip:
        cmd += ["--revcomp"]
    for marker, values in sorted(marker_definitions.items()):
        left_primer = primer_clean(values["left_primer"])
        right_primer = primer_clean(values["right_primer"])
        if not left_primer or not right_primer:
            sys.exit(f"ERROR: Missing primer(s) for marker {marker}")
        cmd += [
            # -a LEFT...RIGHT = left-anchored
            # -g LEFT...RIGHT = non-anchored
            "-g",
            # Here o=... short for min_overlap=...
            f"{marker}={left_primer};o={len(left_primer)}..."
            f"{reverse_complement(right_primer)};o={len(right_primer)}",
        ]
    cmd += [
        "-o",
        out_template,
        long_in,
    ]
    return parse_cutadapt_stdout(run(cmd, debug=debug).stdout)


def parse_flash_stdout(stdout):
    r"""Extract FASTQ pair count before/after running flash.

    >>> parse_flash_stdout("...\n[FLASH] Read combination statistics:[FLASH]     Total pairs:      6105\n[FLASH]     Combined pairs:   5869\n...")
    (6105, 5869)
    """  # noqa: E501
    before = None
    after = None
    for line in stdout.strip().splitlines():
        words = line.strip().split()
        if words[:-1] == ["[FLASH]", "Total", "pairs:"]:
            before = int(words[-1])
        elif words[:-1] == ["[FLASH]", "Combined", "pairs:"]:
            after = int(words[-1])
    if before is None or after is None:
        sys.exit(
            f"ERROR: Could not extract flash before and after pair count:\n\n{stdout}\n"
        )
    return before, after


assert (
    parse_flash_stdout(
        """\
...
[FLASH] Read combination statistics:
[FLASH]     Total pairs:      6105
[FLASH]     Combined pairs:   5869
...
"""
    )
    == (6105, 5869)
)


def run_flash(trimmed_R1, trimmed_R2, output_dir, output_prefix, debug=False, cpu=0):
    """Run FLASH on a pair of trimmed FASTQ files to merge overlapping pairs.

    Returns two integers, FASTQ pair count for input and output files.
    """
    # Note our reads tend to overlap a lot, thus increase max overlap with -M
    # Also, some of our samples are mostly 'outies' rather than 'innies', so -O
    cmd = ["flash", "-O", "-M", "300"]
    if cpu:
        cmd += ["--threads", str(cpu)]
    else:
        cmd += ["-t", "1"]  # Default is all CPUs
    cmd += ["-d", output_dir, "-o", output_prefix, trimmed_R1, trimmed_R2]
    return parse_flash_stdout(run(cmd, debug=debug).stdout)


def save_nr_fasta(
    counts,
    output_fasta,
    min_abundance=0,
    gzipped=False,
    spikes=None,
    header_dict=None,
):
    r"""Save a dictionary of sequences and counts as a FASTA file.

    The output FASTA records are named ``>MD5_abundance\n``, which is the
    default style used in SWARM. This could in future be generalised,
    for example ``>MD5;size=abundance;\n`` for the VSEARCH default.

    If they match a spike-in, they are named ``>MD5_abundance spike-name\n``
    instead.

    Results are sorted by decreasing abundance then alphabetically by
    sequence.

    Returns the total and number of unique sequences accepted (above any
    minimum abundance specified), and a dict of max spike abundances
    (including non-spikes under the empty string).

    Use output_fasta='-' for standard out.
    """
    accepted_total = 0
    accepted_count = 0
    max_spike_abundance = {"": 0}
    if not spikes:
        spikes = []
    # Storing negative count for decreasing sort, then alphabetic sort
    # Could alternatively have used a key function to achieve this
    values = sorted(
        (-count, seq) for seq, count in counts.items() if count >= min_abundance
    )
    if output_fasta == "-":
        if gzipped:
            raise ValueError("Does not support gzipped output to stdout.")
        out_handle = sys.stdout
    elif gzipped:
        out_handle = gzip.open(output_fasta, "wt")
    else:
        out_handle = open(output_fasta, "w")
    if header_dict:
        assert "abundance" not in header_dict
        assert "threshold" not in header_dict
        for key, value in header_dict.items():
            out_handle.write(f"#{key}:{value}\n")
        # Note counts currently negative for sorting requirement
        out_handle.write(f"#abundance:{sum(-count for count, _ in values)}\n")
        out_handle.write(f"#threshold:{min_abundance}\n")
    for spike_name, _, _ in spikes:
        max_spike_abundance[spike_name] = 0
    for count, seq in values:
        count = -count  # was negative for decreasing sorting
        spike_name = is_spike_in(seq, spikes)
        if spike_name:
            out_handle.write(f">{md5seq(seq)}_{count} {spike_name}\n{seq}\n")
        else:
            out_handle.write(f">{md5seq(seq)}_{count}\n{seq}\n")
        accepted_total += count
        accepted_count += 1
        max_spike_abundance[spike_name] = max(
            max_spike_abundance[spike_name],
            count,
        )
    if output_fasta != "-":
        out_handle.close()
    assert accepted_total >= 0 and accepted_count >= 0
    return accepted_total, accepted_count, max_spike_abundance


def kmers(sequence, k=KMER_LENGTH):
    """Make set of all kmers in the given sequence."""
    return {sequence[i : i + k] for i in range(len(sequence) - k + 1)}


def has_enough_kmers(sequence, kmers, threshold=70, k=KMER_LENGTH):
    """Check if given sequence shares at least this many kmers."""
    count = 0
    for i in range(len(sequence) - k + 1):
        if sequence[i : i + k] in kmers:
            count += 1
            if count >= threshold:
                return True
    return False


def is_spike_in(sequence, spikes):
    """Return spike-in name if sequence matches, else empty string."""
    for spike_name, spike_seq, spike_kmers in spikes:
        if sequence == spike_seq:
            return spike_name
        # This will not work when len(spike) <~ kmer length
        # (fail gracefully with an impossible to meet value of 10)
        threshold = min((len(spike_seq) - KMER_LENGTH) / 3, 10)
        if has_enough_kmers(sequence, spike_kmers, threshold):
            return spike_name
    return ""


def make_nr_fasta(
    input_fasta_or_fastq,
    output_fasta,
    min_abundance=0,
    min_len=0,
    max_len=sys.maxsize,
    weighted_input=False,
    fastq=False,
    gzipped=False,
    spikes=None,
    header_dict=None,
    debug=False,
):
    r"""Trim and make non-redundant FASTA/Q file from FASTA input.

    Makes a non-redundant FASTA file with the sequences named
    ``>MD5_abundance\n``.

    For FASTQ files all input reads are treated as abundance one
    (using weighted_input=True gives an error).

    If FASTA input and weighted_input=True, reads must follow
    ``>identifier_abundance\n`` naming and the abundance is used.
    Otherwise all treated as abundance one.

    Makes a non-redundant FASTA file with the sequences named
    ``>MD5_abundance\n``.

    Returns the total number of accepted reads before de-duplication
    (integer), number of those unique (integer), and the total number
    of those which passed the minimum abundance threshold (integer),
    number of those unique (integer), and a dict of maximum abundance
    by spike in (for use with controls for setting the threshold).
    """
    counts = Counter()
    with open(input_fasta_or_fastq) as handle:
        if fastq:
            assert not weighted_input, "Not implemented for FASTQ"
            for _, seq, _ in FastqGeneralIterator(handle):
                if min_len <= len(seq) <= max_len:
                    counts[seq.upper()] += 1
        elif weighted_input:
            for title, seq in SimpleFastaParser(handle):
                assert title.count(" ") == 0 or (
                    title.count(" ") == 1 and title.endswith(" rc")
                ), title
                assert title.count("_") == 1 and title[32] == "_", title
                if min_len <= len(seq) <= max_len:
                    counts[seq.upper()] += abundance_from_read_name(
                        title.split(None, 1)[0]
                    )
        else:
            for _, seq in SimpleFastaParser(handle):
                if min_len <= len(seq) <= max_len:
                    counts[seq.upper()] += 1
    accepted_total, accepted_count, max_spike_abundance = save_nr_fasta(
        counts,
        output_fasta,
        min_abundance,
        gzipped=gzipped,
        spikes=spikes,
        header_dict=header_dict,
    )
    return (
        sum(counts.values()) if counts else 0,
        len(counts),
        accepted_total,
        accepted_count,
        max_spike_abundance,
    )


def merge_paired_reads(
    raw_R1,
    raw_R2,
    merged_fasta_gz,
    tmp,
    debug=False,
    cpu=0,
):
    """Create NR FASTA file by overlap merging the paired FASTQ files."""
    if os.path.isfile(merged_fasta_gz):
        if debug:
            sys.stderr.write(f"DEBUG: Reusing {merged_fasta_gz}\n")
        # Just unzip and read header:
        header = load_fasta_header(merged_fasta_gz, gzipped=True)
        return header["raw_fastq"], header["flash"]
    # flash
    merged_fastq = os.path.join(tmp, "flash.extendedFrags.fastq")
    count_raw, count_flash = run_flash(
        raw_R1, raw_R2, tmp, "flash", debug=debug, cpu=cpu
    )
    if not os.path.isfile(merged_fastq):
        sys.exit(f"ERROR: Expected file {merged_fastq!r} from flash\n")

    if debug:
        sys.stderr.write(f"DEBUG: Caching {merged_fasta_gz}\n")
    tmp_fasta_gz = os.path.join(tmp, "flash.extendedFrags.fasta")
    make_nr_fasta(
        merged_fastq,
        tmp_fasta_gz,
        min_abundance=0,
        fastq=True,
        gzipped=True,
        header_dict={
            # "left_primer": left_primer,
            # "right_primer": right_primer,
            "raw_fastq": count_raw,
            "flash": count_flash,
            # "cutadapt": count_cutadapt,
            # "abundance": accepted_total,
            # "threshold": min_abundance,
        },
    )
    shutil.move(tmp_fasta_gz, merged_fasta_gz)
    del tmp_fasta_gz
    merged_fastq = None
    return count_raw, count_flash


def prepare_sample(
    fasta_name,
    trimmed_fasta,
    marker,
    left_primer,
    right_primer,
    min_len,
    max_len,
    spikes,
    min_abundance,
    min_abundance_fraction,
    control,
    count_raw,
    count_flash,
    tmp,
    debug=False,
    cpu=0,
):
    """Create specified FASTA file for sample from paired FASTQ.

    Does spike-in detection, abundance threshold, and min/max length.

    Returns accepted unique sequence count, accepted total read count, a dict
    of max abundance (keyed by spike-in name), and the absolute abundance
    threshold used.
    """
    if debug:
        sys.stderr.write(f"DEBUG: prepare_sample {trimmed_fasta} --> {fasta_name}\n")
    if os.path.isfile(fasta_name):
        if control:
            uniq_count, total, max_spike_abundance = abundance_values_in_fasta(
                fasta_name
            )
            return uniq_count, total, max_spike_abundance, -1
        else:
            # Don't actually need the max abundance
            return None, None, {}, -1

    if debug:
        sys.stderr.write(
            f"DEBUG: Starting to prepare {'control' if control else 'sample'}"
            f" {fasta_name} (min abundance set to {min_abundance};"
            f" fraction abundance set to {min_abundance_fraction*100}%)\n"
        )

    # count_raw = _header["raw_fastq"]
    # count_flash = _header["flash"]
    # count_tmp, count_cutadapt = run_cutadapt(...)
    if not os.path.isfile(trimmed_fasta):
        sys.exit(f"ERROR: Expected file {trimmed_fasta!r} from cutadapt\n")

    # cutadapt worked on the merged NR reads, so count_cutadapt
    # is currently the unique read count.
    # _unique, _total, _max = abundance_values_in_fasta(merged_fasta_gz, gzipped=True)
    # if _unique != count_tmp:
    #     sys.exit(
    #         f"ERROR: Gave it {_unique} unique sequences, "
    #         f"but cutadapt says saw {count_cutadapt} unique sequences."
    #     )
    _unique, _total, _max = abundance_values_in_fasta(trimmed_fasta)
    count_cutadapt = _total
    del _unique, _total, _max

    if debug:
        sys.stderr.write(
            f"DEBUG: Absolute abundance threshold {min_abundance}, vs"
            f" {min_abundance_fraction*100}% giving"
            f" {ceil(min_abundance_fraction * count_cutadapt)}\n"
        )

    # Using ceiling not floor, as will then take greater-than-or-equal
    min_abundance = max(min_abundance, ceil(min_abundance_fraction * count_cutadapt))

    # deduplicate and apply minimum abundance threshold
    # and tag sequences if they look like spike-ins
    dedup = os.path.join(tmp, "dedup_trimmed.fasta")
    (
        count,
        uniq_count,
        accepted_total,
        accepted_uniq_count,
        max_spike_abundance,
    ) = make_nr_fasta(
        trimmed_fasta,
        dedup,
        min_abundance=min_abundance,
        min_len=min_len,
        max_len=max_len,
        weighted_input=True,
        gzipped=False,
        spikes=spikes,
        header_dict={
            "marker": marker,
            "left_primer": left_primer,
            "right_primer": right_primer,
            "raw_fastq": count_raw,
            "flash": count_flash,
            "cutadapt": count_cutadapt,
            # "abundance": accepted_total,
            # "threshold": min_abundance,
        },
        debug=debug,
    )
    if accepted_total or accepted_uniq_count:
        assert max_spike_abundance, f"Got {max_spike_abundance!r} from {trimmed_fasta}"
    if count_cutadapt < count:
        # Can be less if cutadapt had to use more relaxed global min/max length
        sys.exit(
            f"ERROR: Cutadapt says wrote {count_cutadapt}, but we saw {count} reads."
        )
    assert bool(sum(max_spike_abundance.values())) == bool(accepted_uniq_count)
    if debug:
        sys.stderr.write(
            "DEBUG:"
            f" FASTQ pairs {count_raw}; flash -> {count_flash};"
            f" cutadapt -> {count_cutadapt} [{uniq_count} unique];"
            f" min abundance {min_abundance} -> {accepted_total}"
            f" [{accepted_uniq_count} unique]"
            f", or {accepted_total*100.0/count_raw:0.1f}%\n"
        )
        if accepted_uniq_count:
            sys.stderr.write(
                f"From {count_raw} paired FASTQ reads,"
                f" found {uniq_count} unique sequences,"
                f" {accepted_uniq_count} above min abundance {min_abundance}"
                f" (max abundance {max(max_spike_abundance.values())})\n"
            )

    # File done
    shutil.move(dedup, fasta_name)

    if not accepted_uniq_count:
        if debug:
            sys.stderr.write(
                f"{fasta_name} had {uniq_count} unique marker sequences,"
                f" but none above {'control' if control else 'sample'}"
                f" minimum abundance threshold {min_abundance}\n"
            )
    elif debug:
        sys.stderr.write(
            f"DEBUG: Filtered {fasta_name} down to {accepted_uniq_count} unique"
            f" sequences above {'control' if control else 'sample'}"
            f" min abundance threshold {min_abundance}"
            f" (max abundance {max(max_spike_abundance.values(), default=0)})\n"
        )

    if accepted_uniq_count or accepted_total:
        assert max_spike_abundance, f"Got {max_spike_abundance!r} from {trimmed_fasta}"
    return accepted_uniq_count, accepted_total, max_spike_abundance, min_abundance


def marker_cut(
    marker_definitions,
    file_pairs,
    out_dir,
    merged_cache,
    tmp,
    flip,
    min_abundance,
    min_abundance_fraction,
    debug=False,
    cpu=0,
):
    """Apply primer-trimming for given markers."""
    sys.stderr.write(
        f"Looking for {len(marker_definitions)} markers in {len(file_pairs)} samples\n"
    )

    time_flash = time_cutadapt = time_abundance = 0

    for marker in marker_definitions:
        if not os.path.isdir(os.path.join(out_dir, marker)):
            if debug:
                sys.stderr.write(f"Making {marker} output sub-directory\n")
            os.mkdir(os.path.join(out_dir, marker))

    pool_worst_control = {}  # (marker, folder) as key, max abundance as value
    skipped_samples = set()  # marker specific
    fasta_files_prepared = []  # return value

    for control, stem, raw_R1, raw_R2 in file_pairs:
        if debug:
            sys.stderr.write(f"Preparing {'control' if control else 'sample'} {stem}\n")

        sys.stdout.flush()
        sys.stderr.flush()

        pool_path = os.path.abspath(os.path.split(stem)[0])
        stem = os.path.split(stem)[1]
        if merged_cache:
            merged_fasta_gz = os.path.join(merged_cache, f"{stem}.fasta.gz")
        else:
            # Not told to keep it, just use a temp folder for the merged reads
            if not os.path.isdir(os.path.join(tmp, "merged")):
                os.mkdir(os.path.join(tmp, "merged"))
            merged_fasta_gz = os.path.join(tmp, "merged", f"{stem}.fasta.gz")

        count_raw = count_flash = None  # Won't need if just parsing a control
        if any(
            not os.path.isfile(os.path.join(out_dir, marker, f"{stem}.fasta"))
            for marker in marker_definitions
        ):
            # Run flash to merge reads; or parse pre-existing files
            start = time()
            count_raw, count_flash = merge_paired_reads(
                raw_R1, raw_R2, merged_fasta_gz, tmp, debug=debug, cpu=cpu
            )
            time_flash += time() - start

            # Run cutadapt to cut primers (giving one output per marker)
            start = time()
            unique_merged_, unique_cutadapt_ = run_cutadapt(
                merged_fasta_gz,
                os.path.join(
                    tmp, stem + ".{name}.fasta"
                ),  # template - leave {name} as is!
                marker_definitions,
                flip=flip,
                debug=debug,
                cpu=cpu,
            )
            time_cutadapt += time() - start

        # Apply abundance thresholds
        start = time()
        for marker, marker_values in marker_definitions.items():
            sys.stdout.flush()
            sys.stderr.flush()
            fasta_name = os.path.join(out_dir, marker, f"{stem}.fasta")
            if fasta_name in fasta_files_prepared:
                sys.exit(f"ERROR: Multiple files named {fasta_name}")
            pool_key = (marker, pool_path)
            # Assumes controls before samples in input file list!
            min_a = (
                min_abundance
                if control
                else max(min_abundance, pool_worst_control.get(pool_key, 0))
            )
            if debug and control:
                assert min_a == min_abundance
                sys.stderr.write(
                    f"DEBUG: Control sample so keeping {min_a} as min abundance\n"
                )
            # Will parse pre-existing control file, skips pre-existing samples
            uniq_count, total, max_abundance_by_spike, min_a = prepare_sample(
                fasta_name,
                os.path.join(tmp, f"{stem}.{marker}.fasta"),
                marker,
                marker_values["left_primer"],
                marker_values["right_primer"],
                marker_values["min_length"],
                marker_values["max_length"],
                marker_values["spike_kmers"],
                min_a,
                0  # Not applied to negative controls
                if control
                else min_abundance_fraction,
                control,
                count_raw,
                count_flash,
                tmp,
                debug=debug,
                cpu=cpu,
            )
            fasta_files_prepared.append(fasta_name)
            if uniq_count:
                assert (
                    max_abundance_by_spike
                ), f"Got {max_abundance_by_spike!r} from {fasta_name}"
            # Any spike-in is assumed to be a synthetic control, rest assumed biological
            max_non_spike_abundance = max_abundance_by_spike.get("", 0)
            if control:
                if debug or max_non_spike_abundance > min_abundance:
                    sys.stderr.write(
                        f"Control {stem} max {marker} abundance"
                        f" {max_non_spike_abundance} ({uniq_count} unique"
                        f" sequences, {total} reads, over default"
                        f" threshold {min_abundance})\n"
                    )
                if max_non_spike_abundance > pool_worst_control.get(pool_key, -1):
                    # Record even if zero, nice to have for summary later
                    pool_worst_control[pool_key] = max_non_spike_abundance
                if debug:
                    sys.stderr.write(
                        "Control %s max %s abundance breakdown %s\n"
                        % (
                            stem,
                            marker,
                            ", ".join(
                                f"{k}: {v}"
                                for k, v in sorted(max_abundance_by_spike.items())
                            ),
                        )
                    )
            elif uniq_count is None:
                skipped_samples.add(stem)
                if debug:
                    sys.stderr.write(f"Skipping {fasta_name} as already done\n")
            else:
                sys.stderr.write(
                    f"Sample {stem} has {uniq_count} unique {marker} sequences,"
                    f" {total} reads, over abundance threshold {min_a}"
                    f" (max marker abundance {max_non_spike_abundance})\n"
                )
        time_abundance += time() - start
        if debug:
            sys.stderr.write(
                f"Thus far, {time_flash:0.1f}s running flash and making NR,"
                f" {time_cutadapt:0.1f}s on cutadapt,"
                f" and {time_abundance:0.1f}s applying abundance thresholds\n"
            )

    # Finished all files
    if skipped_samples:
        sys.stderr.write(
            f"Skipped {len(skipped_samples)} previously prepared {marker} samples\n"
        )
    for (marker, pool_path), a in sorted(pool_worst_control.items()):
        if a > min_abundance:
            sys.stderr.write(
                (pool_path if os.path.isabs(pool_path) else os.path.relpath(pool_path))
                + f" {marker} abundance threshold raised to {a}\n"
            )
        else:
            sys.stderr.write(
                (pool_path if os.path.isabs(pool_path) else os.path.relpath(pool_path))
                + f" {marker} negative control abundance {a} (good)\n"
            )

    sys.stderr.write(
        f"Spent {time_flash:0.1f}s running flash and making NR,"
        f" {time_cutadapt:0.1f}s on cutadapt,"
        f" and {time_abundance:0.1f}s applying abundance thresholds\n"
    )

    sys.stdout.flush()
    sys.stderr.flush()

    return fasta_files_prepared


def main(
    fastq,
    negative_controls,
    out_dir,
    session,
    spike_genus,
    flip=False,
    min_abundance=100,
    min_abundance_fraction=0.001,
    ignore_prefixes=None,
    merged_cache=None,
    tmp_dir=None,
    debug=False,
    cpu=0,
):
    """Implement the ``thapbi_pict prepare-reads`` command.

    If there are controls, they will be used to potentially increase
    the minimum abundance threshold used for the non-control files.

    For use in the pipeline command, returns a filename listing of the
    FASTA files created.
    """
    assert isinstance(fastq, list)

    if merged_cache and not os.path.isdir(merged_cache):
        sys.exit(f"ERROR: {merged_cache} for merged cache is not a directory.")

    if out_dir == "-":
        # Can't put files next to FASTQ input when have multiple markers
        # (and mixing raw data with intermediates not a great idea anyway)
        sys.exit("ERROR: Use of output directory '-' no longer supported.")

    check_tools(["flash", "cutadapt"], debug)

    if fastq:
        # Possible in a pipeline setting may need to pass a null value,
        # e.g. -i "" or -i "-" when only have negatives?
        fastq = [_ for _ in fastq if _ and _ != "-"]
    if negative_controls:
        # Possible in a pipeline setting may need to pass a null value,
        # e.g. -n "" or -n "-"
        negative_controls = [_ for _ in negative_controls if _ and _ != "-"]

    # Split on commas, strip white spaces
    spike_genus = [_.strip() for _ in spike_genus.strip().split(",") if _.strip()]
    for x in spike_genus:
        if not session.query(Taxonomy).filter_by(genus=x).count():
            sys.stderr.write(f"WARNING: Spike-in genus {x!r} not in database\n")

    if negative_controls:
        control_file_pairs = find_fastq_pairs(
            negative_controls, ignore_prefixes=ignore_prefixes, debug=debug
        )
    else:
        control_file_pairs = []

    fastq_file_pairs = find_fastq_pairs(
        fastq, ignore_prefixes=ignore_prefixes, debug=debug
    )
    fastq_file_pairs = [_ for _ in fastq_file_pairs if _ not in control_file_pairs]

    if min_abundance < 2:
        # Warning only if on a single file, or on negative controls only
        if len(control_file_pairs) > 1:
            sys.exit(
                "ERROR: Singletons should never be accepted"
                " (except perhaps for debugging controls)"
            )
        else:
            sys.stderr.write(
                "STRONG WARNING: Singletons should never be accepted"
                " (except perhaps for debugging controls).\n"
            )
    elif min_abundance < 10:
        sys.stderr.write(
            "STRONG WARNING: Setting the minimum abundance threshold below 10"
            " is not  advised. You will accept many erroneous reads, and also"
            " slow down the pipeline.\n"
        )
    elif min_abundance < 50 and not min_abundance_fraction:
        sys.stderr.write(
            "WARNING: Only set the minimum abundance threshold below 50 if"
            " you have negative controls to justify this.\n"
        )
    if min_abundance_fraction < 0 or 1 < min_abundance_fraction:
        sys.exit(
            f"ERROR: Minimum abundance fraction {min_abundance_fraction}"
            " should be between zero (no effect) and one (all reads must match)."
        )
    elif 0.1 < min_abundance_fraction:
        sys.stderr.write("WARNING: Minimum abundance fraction should be small\n")

    # Make a unified file list, with control flag
    file_pairs = [
        (True, stem, raw_R1, raw_R2) for stem, raw_R1, raw_R2 in control_file_pairs
    ] + [(False, stem, raw_R1, raw_R2) for stem, raw_R1, raw_R2 in fastq_file_pairs]

    if debug:
        sys.stderr.write(
            f"Preparing {len(fastq_file_pairs)} data FASTQ pairs,"
            f" and {len(control_file_pairs)} control FASTQ pairs\n"
        )
    if control_file_pairs and not fastq_file_pairs:
        sys.stderr.write(
            f"WARNING: {len(control_file_pairs)} control FASTQ pairs,"
            " no non-control reads!\n"
        )

    if out_dir and not os.path.isdir(out_dir):
        sys.stderr.write(f"Making output directory {out_dir!r}\n")
        os.mkdir(out_dir)

    if tmp_dir:
        # Up to the user to remove the files
        tmp_obj = None
        shared_tmp = tmp_dir
    else:
        tmp_obj = tempfile.TemporaryDirectory()
        shared_tmp = tmp_obj.name

    if debug:
        sys.stderr.write(f"DEBUG: Shared temp folder {shared_tmp}\n")

    marker_definitions = {}
    for reference_marker in session.query(MarkerDef).order_by(MarkerDef.name):
        if not reference_marker.left_primer or not reference_marker.right_primer:
            # TODO - ERROR if more than one marker? Always an error?
            sys.exit(f"ERROR: Missing primer(s) for {reference_marker.name}")
        if reference_marker.min_length > reference_marker.max_length:
            sys.exit(
                f"ERROR: Marker {reference_marker.name}"
                f" min length {reference_marker.min_length}"
                f" but max length {reference_marker.max_length}"
            )

        # Spike-in negative controls are marker specific
        spikes = []
        if negative_controls and spike_genus:
            if debug:
                sys.stderr.write(
                    f"DEBUG: Loading any {reference_marker.name} spike-in controls.\n"
                )
            # Doing a join to pull in the marker and taxonomy tables too:
            cur_tax = aliased(Taxonomy)
            marker_seq = aliased(MarkerSeq)
            for seq_source in (
                session.query(SeqSource)
                .join(marker_seq, SeqSource.marker_seq)
                .join(cur_tax, SeqSource.taxonomy)
                .options(contains_eager(SeqSource.marker_seq, alias=marker_seq))
                .options(contains_eager(SeqSource.taxonomy, alias=cur_tax))
                .filter(SeqSource.marker_definition_id == reference_marker.id)
                .order_by(marker_seq.sequence, SeqSource.id)
                .filter(cur_tax.genus.in_(spike_genus))
            ):
                spikes.append(
                    (
                        seq_source.source_accession,
                        seq_source.marker_seq.sequence,
                        kmers(seq_source.marker_seq.sequence),
                    )
                )
            sys.stderr.write(f"Loaded {len(spikes)} spike-in control sequences.\n")

        marker_definitions[reference_marker.name] = {
            "left_primer": reference_marker.left_primer,
            "right_primer": reference_marker.right_primer,
            "min_length": reference_marker.min_length,
            "max_length": reference_marker.max_length,
            "spike_kmers": spikes,
        }

    # Run flash & cutadapt (once doing demultiplexing), apply abundance thresholds
    fasta_files_prepared = marker_cut(
        marker_definitions,
        file_pairs,
        out_dir,
        merged_cache,
        shared_tmp,
        flip,
        min_abundance,
        min_abundance_fraction,
        debug=debug,
        cpu=cpu,
    )

    if tmp_dir:
        sys.stderr.write(
            f"WARNING: Please remove temporary files written to {tmp_dir}\n"
        )
    else:
        tmp_obj.cleanup()

    if debug:
        sys.stderr.write(f"Prepared {len(fasta_files_prepared)} FASTA files\n")
    sys.stdout.flush()
    sys.stderr.flush()
    return fasta_files_prepared
