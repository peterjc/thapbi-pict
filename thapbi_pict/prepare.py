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

from Bio.Seq import reverse_complement
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from sqlalchemy.orm import aliased
from sqlalchemy.orm import contains_eager

from .db_orm import connect_to_db
from .db_orm import MarkerDef
from .db_orm import MarkerSeq
from .db_orm import SeqSource
from .db_orm import Taxonomy
from .utils import abundance_from_read_name
from .utils import abundance_values_in_fasta
from .utils import load_fasta_header
from .utils import md5seq
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
    trimmed_out,
    left_primer,
    right_primer,
    min_len=None,
    max_len=None,
    flip=False,
    debug=False,
    cpu=0,
):
    """Run cutadapt on a single file (i.e. after merging paired FASTQ).

    The input and/or output files may be compressed as long as they
    have an appropriate suffix (e.g. gzipped with ``.gz`` suffix).

    Returns FASTA count before and after cutadapt.
    """
    if not left_primer or not right_primer:
        sys.exit("ERROR: Can't run cutadapt without two primers")
    cmd = ["cutadapt", "--fasta", "--discard-untrimmed"]
    if cpu:
        cmd += ["-j", str(cpu)]
    if min_len:
        cmd += ["-m", str(min_len)]
    if max_len:
        cmd += ["-M", str(max_len)]
    if flip:
        cmd += ["--revcomp"]
    cmd += [
        # -a LEFT...RIGHT = left-anchored
        # -g LEFT...RIGHT = non-anchored
        "-g",
        # Here o=... short for min_overlap=...
        f"{left_primer};o={len(left_primer)}..."
        f"{reverse_complement(right_primer)};o={len(right_primer)}",
        "-o",
        trimmed_out,
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
    """Run FLASH on a pair of trimmed FASTQ files to merge overlaping pairs.

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
    (incuding non-spikes under the empty string).

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
                assert min_len <= len(seq) <= max_len, f"{_} len {len(seq)}"
                counts[seq.upper()] += 1
        elif weighted_input:
            for title, seq in SimpleFastaParser(handle):
                assert min_len <= len(seq) <= max_len, f"{_} len {len(seq)}"
                assert title.count(" ") == 0 or (
                    title.count(" ") == 1 and title.endswith(" rc")
                ), title
                assert title.count("_") == 1 and title[32] == "_", title
                counts[seq.upper()] += abundance_from_read_name(title.split(None, 1)[0])
        else:
            for _, seq in SimpleFastaParser(handle):
                assert min_len <= len(seq) <= max_len, f"{_} len {len(seq)}"
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


def annotate_fasta_with_spike_and_header(
    input_fasta,
    output_fasta,
    stem,
    shared_tmp_dir,
    spikes=None,
    header_dict=None,
    debug=False,
):
    """Annotate FASTA file with header and any spike-in matches.

    Assumes you have already applied trimming.

    Assumes the SWARM naming convention.

    Returns the number of unique sequences (integer), total read count
    (integer) and maximum abundance (dict keyed by spike name).
    """
    if not spikes:
        spikes = []

    max_spike_abundance = Counter()
    # This could be generalised if need something else, e.g.
    # >name;size=6; for VSEARCH.
    count = 0
    total = 0
    with open(output_fasta, "w") as out_handle:
        if header_dict:
            for key, value in header_dict.items():
                out_handle.write(f"#{key}:{value}\n")
        with open(input_fasta) as handle:
            for title, full_seq in SimpleFastaParser(handle):
                spike_name = is_spike_in(full_seq, spikes)
                if spike_name:
                    out_handle.write(f">{title} {spike_name}\n{full_seq}\n")
                else:
                    out_handle.write(f">{title}\n{full_seq}\n")
                count += 1
                abundance = abundance_from_read_name(title.split(None, 1)[0])
                total += abundance
                max_spike_abundance[spike_name] = max(
                    max_spike_abundance[spike_name],
                    abundance,
                )

    return count, total, max_spike_abundance


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
        return
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


def prepare_sample(
    fasta_name,
    merged_fasta_gz,
    left_primer,
    right_primer,
    flip,
    min_len,
    max_len,
    spikes,
    min_abundance,
    control,
    tmp,
    debug=False,
    cpu=0,
):
    """Create FASTA file for sample from paired FASTQ.

    Runs flash, does primer filtering, does spike-in detection.

    Returns fasta filename, accepted unique sequence count, accepted total
    read count, and a dict of max abundance (keyed by spike-in name).
    """
    if os.path.isfile(fasta_name):
        if control:
            uniq_count, total, max_spike_abundance = abundance_values_in_fasta(
                fasta_name
            )
            return fasta_name, uniq_count, total, max_spike_abundance
        else:
            # Don't actually need the max abundance
            return fasta_name, None, None, {}

    if debug:
        sys.stderr.write(
            f"DEBUG: Starting to prepare {'control' if control else 'sample'}"
            f" {fasta_name} (min abundance set to {min_abundance})\n"
        )

    # Just unzip and read header:
    _header = load_fasta_header(merged_fasta_gz, gzipped=True)
    count_raw = _header["raw_fastq"]
    count_flash = _header["flash"]
    del _header

    # trim
    trimmed_fasta = os.path.join(tmp, "cutadapt.fasta")
    # These counts will be of NR FASTA if using merged cache...
    count_tmp, count_cutadapt = run_cutadapt(
        merged_fasta_gz,
        trimmed_fasta,
        left_primer,
        right_primer,
        min_len=min_len,
        max_len=max_len,
        flip=flip,
        debug=debug,
        cpu=cpu,
    )
    if not os.path.isfile(trimmed_fasta):
        sys.exit(f"ERROR: Expected file {trimmed_fasta!r} from cutadapt\n")

    # cutadapt worked on the merged NR reads, so count_cutadapt
    # is currently the unique read count.
    _unique, _total, _max = abundance_values_in_fasta(merged_fasta_gz, gzipped=True)
    if _unique != count_tmp:
        sys.exit(
            f"ERROR: Gave it {_unique} unique sequences, "
            f"but cutadapt says saw {count_cutadapt} unique sequences."
        )
    _unique, _total, _max = abundance_values_in_fasta(trimmed_fasta)
    if _unique != count_cutadapt:
        sys.exit(
            f"ERROR: Cutadapt says wrote {count_cutadapt} unique sequences, "
            f"but we see {_unique} unique sequences."
        )
    count_cutadapt = _total
    del _unique, _total, _max

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
        weighted_input=bool(merged_fasta_gz),
        gzipped=False,
        spikes=spikes,
        header_dict={
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
    if count_cutadapt != count:
        sys.exit(
            f"ERROR: Cutadapt says wrote {count_cutadapt}, but we saw {count} reads."
        )
    if debug:
        sys.stderr.write(
            f"DEBUG: FASTQ pairs {count_raw};"
            f" flash -> {count_flash};"
            f" cutadapt -> {count_cutadapt} [{uniq_count} unique];"
            f" abundance -> {accepted_total} [{accepted_uniq_count} unique]"
            f", or {accepted_total*100.0/count_raw:0.1f}%\n"
        )
        sys.stderr.write(
            f"From {count_raw} paired FASTQ reads,"
            f" found {uniq_count} unique sequences,"
            f" {accepted_uniq_count} above min abundance {min_abundance}"
            f" (max abundance {max(max_spike_abundance.values())})\n"
        )

    if not accepted_uniq_count:
        if debug:
            sys.stderr.write(
                f"{fasta_name} had {uniq_count} unique sequences,"
                f" but none above {'control' if control else 'sample'}"
                f" minimum abundance threshold {min_abundance}\n"
            )

    # File done
    shutil.move(dedup, fasta_name)

    if debug:
        sys.stderr.write(
            f"DEBUG: Filtered {fasta_name} down to {uniq_count} unique sequences"
            f" above {'control' if control else 'sample'}"
            f" min abundance threshold {min_abundance}"
            f" (max abundance {max(max_spike_abundance.values(), default=0)})\n"
        )

    if accepted_uniq_count or accepted_total:
        assert max_spike_abundance, f"Got {max_spike_abundance!r} from {trimmed_fasta}"
    return fasta_name, accepted_uniq_count, accepted_total, max_spike_abundance


def marker_cut(
    file_pairs,
    out_dir,
    merged_cache,
    tmp,
    marker,
    left_primer,
    right_primer,
    min_length,
    max_length,
    flip,
    min_abundance,
    spikes,
    debug=False,
    cpu=0,
):
    """Apply primer-trimming for a given marker."""
    sys.stderr.write(f"Looking for {marker} sequences...\n")
    if not left_primer and not right_primer:
        # TODO - ERROR if more than one marker? Always an error?
        sys.stderr.write(f"WARNING: {marker} primer trimming disabled\n")
    if min_length > max_length:
        sys.exit(
            f"ERROR: Marker {marker} has "
            f"min length {min_length} but max length {max_length}"
        )

    fasta_files_prepared = []  # return value
    skipped_samples = set()
    pool_worst_control = {}  # folder as key, max abundance as value
    for control, stem, raw_R1, raw_R2 in file_pairs:
        sys.stdout.flush()
        sys.stderr.flush()

        pool_key = os.path.abspath(os.path.split(stem)[0])
        stem = os.path.split(stem)[1]
        fasta_name = os.path.join(out_dir, f"{stem}.fasta")  # insert marker here
        if merged_cache:
            merged_fasta_gz = os.path.join(merged_cache, f"{stem}.fasta.gz")
        else:
            # Not told to keep it, just use a temp folder for the merged reads
            if not os.path.isdir(os.path.join(tmp, "merged")):
                os.mkdir(os.path.join(tmp, "merged"))
            merged_fasta_gz = os.path.join(tmp, "merged", f"{stem}.fasta.gz")
        merge_paired_reads(raw_R1, raw_R2, merged_fasta_gz, tmp, debug=debug, cpu=cpu)

        min_a = (
            min_abundance
            if control
            else max(min_abundance, pool_worst_control.get(pool_key, 0))
        )
        fasta_file, uniq_count, total, max_abundance_by_spike = prepare_sample(
            fasta_name,
            merged_fasta_gz,
            left_primer,
            right_primer,
            flip,
            min_length,
            max_length,
            spikes,
            min_a,
            control,
            tmp,
            debug=debug,
            cpu=cpu,
        )
        if fasta_file in fasta_files_prepared:
            sys.exit(f"ERROR: Multiple files named {fasta_file}")
        fasta_files_prepared.append(fasta_file)
        if uniq_count:
            assert (
                max_abundance_by_spike
            ), f"Got {max_abundance_by_spike!r} from {fasta_file}"
        # Any spike-in is assumed to be a synthetic control, rest assumed biological
        max_non_spike_abundance = max_abundance_by_spike.get("", 0)
        if control:
            if debug or max_non_spike_abundance > min_abundance:
                sys.stderr.write(
                    f"Control {stem} has max marker abundance {max_non_spike_abundance}"
                    f" ({uniq_count} unique sequences, {total} reads, over default"
                    f" threshold {min_abundance})\n"
                )
            if max_non_spike_abundance > pool_worst_control.get(pool_key, -1):
                # Record even if zero, nice to have for summary later
                pool_worst_control[pool_key] = max_non_spike_abundance
            if debug:
                sys.stderr.write(
                    "Control %s max abundance breakdown %s\n"
                    % (
                        stem,
                        ", ".join(
                            f"{k}: {v}"
                            for k, v in sorted(max_abundance_by_spike.items())
                        ),
                    )
                )
        elif uniq_count is None:
            skipped_samples.add(stem)
            if debug:
                sys.stderr.write(f"Skipping {fasta_file} as already done\n")
        else:
            sys.stderr.write(
                f"Sample {stem} has {uniq_count} unique sequences, {total}"
                f" reads, over abundance threshold {min_a}"
                f" (max marker abundance {max_non_spike_abundance})\n"
            )

    if skipped_samples:
        sys.stderr.write(
            f"Skipped {len(skipped_samples)} previously prepared samples\n"
        )
    for pool_key in sorted(pool_worst_control):
        if pool_worst_control[pool_key] > min_abundance:
            sys.stderr.write(
                (pool_key if os.path.isabs(pool_key) else os.path.relpath(pool_key))
                + f" abundance threshold raised to {pool_worst_control[pool_key]}\n"
            )
        else:
            sys.stderr.write(
                (pool_key if os.path.isabs(pool_key) else os.path.relpath(pool_key))
                + f" negative control abundance {pool_worst_control[pool_key]} (good)\n"
            )

    return fasta_files_prepared


def main(
    fastq,
    negative_controls,
    out_dir,
    db_url,
    spike_genus,
    flip=False,
    min_abundance=100,
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

    # Connect to the DB,
    Session = connect_to_db(db_url)  # echo=debug
    session = Session()

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
    elif min_abundance < 50:
        sys.stderr.write(
            "WARNING: Only set the minimum abundance threshold below 50 if"
            " you have negative controls to justify this.\n"
        )

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

    if not os.path.isdir(out_dir):
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

    fasta_files_prepared = []
    for reference_marker in session.query(MarkerDef).order_by(MarkerDef.name):
        if not os.path.isdir(os.path.join(out_dir, reference_marker.name)):
            sys.stderr.write(f"Making {reference_marker.name} output sub-directory\n")
            os.mkdir(os.path.join(out_dir, reference_marker.name))
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

        fasta_files_prepared.extend(
            marker_cut(
                file_pairs,
                os.path.join(out_dir, reference_marker.name),
                merged_cache,
                shared_tmp,
                reference_marker.name,
                reference_marker.left_primer,
                reference_marker.right_primer,
                reference_marker.min_length,
                reference_marker.max_length,
                flip,
                min_abundance,
                spikes,
                debug=debug,
                cpu=cpu,
            )
        )

    if tmp_dir:
        sys.stderr.write(
            f"WARNING: Please remove temporary files written to {tmp_dir}\n"
        )
    else:
        tmp_obj.cleanup()

    session.close()

    sys.stdout.flush()
    sys.stderr.flush()
    return fasta_files_prepared
