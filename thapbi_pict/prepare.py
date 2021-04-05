# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Prepare raw ITS1 sequencing reads (trimming, merging, etc).

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

from .hmm import filter_for_hmm
from .utils import abundance_from_read_name
from .utils import abundance_values_in_fasta
from .utils import load_fasta_header
from .utils import md5seq
from .utils import run
from .versions import check_tools


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
        x = os.path.normpath(os.path.relpath(x))
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


def parse_cutadapt_stdout(stdout, bad_primers=False):
    r"""Extract FASTA count before and after cutadapt.

    >>> parse_cutadapt_stdout("...\nTotal reads processed: 5,869\n...\nReads written (passing filters): 5,861 (99.9%)\n...")
    (5869, 5861)
    """  # noqa: E501
    before = None
    after = None
    for line in stdout.strip().split("\n"):
        words = line.strip().split()
        if words == ["No", "reads", "processed!"]:
            before = after = 0
        elif words[:-1] == ["Total", "reads", "processed:"]:
            before = int(words[3].replace(",", ""))
        elif bad_primers and words[:-2] == ["Reads", "with", "adapters:"]:
            after = int(words[3].replace(",", ""))
        elif not bad_primers and words[:4] == [
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
        bad_primers=False,
    )
    == (106089, 1471)
)

# Example recording the bad primers, want lower number
assert (
    parse_cutadapt_stdout(
        """\
=== Summary ===

Total reads processed:                  51,804
Reads with adapters:                    51,601 (99.6%)
Reads that were too short:                   0 (0.0%)
Reads that were too long:                    0 (0.0%)
Reads written (passing filters):        51,804 (100.0%)
""",
        bad_primers=True,
    )
    == (51804, 51601)
)


def run_cutadapt(
    long_in,
    trimmed_out,
    bad_out,
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
    if not left_primer and not right_primer:
        # special case!

        if long_in.endswith(".fasta.gz"):
            # Should the main intermediates be gzipped?
            run(f"cat {long_in} | gunzip > {trimmed_out}", debug=debug)
            _, total, _ = abundance_values_in_fasta(trimmed_out)
            return total, total
        elif long_in.endswith(".fastq"):
            from Bio.SeqIO import convert

            # don't need to parse the names as FASTQ files not using
            # MD5_abundance yet
            total = convert(long_in, "fastq", trimmed_out, "fasta")
            return total, total
        else:
            sys.exit(f"ERROR: called on {long_in} with no primers")
    cmd = ["cutadapt", "--fasta"]
    if bad_out:
        cmd += ["--untrimmed-output", bad_out]
    else:
        cmd += ["--discard-untrimmed"]
        if cpu:
            # Not compatible with --untrimmed-output
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
        f"{left_primer}...{reverse_complement(right_primer)}",
        "-o",
        trimmed_out,
        long_in,
    ]
    return parse_cutadapt_stdout(run(cmd, debug=debug).stdout, bad_out)


def parse_flash_stdout(stdout):
    r"""Extract FASTQ pair count before/after running flash.

    >>> parse_flash_stdout("...\n[FLASH] Read combination statistics:[FLASH]     Total pairs:      6105\n[FLASH]     Combined pairs:   5869\n...")
    (6105, 5869)
    """  # noqa: E501
    before = None
    after = None
    for line in stdout.strip().split("\n"):
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
    counts, output_fasta, min_abundance=0, gzipped=False, header_dict=None
):
    r"""Save a dictionary of sequences and counts as a FASTA file.

    The output FASTA records are named ``>MD5_abundance\n``, which is the
    default style used in SWARM. This could in future be generalised,
    for example ``>MD5;size=abundance;\n`` for the VSEARCH default.

    Results are sorted by decreasing abundance then alphabetically by
    sequence.

    Returns the total and number of unique sequences accepted (above any
    minimum abundance specified).

    Use output_fasta='-' for standard out.
    """
    accepted_total = 0
    accepted_count = 0
    values = sorted((-count, seq) for seq, count in counts.items())
    if output_fasta == "-":
        if gzipped:
            raise ValueError("Does not support gzipped output to stdout.")
        out_handle = sys.stdout
    elif gzipped:
        out_handle = gzip.open(output_fasta, "wt")
    else:
        out_handle = open(output_fasta, "w")
    if header_dict:
        for key, value in header_dict.items():
            out_handle.write(f"#{key}:{value}\n")
    for count, seq in values:
        if -count < min_abundance:
            # Sorted, so everything hereafter is too rare
            break
        out_handle.write(f">{md5seq(seq)}_{-count}\n{seq}\n")
        accepted_total -= count
        accepted_count += 1
    if output_fasta != "-":
        out_handle.close()
    return accepted_total, accepted_count


def make_nr_fasta(
    input_fasta_or_fastq,
    output_fasta,
    min_abundance=0,
    min_len=0,
    max_len=sys.maxsize,
    weighted_input=False,
    fastq=False,
    gzipped=False,
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
    number of those unique (integer), and the maximum abundance of
    any one read (for use with controls for setting the threshold).
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
    accepted_total, accepted_count = save_nr_fasta(
        counts, output_fasta, min_abundance, gzipped=gzipped, header_dict=header_dict
    )
    return (
        sum(counts.values()) if counts else 0,
        len(counts),
        accepted_total,
        accepted_count,
        max(counts.values()) if counts else 0,
    )


def annotate_fasta_with_hmm_and_header(
    input_fasta,
    output_fasta,
    stem,
    shared_tmp_dir,
    hmm_stem=None,
    header_dict=None,
    debug=False,
):
    """Filter for ITS1 regions.

    Assumes you have already applied trimming (and are not using
    the HMM results for this).

    Assumes the SWARM naming convention.

    Returns the number of unique sequences (integer), and
    maximum abundance (dict keyed by HMM name).
    """
    max_hmm_abundance = Counter()
    # This could be generalised if need something else, e.g.
    # >name;size=6; for VSEARCH.
    count = 0
    with open(output_fasta, "w") as out_handle:
        if header_dict:
            for key, value in header_dict.items():
                out_handle.write(f"#{key}:{value}\n")
        for title, full_seq, hmm_name in filter_for_hmm(
            input_fasta, shared_tmp_dir, hmm=hmm_stem, debug=debug
        ):
            if hmm_name is None:
                # Using HMM match(es) to spot control reads
                hmm_name = ""
            if hmm_name:
                out_handle.write(f">{title} {hmm_name}\n{full_seq}\n")
            else:
                out_handle.write(f">{title}\n{full_seq}\n")
            count += 1
            max_hmm_abundance[hmm_name] = max(
                max_hmm_abundance[hmm_name],
                abundance_from_read_name(title.split(None, 1)[0]),
            )

    return count, max_hmm_abundance


def prepare_sample(
    stem,
    raw_R1,
    raw_R2,
    out_dir,
    primer_dir,
    left_primer,
    right_primer,
    flip,
    min_len,
    max_len,
    hmm_stem,
    min_abundance,
    control,
    merged_cache,
    shared_tmp,
    debug=False,
    cpu=0,
):
    """Create FASTA file for sample from paired FASTQ.

    Runs flash, does primer filtering, does HMM filtering.

    Returns fasta filename, unique sequence count, and max abundance
    (dict keyed by HMM name).
    """
    folder, stem = os.path.split(stem)
    if out_dir and out_dir != "-":
        folder = out_dir
    fasta_name = os.path.join(folder, f"{stem}.fasta")
    if primer_dir:
        failed_primer_name = os.path.join(primer_dir, f"{stem}.failed-primers.fasta")
    else:
        failed_primer_name = None

    if os.path.isfile(fasta_name) and (
        failed_primer_name is None or os.path.isfile(failed_primer_name)
    ):
        if control:
            uniq_count, _, max_hmm_abundance = abundance_values_in_fasta(fasta_name)
            return fasta_name, uniq_count, max_hmm_abundance
        else:
            return fasta_name, None, {}

    if debug:
        sys.stderr.write(
            f"DEBUG: Starting to prepare {'control' if control else 'sample'}"
            f" {fasta_name} (min abundance set to {min_abundance})\n"
        )

    tmp = os.path.join(shared_tmp, stem)
    if not os.path.isdir(tmp):
        # If using tempfile.TemporaryDirectory() for shared_tmp
        # this will be deleted automatically, otherwise user must:
        os.mkdir(tmp)

    if debug:
        sys.stderr.write(f"DEBUG: Temp folder of {stem} is {tmp}\n")

    if merged_cache:
        merged_fasta_gz = os.path.join(merged_cache, f"{stem}.fasta.gz")
    else:
        merged_fasta_gz = None

    if merged_fasta_gz and os.path.isfile(merged_fasta_gz):
        if debug:
            sys.stderr.write(f"DEBUG: Reusing {merged_fasta_gz}\n")
        merged_fastq = None
        _header = load_fasta_header(merged_fasta_gz, gzipped=True)
        count_raw = _header["raw_fastq"]
        count_flash = _header["flash"]
        del _header
    else:
        # flash
        merged_fastq = os.path.join(tmp, "flash.extendedFrags.fastq")
        count_raw, count_flash = run_flash(
            raw_R1, raw_R2, tmp, "flash", debug=debug, cpu=cpu
        )
        if not os.path.isfile(merged_fastq):
            sys.exit(f"ERROR: Expected file {merged_fastq!r} from flash\n")

        if merged_fasta_gz:
            if debug:
                sys.stderr.write(f"DEBUG: Caching {merged_fasta_gz}\n")
            tmp_fasta_gz = os.path.join(tmp, "flash.extendedFrags.fasta")
            make_nr_fasta(
                merged_fastq,
                tmp_fasta_gz,
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

    # trim
    trimmed_fasta = os.path.join(tmp, "cutadapt.fasta")
    if failed_primer_name:
        bad_primer_fasta = os.path.join(tmp, "bad_primers.fasta")
    else:
        bad_primer_fasta = None
    # These counts will be of NR FASTA if using merged cache...
    count_tmp, count_cutadapt = run_cutadapt(
        merged_fastq or merged_fasta_gz,
        trimmed_fasta,
        bad_primer_fasta,
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
    if merged_fasta_gz:
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
    else:
        # cutadapt worked on redundant reads
        if count_tmp != count_flash:
            sys.exit(
                f"ERROR: Flash says wrote {count_flash},"
                f" but cutadapt saw {count_tmp} reads."
            )

    # deduplicate and apply minimum abundance threshold
    merged_fasta = os.path.join(tmp, "dedup_trimmed.fasta")
    (
        count,
        uniq_count,
        accepted_total,
        accepted_uniq_count,
        max_hmm_abundance,
    ) = make_nr_fasta(
        trimmed_fasta,
        merged_fasta,
        min_abundance=min_abundance,
        min_len=min_len,
        max_len=max_len,
        weighted_input=bool(merged_fasta_gz),
        gzipped=False,
        debug=debug,
    )
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
            f" (max abundance {max_hmm_abundance})\n"
        )

    if not accepted_uniq_count:
        if debug:
            sys.stderr.write(
                f"{stem} had {uniq_count} unique sequences,"
                f" but none above {'control' if control else 'sample'}"
                f" minimum abundance threshold {min_abundance}\n"
            )
        # Effectively empty FASTA file, just header
        save_nr_fasta(
            {},
            fasta_name,
            header_dict={
                "left_primer": left_primer,
                "right_primer": right_primer,
                "raw_fastq": count_raw,
                "flash": count_flash,
                "cutadapt": count_cutadapt,
                "abundance": accepted_total,
                "threshold": min_abundance,
            },
        )
        if failed_primer_name:
            shutil.move(bad_primer_fasta, failed_primer_name)
        return fasta_name, 0, {}

    # Determine if synthetic controls are present using hmmscan,
    dedup = os.path.join(tmp, "dedup_its1.fasta")
    uniq_count, max_hmm_abundance = annotate_fasta_with_hmm_and_header(
        merged_fasta,
        dedup,
        stem,
        shared_tmp,
        hmm_stem=hmm_stem,
        header_dict={
            "left_primer": left_primer,
            "right_primer": right_primer,
            "raw_fastq": count_raw,
            "flash": count_flash,
            "cutadapt": count_cutadapt,
            "abundance": accepted_total,
            "threshold": min_abundance,
        },
        debug=debug,
    )

    # File done
    shutil.move(dedup, fasta_name)
    if failed_primer_name:
        shutil.move(bad_primer_fasta, failed_primer_name)

    if debug:
        sys.stderr.write(
            f"DEBUG: Filtered {stem} down to {uniq_count} unique sequences"
            f" above {'control' if control else 'sample'}"
            f" min abundance threshold {min_abundance}"
            f" (max abundance {max(max_hmm_abundance.values(), default=0)})\n"
        )

    return fasta_name, uniq_count, max_hmm_abundance


def main(
    fastq,
    negative_controls,
    out_dir,
    hmm_stem,
    primer_dir,
    left_primer,
    right_primer,
    flip=False,
    min_abundance=100,
    min_length=0,
    max_length=sys.maxsize,
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

    if negative_controls:
        # Possible in a pipeline setting may need to pass a null value,
        # e.g. -n "" or -n "-"
        negative_controls = [_ for _ in negative_controls if _ and _ != "-"]

    if negative_controls and not hmm_stem:
        sys.exit("ERROR: If using negative controls, must use --hmm too.")

    check_tools(["flash", "cutadapt"], debug)
    if hmm_stem:
        check_tools(["hmmscan"], debug)

    if negative_controls:
        control_file_pairs = find_fastq_pairs(
            negative_controls, ignore_prefixes=ignore_prefixes, debug=debug
        )
    else:
        control_file_pairs = []

    if not left_primer and not right_primer:
        sys.stderr.write("WARNING: Primer trimming disabled\n")

    fastq_file_pairs = find_fastq_pairs(
        fastq, ignore_prefixes=ignore_prefixes, debug=debug
    )
    fastq_file_pairs = [_ for _ in fastq_file_pairs if _ not in control_file_pairs]

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

    if out_dir and out_dir != "-" and not os.path.isdir(out_dir):
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

    fasta_files_prepared = []  # return value

    skipped_samples = set()
    pool_worst_control = {}  # folder as key, max abundance as value
    for control, stem, raw_R1, raw_R2 in file_pairs:
        sys.stdout.flush()
        sys.stderr.flush()

        pool_key = os.path.abspath(os.path.split(stem)[0])
        min_a = (
            min_abundance
            if control
            else max(min_abundance, pool_worst_control.get(pool_key, 0))
        )
        fasta_file, uniq_count, max_abundance_by_hmm = prepare_sample(
            stem,
            raw_R1,
            raw_R2,
            out_dir,
            primer_dir,
            left_primer,
            right_primer,
            flip,
            min_length,
            max_length,
            hmm_stem,
            min_a,
            control,
            merged_cache,
            shared_tmp,
            debug=debug,
            cpu=cpu,
        )
        if fasta_file in fasta_files_prepared:
            sys.exit(f"ERROR: Multiple files named {fasta_file}")
        fasta_files_prepared.append(fasta_file)
        if uniq_count:
            assert max_abundance_by_hmm, max_abundance_by_hmm
        # Any HMM is assumed to be a synthetic control, no HMM means biological
        max_its1_abundance = max_abundance_by_hmm.get("", 0)
        if control:
            if debug or max_its1_abundance > min_abundance:
                sys.stderr.write(
                    f"Control {stem} has max marker abundance {max_its1_abundance}"
                    f" ({uniq_count} unique sequences over default threshold"
                    f" {min_abundance})\n"
                )
            if max_its1_abundance > pool_worst_control.get(pool_key, -1):
                # Record even if zero, nice to have for summary later
                pool_worst_control[pool_key] = max_its1_abundance
            if debug:
                sys.stderr.write(
                    "Control %s max abundance breakdown %s\n"
                    % (
                        stem,
                        ", ".join(
                            f"{k}: {v}" for k, v in sorted(max_abundance_by_hmm.items())
                        ),
                    )
                )
        elif uniq_count is None:
            skipped_samples.add(stem)
            if debug:
                sys.stderr.write(f"Skipping {fasta_file} as already done\n")
        else:
            sys.stderr.write(
                f"Sample {stem} has {uniq_count} unique sequences over"
                f" abundance threshold {min_a}"
                f" (max marker abundance {max_its1_abundance})\n"
            )

    if skipped_samples:
        sys.stderr.write(
            f"Skipped {len(skipped_samples)} previously prepared samples\n"
        )
    for pool_key in sorted(pool_worst_control):
        if pool_worst_control[pool_key] > min_abundance:
            sys.stderr.write(
                os.path.relpath(pool_key)
                + f" abundance threshold raised to {pool_worst_control[pool_key]}\n"
            )
        else:
            sys.stderr.write(
                os.path.relpath(pool_key)
                + f" negative control abundance {pool_worst_control[pool_key]} (good)\n"
            )

    if tmp_dir:
        sys.stderr.write(
            f"WARNING: Please remove temporary files written to {tmp_dir}\n"
        )
    else:
        tmp_obj.cleanup()

    sys.stdout.flush()
    sys.stderr.flush()
    return fasta_files_prepared
