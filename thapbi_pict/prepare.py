# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

"""Prepare raw ITS1 sequencing reads (trimming, merging, etc).

This implements the ``thapbi_pict prepare-reads ...`` command.
"""

import os
import subprocess
import shutil
import sys
import tempfile

from collections import Counter

from Bio.Seq import reverse_complement
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .hmm import filter_for_ITS1
from .utils import abundance_from_read_name
from .utils import abundance_values_in_fasta
from .utils import md5seq
from .utils import run
from .versions import check_tools


hmm_cropping_warning = 0  # global variable for warning msg


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
                sys.stderr.write("Walking directory %r\n" % x)
            for root, _, files in os.walk(x, followlinks=True):
                for f in files:
                    if f.endswith(ext):
                        if ignore_prefixes and f.startswith(ignore_prefixes):
                            sys.stderr.write(
                                "WARNING: Ignoring %s due to prefix\n"
                                % os.path.join(root, f)
                            )
                            continue
                        answer.append(os.path.join(root, f))
        elif os.path.isfile(x):
            if ignore_prefixes and os.path.split(x)[1].startswith(ignore_prefixes):
                sys.stderr.write("WARNING: Ignoring %s due to prefix\n" % x)
                continue
            answer.append(x)
        else:
            sys.exit("ERROR: %r is not a file or a directory\n" % x)
    # Warn if there were duplicates?
    answer = sorted(set(answer))
    if len(answer) % 2:
        sys.exit("ERROR: Found %i FASTQ files, expected pairs\n" % len(answer))

    pairs = []
    while answer:
        left = answer.pop(0)
        right = answer.pop(0)
        if left.endswith(tuple("_I1_001" + _ for _ in ext)) and right.endswith(
            tuple("_I2_001" + _ for _ in ext)
        ):
            if debug:
                sys.stderr.write(
                    "WARNING: Ignoring %r and %r due to suffix\n" % (left, right)
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
                        "ERROR: Did not recognise %r and %r as a pair\n" % (left, right)
                    )
        if not os.path.isfile(left) and os.path.islink(left):
            sys.stderr.write("WARNING: Ignoring %s as symlink broken\n" % left)
            continue
        if not os.path.isfile(right) and os.path.islink(right):
            sys.stderr.write("WARNING: Ignoring %s as symlink broken\n" % right)
            continue
        for filename in (left, right):
            if not os.path.isfile(filename):
                if os.path.islink(filename):
                    sys.exit("ERROR: Broken symlink: %s" % filename)
                else:
                    # What might cause this - other than deletion after our dir listing?
                    sys.exit("ERROR: Missing file: %s" % filename)
            if not os.stat(filename).st_size:
                if filename.endswith(".gz"):
                    sys.exit("ERROR: Empty gzip file %s" % filename)
                else:
                    sys.stderr.write("WARNING: Empty file %s\n" % filename)
        if not stem:
            sys.exit(
                "ERROR: Did not recognise pair naming for %r and %r\n" % (left, right)
            )
        pairs.append((stem, left, right))

    return pairs


def find_trimmomatic_adapters(fasta_name="TruSeq3-PE.fa"):
    """Locate Illumina adapter FASTA file bundled with trimmomatic."""
    # This works on a bioconda installed trimmomatic,
    # which gives .../conda/bin/trimmomatic and we want
    # .../conda/share/trimmomatic/adapters/TruSeq3-PE.fa
    bin_path = os.path.split(subprocess.getoutput("which trimmomatic"))[0]
    filename = os.path.join(
        bin_path, "..", "share", "trimmomatic", "adapters", fasta_name
    )
    if os.path.isfile(filename):
        return filename
    sys.exit("ERROR: Could not find %s installed with trimmomatic." % fasta_name)


def run_trimmomatic(
    left_in, right_in, left_out, right_out, adapters=None, debug=False, cpu=0
):
    """Run trimmomatic on a pair of FASTQ files.

    The input files may be gzipped.

    If the FASTA adapters file is not specified, will try to use
    ``TruSeq3-PE.fa`` as bundled with a BioConda install of trimmomatic.
    """
    if not adapters:
        adapters = find_trimmomatic_adapters()
    if not os.path.isfile(adapters):
        sys.exit(
            "ERROR: Missing Illumina adapters file for trimmomatic: %s\n" % adapters
        )
    if " " in adapters:
        # Can we do this with slash escaping? Clever quoting?
        sys.exit(
            "ERROR: Spaces in the adapter filename are a bad idea: %s\n" % adapters
        )
    cmd = ["trimmomatic", "PE"]
    if cpu:
        cmd += ["-threads", str(cpu)]
    cmd += ["-phred33"]
    # We don't want the unpaired left and right output files, so /dev/null
    cmd += [left_in, right_in, left_out, os.devnull, right_out, os.devnull]
    cmd += ["ILLUMINACLIP:%s:2:30:10" % adapters]
    return run(cmd, debug=debug)


def run_cutadapt(
    long_in,
    trimmed_out,
    bad_out,
    left_primer,
    right_primer,
    min_len=None,
    max_len=None,
    debug=False,
    cpu=0,
):
    """Run cutadapt on a single file (i.e. after merging paired FASTQ).

    The input and/or output files may be compressed as long as they
    have an appropriate suffix (e.g. gzipped with ``.gz`` suffix).
    """
    cmd = ["cutadapt"]
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
    cmd += [
        # -a LEFT...RIGHT = left-anchored
        # -g LEFT...RIGHT = non-anchored
        "-g",
        "%s...%s" % (left_primer, reverse_complement(right_primer)),
        "-o",
        trimmed_out,
        long_in,
    ]
    return run(cmd, debug=debug)


def run_flash(trimmed_R1, trimmed_R2, output_dir, output_prefix, debug=False, cpu=0):
    """Run pear on a pair of trimmed FASTQ files."""
    # Note our reads tend to overlap a lot, thus increase max overlap with -M
    # Also, some of our samples are mostly 'outies' rather than 'innies', so -O
    cmd = ["flash", "-O", "-M", "300"]
    if cpu:
        cmd += ["--threads", str(cpu)]
    else:
        cmd += ["-t", "1"]  # Default is all CPUs
    cmd += ["-d", output_dir, "-o", output_prefix, trimmed_R1, trimmed_R2]
    return run(cmd, debug=debug)


def save_nr_fasta(counts, output_fasta, min_abundance=0):
    r"""Save a dictionary of sequences and counts as a FASTA file.

    The output FASTA records are named ``>MD5_abundance\n``, which is the
    default style used in SWARM. This could in future be generalised,
    for example ``>MD5;size=abundance;\n`` for the VSEARCH default.

    Results are sorted by decreasing abundance then alphabetically by
    sequence.

    Returns the number of sequences accepted (above any minimum
    abundance specified).

    Use output_fasta='-' for standard out.
    """
    accepted = 0
    values = sorted((-count, seq) for seq, count in counts.items())
    if output_fasta == "-":
        out_handle = sys.stdout
    else:
        out_handle = open(output_fasta, "w")
    for count, seq in values:
        if -count < min_abundance:
            # Sorted, so everything hereafter is too rare
            break
        out_handle.write(">%s_%i\n%s\n" % (md5seq(seq), -count, seq))
        accepted += 1
    if output_fasta != "-":
        out_handle.close()
    return accepted


def make_nr_fasta(
    input_fasta,
    input_rc,
    output_fasta,
    min_abundance=0,
    min_len=0,
    max_len=sys.maxsize,
    debug=False,
):
    r"""Trim and make non-redundant FASTA file from FASTA input.

    The read names are ignored and treated as abundance one!
    Makes a non-redundant FASTA file with the sequences named
    ``>MD5_abundance\n``.

    Returns the number of accepted sequences before de-duplication
    (integer), number of unique accepted sequences (integer), and
    the number of those which passed the minimum abundance threshold,
    and the maximum abundance of any one read (for use with controls
    for setting the threshold).
    """
    counts = Counter()
    with open(input_fasta) as handle:
        for _, seq in SimpleFastaParser(handle):
            assert min_len <= len(seq) <= max_len, "%s len %s" % (_, len(seq))
            counts[seq.upper()] += 1
    if input_rc:
        if debug:
            sys.stderr.write(
                "DEBUG: Combining %r and %r (RC) for unique sequences\n"
                % (input_fasta, input_rc)
            )
        with open(input_rc) as handle:
            for _, seq in SimpleFastaParser(handle):
                assert min_len <= len(seq) <= max_len, "%s len %s (RC)" % (_, len(seq))
                counts[reverse_complement(seq.upper())] += 1
    return (
        sum(counts.values()) if counts else 0,
        len(counts),
        save_nr_fasta(counts, output_fasta, min_abundance),
        max(counts.values()) if counts else 0,
    )


def filter_fasta_for_its1(
    input_fasta, output_fasta, stem, shared_tmp_dir, hmm_stem=None, debug=False
):
    """Filter for ITS1 regions.

    Assumes you have already applied trimming (and are not using
    the HMM results for this).

    Assumes the SWARM naming convention.

    Returns the number of unique sequences (integer),
    maximum abundance (dict keyed by HMM name), and number
    of cropping warnings.
    """
    exp_left = 0
    exp_right = 0
    margin = 10
    cropping_warning = 0
    max_hmm_abundance = Counter()
    # This could be generalised if need something else, e.g.
    # >name;size=6; for VSEARCH.
    count = 0
    with open(output_fasta, "w") as out_handle:
        for title, full_seq, hmm_name, hmm_seqs in filter_for_ITS1(
            input_fasta, shared_tmp_dir, hmm=hmm_stem, debug=debug
        ):
            if not hmm_seqs:
                # Using HMM match(es) as a presense/absense filter
                continue
            assert hmm_name or not hmm_stem, hmm_name
            out_handle.write(">%s %s\n%s\n" % (title, hmm_name, full_seq))
            count += 1
            max_hmm_abundance[hmm_name] = max(
                max_hmm_abundance[hmm_name],
                abundance_from_read_name(title.split(None, 1)[0]),
            )

            # Now some debugging...
            if len(hmm_seqs) == 1:
                hmm_seq = hmm_seqs[0]
                left = full_seq.index(hmm_seq)
                right = len(full_seq) - left - len(hmm_seq)
                if not (
                    exp_left - margin < left < exp_left + margin
                    and exp_right - margin < right < exp_right + margin
                ):
                    cropping_warning += 1
                    sys.stderr.write(
                        "WARNING: %r has HMM cropping %i left, %i right, "
                        "giving %i, vs %i bp from fixed trimming\n"
                        % (
                            title.split(None, 1)[0],
                            left,
                            right,
                            len(hmm_seq),
                            len(full_seq) - exp_left - exp_right,
                        )
                    )
            else:
                sys.stderr.write(
                    "WARNING: %i %s HMM matches from %s\n"
                    % (len(hmm_seqs), hmm_name, title)
                )

    if cropping_warning:
        sys.stderr.write(
            "WARNING: HMM cropping very different from fixed trimming "
            "in %i sequences in %s\n" % (cropping_warning, stem)
        )

    return count, max_hmm_abundance, cropping_warning


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
    shared_tmp,
    debug=False,
    cpu=0,
):
    """Create FASTA file for sample from paired FASTQ.

    Runs trimmomatic, pear, does primer filtering, does HMM filtering.

    Returns fasta filename,  unique sequence count, max abundance
    (dict keyed by HMM name), and hmm_cropping_warning count.
    """
    folder, stem = os.path.split(stem)
    if out_dir and out_dir != "-":
        folder = out_dir
    fasta_name = os.path.join(folder, "%s.fasta" % stem)
    if primer_dir:
        failed_primer_name = os.path.join(primer_dir, "%s.failed-primers.fasta" % stem)
    else:
        failed_primer_name = None

    if os.path.isfile(fasta_name) and (
        failed_primer_name is None or os.path.isfile(failed_primer_name)
    ):
        if control:
            uniq_count, max_hmm_abundance = abundance_values_in_fasta(fasta_name)
            return fasta_name, uniq_count, max_hmm_abundance, 0
        else:
            return fasta_name, None, {}, 0

    if debug:
        sys.stderr.write(
            "DEBUG: Starting to prepare %s %s (min abundance set to %i)\n"
            % ("control" if control else "sample", fasta_name, min_abundance)
        )

    tmp = os.path.join(shared_tmp, stem)
    if not os.path.isdir(tmp):
        # If using tempfile.TemporaryDirectory() for shared_tmp
        # this will be deleted automatically, otherwise user must:
        os.mkdir(tmp)

    if debug:
        sys.stderr.write("DEBUG: Temp folder of %s is %s\n" % (stem, tmp))

    # trimmomatic
    trim_R1 = os.path.join(tmp, "trimmomatic_R1.fastq")
    trim_R2 = os.path.join(tmp, "trimmomatic_R2.fastq")
    run_trimmomatic(raw_R1, raw_R2, trim_R1, trim_R2, debug=debug, cpu=cpu)
    for _ in (trim_R1, trim_R2):
        if not os.path.isfile(_):
            sys.exit("ERROR: Expected file %r from trimmomatic\n" % _)

    # flash
    merged_fastq = os.path.join(tmp, "flash.extendedFrags.fastq")
    run_flash(trim_R1, trim_R2, tmp, "flash", debug=debug, cpu=cpu)
    if not os.path.isfile(merged_fastq):
        sys.exit("ERROR: Expected file %r from flash\n" % merged_fastq)

    # trim
    trimmed_fasta = os.path.join(tmp, "cutadapt.fasta")
    if flip or failed_primer_name:
        bad_primer_fasta = os.path.join(tmp, "bad_primers.fasta")
    else:
        bad_primer_fasta = None
    run_cutadapt(
        merged_fastq,
        trimmed_fasta,
        bad_primer_fasta,
        left_primer,
        right_primer,
        min_len=min_len,
        max_len=max_len,
        debug=debug,
        cpu=cpu,
    )
    if not os.path.isfile(trimmed_fasta):
        sys.exit("ERROR: Expected file %r from cutadapt\n" % trimmed_fasta)

    if flip:
        # Call cutadapt again - with primers reversed, and flip output,
        flipped_fasta = os.path.join(tmp, "cutadapt_flipped.fasta")
        if failed_primer_name:
            bad_primer_fasta2 = os.path.join(tmp, "bad_primers2.fasta")
        else:
            bad_primer_fasta2 = None
        run_cutadapt(
            bad_primer_fasta,
            flipped_fasta,
            bad_primer_fasta2,
            right_primer,
            left_primer,
            min_len=min_len,
            max_len=max_len,
            debug=debug,
            cpu=cpu,
        )
        if not os.path.isfile(trimmed_fasta):
            sys.exit("ERROR: Expected file %r from cutadapt\n" % flipped_fasta)
    else:
        flipped_fasta = None

    # deduplicate and apply minimum abundance threshold
    merged_fasta = os.path.join(tmp, "dedup_trimmed.fasta")
    (count, uniq_count, acc_uniq_count, max_hmm_abundance) = make_nr_fasta(
        trimmed_fasta,
        flipped_fasta,
        merged_fasta,
        min_abundance=min_abundance,
        min_len=min_len,
        max_len=max_len,
        debug=debug,
    )
    if debug:
        sys.stderr.write(
            "Merged %i paired FASTQ reads into %i unique sequences, "
            "%i above min abundance %i (max abundance %i)\n"
            % (count, uniq_count, acc_uniq_count, min_abundance, max_hmm_abundance)
        )

    if not acc_uniq_count:
        if debug:
            sys.stderr.write(
                "%s had %i unique sequences, "
                "but none above %s minimum abundance threshold %i\n"
                % (stem, uniq_count, "control" if control else "sample", min_abundance)
            )
        with open(fasta_name, "w"):
            # Write empty file
            pass
        if failed_primer_name:
            shutil.move(bad_primer_fasta, failed_primer_name)
        return fasta_name, 0, {}, 0

    if debug:
        sys.stderr.write(
            "Merged %s %i paired FASTQ reads into %i unique sequences, "
            "%i above %s min abundance %i (max abundance %i)\n"
            % (
                stem,
                count,
                uniq_count,
                acc_uniq_count,
                "control" if control else "sample",
                min_abundance,
                max_hmm_abundance,
            )
        )

    # Determine if ITS1 region is present using hmmscan,
    dedup = os.path.join(tmp, "dedup_its1.fasta")
    uniq_count, max_hmm_abundance, cropping = filter_fasta_for_its1(
        merged_fasta, dedup, stem, shared_tmp, hmm_stem=hmm_stem, debug=debug
    )

    # File done
    shutil.move(dedup, fasta_name)
    if failed_primer_name:
        shutil.move(bad_primer_fasta, failed_primer_name)

    if debug:
        sys.stderr.write(
            "DEBUG: Filtered %s down to %i unique sequences "
            "above %s min abundance threshold %i (max abundance %i)\n"
            % (
                stem,
                uniq_count,
                "control" if control else "sample",
                min_abundance,
                max(max_hmm_abundance.values(), default=0),
            )
        )

    return fasta_name, uniq_count, max_hmm_abundance, cropping


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
    hmm_cropping_warning = 0

    assert isinstance(fastq, list)

    check_tools(["trimmomatic", "flash", "cutadapt"], debug)
    if hmm_stem:
        check_tools(["hmmscan"], debug)

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

    # Make a unified file list, with control flag
    file_pairs = [
        (True, stem, raw_R1, raw_R2) for stem, raw_R1, raw_R2 in control_file_pairs
    ] + [(False, stem, raw_R1, raw_R2) for stem, raw_R1, raw_R2 in fastq_file_pairs]

    if debug:
        sys.stderr.write(
            "Preparing %i data FASTQ pairs, and %i control FASTQ pairs\n"
            % (len(fastq_file_pairs), len(control_file_pairs))
        )
    if control_file_pairs and not fastq_file_pairs:
        sys.stderr.write(
            "WARNING: %i control FASTQ pairs, no non-control reads!\n"
            % len(control_file_pairs)
        )

    if out_dir and out_dir != "-" and not os.path.isdir(out_dir):
        sys.stderr.write("Making output directory %r\n" % out_dir)
        os.mkdir(out_dir)

    if tmp_dir:
        # Up to the user to remove the files
        tmp_obj = None
        shared_tmp = tmp_dir
    else:
        tmp_obj = tempfile.TemporaryDirectory()
        shared_tmp = tmp_obj.name

    if debug:
        sys.stderr.write("DEBUG: Shared temp folder %s\n" % shared_tmp)

    fasta_files_prepared = []  # return value

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
        fasta_file, uniq_count, max_abundance_by_hmm, hmm_cropping = prepare_sample(
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
            shared_tmp,
            debug=debug,
            cpu=cpu,
        )
        if fasta_file in fasta_files_prepared:
            sys.exit("ERROR: Multiple files named %s" % fasta_file)
        fasta_files_prepared.append(fasta_file)
        if uniq_count:
            assert max_abundance_by_hmm, max_abundance_by_hmm
        hmm_cropping_warning += hmm_cropping
        max_its1_abundance = max(
            (
                max_abundance_by_hmm[_]
                for _ in max_abundance_by_hmm
                if not _.startswith("SynCtrl")
            ),
            default=0,
        )
        if control:
            sys.stderr.write(
                "Control %s has %i unique sequences "
                "over control abundance threshold %i (max marker abundance %i)\n"
                % (stem, uniq_count, min_abundance, max_its1_abundance)
            )
            if max_its1_abundance > pool_worst_control.get(pool_key, 0):
                pool_worst_control[pool_key] = max_its1_abundance
            if debug:
                sys.stderr.write(
                    "Control %s max abundance breakdown %s\n"
                    % (
                        stem,
                        ", ".join(
                            "%s: %i" % (k, v)
                            for k, v in sorted(max_abundance_by_hmm.items())
                        ),
                    )
                )
        elif uniq_count is None:
            sys.stderr.write("Sample %s already done\n" % stem)
        else:
            sys.stderr.write(
                "Sample %s has %i unique sequences "
                "over abundance threshold %i (max marker abundance %i)\n"
                % (stem, uniq_count, min_a, max_its1_abundance)
            )

    if tmp_dir:
        sys.stderr.write(
            "WARNING: Please remove temporary files written to %s\n" % tmp_dir
        )
    else:
        tmp_obj.cleanup()

    if hmm_cropping_warning:
        sys.stderr.write(
            "WARNING: HMM cropping very different from fixed trimming "
            "in %i sequences over %i paired FASTQ files\n"
            % (hmm_cropping_warning, len(file_pairs))
        )

    sys.stdout.flush()
    sys.stderr.flush()
    return fasta_files_prepared
