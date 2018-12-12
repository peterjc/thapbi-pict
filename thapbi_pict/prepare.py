"""Prepare raw ITS1 sequencing reads (trimming, merging, etc).

This implementes the ``thapbi_pict prepare-reads ...`` command.
"""

import os
import sys
import tempfile

from .utils import run


def find_fastq_pairs(filenames_or_folders, ext=(".fastq", ".fastq.gz"),
                     debug=False):
    """Interpret a list of filenames and/or foldernames.

    Returns a list of tuples (stem, left filename, right filename)
    where stem is intended for use in logging and output naming,
    and may include a directory name.
    """
    answer = []
    for x in filenames_or_folders:
        if os.path.isdir(x):
            if debug:
                sys.stderr.write("Walking directory %r\n" % x)
            for f in os.listdir(x):
                if f.endswith(ext):
                    # Check not a directory?
                    answer.append(os.path.join(x, f))
        elif os.path.isfile(x):
            answer.append(x)
        else:
            sys.exit("ERROR: %r is not a file or a directory\n" % x)
    # Warn if there were duplicates?
    answer = sorted(set(answer))
    if len(answer) % 2:
        sys.exit(
            "ERROR: Found %i FASTQ files, expected pairs\n"
            % len(answer))

    pairs = []
    while answer:
        left = answer.pop(0)
        right = answer.pop(0)
        if (left.endswith(tuple("_I1_001" + _ for _ in ext))
                and right.endswith(tuple("_I2_001" + _ for _ in ext))):
            if debug:
                sys.stderr.write(
                    "WARNING: Ignoring %r and %r\n" % (left, right))
            continue
        if not (left.endswith(tuple("_R1_001" + _ for _ in ext)) and
                right.endswith(tuple("_R2_001" + _ for _ in ext))):
            sys.exit(
                "ERROR: Did not recognise pair naming for %r and %r\n"
                % (left, right))
        stem = left.split("_R1_001")[0]
        if stem != right.split("_R2_001")[0]:
            sys.exit(
                "ERROR: Did not recognise %r and %r as a pair\n"
                % (left, right))
        pairs.append((stem, left, right))

    return pairs


def run_trimmomatic(left_in, right_in, left_out, right_out,
                    adapters="TruSeq3-PE.fa",
                    debug=False, cpu=0):
    """Run trimmomatic on a pair of FASTQ files.

    The input files may be gzipped.
    """
    if not os.path.isfile(adapters):
        sys.exit("ERROR: Missing Illumina adapters file for trimmomatic: %s\n"
                 % adapters)
    if " " in adapters:
        # Can we do this with slash escaping? Clever quoting?
        sys.exit("ERROR: Spaces in the adapter filename are a bad idea: %s\n"
                 % adapters)
    cmd = ["trimmomatic", "PE"]
    if cpu:
        cmd += ["-threads", str(cpu)]
    cmd += ["-phred33"]
    # We don't want the unpaired left and right output files, so /dev/null
    cmd += [left_in, right_in, left_out, os.devnull, right_out, os.devnull]
    cmd += ['ILLUMINACLIP:%s:2:30:10' % adapters]
    return run(cmd, debug=debug)


def run_pear(trimmed_R1, trimmed_R2, output_prefix,
             debug=False, cpu=0):
    """Run pear on a pair of trimmed FASTQ files."""
    cmd = ["pear", "-f", trimmed_R1, "-r", trimmed_R2, "-o", output_prefix]
    if cpu:
        cmd += ["--threads", str(cpu)]
    return run(cmd, debug=debug)


def main(fastq, out_dir, debug=False, cpu=0):
    """Implement the thapbi_pict prepare-reads command."""
    assert isinstance(fastq, list)

    fastq_file_pairs = find_fastq_pairs(fastq, debug=debug)
    if debug:
        sys.stderr.write(
            "Preparing %i paired FASTQ files\n"
            % len(fastq_file_pairs))

    for stem, raw_R1, raw_R2 in fastq_file_pairs:
        folder, stem = os.path.split(stem)
        if out_dir and out_dir != "-":
            folder = out_dir
        fasta_name = os.path.join(
            folder, "%s.prepared.fasta" % stem)
        if os.path.isfile(fasta_name):
            sys.stderr.write(
                "WARNING: Skipping %s as already exist\n" % fasta_name)
            continue

        # Context manager should remove the temp dir:
        with tempfile.TemporaryDirectory() as tmp:
            if debug:
                sys.stderr.write(
                    "DEBUG: Temp folder of %s is %s\n" % (stem, tmp))

            # trimmomatic
            trim_R1 = os.path.join(tmp, "trimmomatic_R1.fastq")
            trim_R2 = os.path.join(tmp, "trimmomatic_R2.fastq")
            run_trimmomatic(
                raw_R1, raw_R2, trim_R1, trim_R2,
                debug=debug, cpu=cpu)
            for _ in (trim_R1, trim_R2):
                if not os.path.isfile(_):
                    sys.exit("ERROR: Expected file %r from trimmomatic\n" % _)

            # pear
            pear_prefix = os.path.join(tmp, "pear")
            merged = os.path.join(tmp, "pear.assembled.fastq")
            run_pear(
                trim_R1, trim_R2, pear_prefix,
                debug=debug, cpu=cpu)
            if not os.path.isfile(merged):
                sys.exit("ERROR: Expected file %r from pear\n" % merged)

            # TODO - Turn the merged FASTQ into a FASTA files...

        sys.stderr.write("Prepared %s\n" % stem)
    return 0
