"""Prepare raw ITS1 sequencing reads (trimming, merging, etc).

This implementes the ``thapbi_pict prepare-reads ...`` command.
"""

import os
import sys


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


def main(fastq, out_dir, debug=False):
    """Implement the thapbi_pict prepare-reads command."""
    assert isinstance(fastq, list)

    fastq_file_pairs = find_fastq_pairs(fastq, debug=debug)
    if debug:
        sys.stderr.write(
            "Preparing %i paired FASTQ files\n"
            % len(fastq_file_pairs))

    for stem, left, right in fastq_file_pairs:
        folder, stem = os.path.split(stem)
        if out_dir and out_dir != "-":
            folder = out_dir
        fasta_name = os.path.join(
            folder, "%s.prepared.fasta" % stem)
        if os.path.isfile(fasta_name):
            sys.stderr.write(
                "WARNING: Skipping %s as already exist\n" % fasta_name)
            continue
        sys.stderr.write("Preparing %s\n" % stem)
    return 0
