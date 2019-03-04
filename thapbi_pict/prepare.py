"""Prepare raw ITS1 sequencing reads (trimming, merging, etc).

This implementes the ``thapbi_pict prepare-reads ...`` command.
"""

import os
import subprocess
import shutil
import sys
import tempfile

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from .hmm import filter_for_ITS1
from .utils import abundance_from_read_name
from .utils import abundance_values_in_fasta
from .utils import md5seq
from .utils import run


hmm_cropping_warning = 0  # global variable for warning msg


def find_fastq_pairs(filenames_or_folders, ext=(".fastq", ".fastq.gz"), debug=False):
    """Interpret a list of filenames and/or foldernames.

    Returns a list of tuples (stem, left filename, right filename)
    where stem is intended for use in logging and output naming,
    and may include a directory name.

    The filenames will be normalised relative to the current directory
    (so that we can directly compare file lists which could have been
    defined inconsistently by the user).
    """
    answer = []
    for x in filenames_or_folders:
        x = os.path.normpath(os.path.relpath(x))
        if os.path.isdir(x):
            if debug:
                sys.stderr.write("Walking directory %r\n" % x)
            for root, _, files in os.walk(x):
                for f in files:
                    if f.endswith(ext):
                        answer.append(os.path.join(root, f))
        elif os.path.isfile(x):
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
                sys.stderr.write("WARNING: Ignoring %r and %r\n" % (left, right))
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
    sys.exit("Could not find %s installed with trimmomatic." % fasta_name)


def run_trimmomatic(
    left_in, right_in, left_out, right_out, adapters=None, debug=False, cpu=0
):
    """Run trimmomatic on a pair of FASTQ files.

    The input files may be gzipped.

    If the FASTA adapters file is not specified, will try to use
    TruSeq3-PE.fa as bundled with a BioConda install of trimmomatic.
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


def run_pear(trimmed_R1, trimmed_R2, output_prefix, debug=False, cpu=0):
    """Run pear on a pair of trimmed FASTQ files."""
    cmd = ["pear", "-f", trimmed_R1, "-r", trimmed_R2, "-o", output_prefix]
    if cpu:
        cmd += ["--threads", str(cpu)]
    return run(cmd, debug=debug)


def save_nr_fasta(counts, output_fasta, min_abundance=0):
    r"""Save a dictionary of sequences and counts as a FASTA file.

    The output FASTA records are named ">MD5_abundance\n", which is the
    default style used in SWARM. This could in future be generalised,
    for example ">MD5;size=abundance;\n" for the VSEARCH default.

    Results are sorted by decreasing abundance then alphabetically by
    sequence.
    """
    accepted = 0
    values = sorted((-count, seq) for seq, count in counts.items())
    with open(output_fasta, "w") as out_handle:
        for count, seq in values:
            if -count < min_abundance:
                # Sorted, so everything hereafter is too rare
                break
            out_handle.write(">%s_%i\n%s\n" % (md5seq(seq), -count, seq))
            accepted += 1
    return len(counts), accepted  # number of unique seqs, accepted


def make_nr_fastq_to_fasta(input_fastq, output_fasta, min_abundance=0):
    """Make non-redundant FASTA file from FASTQ inputs, named MD5_abundance.

    The FASTQ read names are ignored and treated as abundance one!

    Returns the number of unique sequences (integer), and the number
    of those which passed the minimum abundance threshold.
    """
    counts = dict()  # OrderedDict on older Python?
    with open(input_fastq) as handle:
        for _, seq, _ in FastqGeneralIterator(handle):
            seq = seq.upper()
            try:
                counts[seq] += 1
            except KeyError:
                counts[seq] = 1
    return save_nr_fasta(counts, output_fasta, min_abundance)


def make_nr_its1(input_fasta, output_fasta, stem, min_abundance=0, debug=False):
    """Make non-redundant FASTA of ITS1 regions, named MD5_abundance.

    Applies HMM with hmmscan to identify any ITS1 region in the
    input reads, applies trimming, then makes these into a
    non-redundant output FASTA file names using a checksum.

    Assumes the input reads also follow this naming, and takes
    the infered abundance into account.

    This naming convention is suitable for SWARM.

    Returns the number of unique ITS1 sequences (integer),
    and the number of those which passed the minimum abundance
    threshold, and the maximum abundance of any one read (for
    use with controls for setting the threshold).
    """
    exp_left = 53
    exp_right = 20
    margin = 10
    file_hmm_cropping_warning = 0
    global hmm_cropping_warning
    # This could be generalised if need something else, e.g.
    # >name;size=6; for VSEARCH.
    counts = dict()  # OrderedDict on older Python?
    for title, full_seq, hmm_seq in filter_for_ITS1(input_fasta, debug=debug):
        if not hmm_seq:
            # Using HMM match as a presense/absense filter
            continue
        seq = full_seq[exp_left:-exp_right].upper()  # using fix trimming
        abundance = abundance_from_read_name(title.split(None, 1)[0])
        left = full_seq.index(hmm_seq)
        right = len(full_seq) - left - len(hmm_seq)
        if not (
            exp_left - margin < left < exp_left + margin
            and exp_right - margin < right < exp_right + margin
        ):
            file_hmm_cropping_warning += 1
            if debug:
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
                sys.stderr.write("Full:  %s (len %i)\n" % (full_seq, len(full_seq)))
                sys.stderr.write(
                    "HMM:   %s%s%s (len %i)\n"
                    % ("-" * left, hmm_seq, "-" * right, len(hmm_seq))
                )
                sys.stderr.write(
                    "Fixed: %s%s%s (len %i)\n"
                    % ("-" * exp_left, seq, "-" * exp_right, len(seq))
                )
        try:
            counts[seq] += abundance
        except KeyError:
            counts[seq] = abundance

    if file_hmm_cropping_warning:
        sys.stderr.write(
            "WARNING: HMM cropping very different from fixed trimming "
            "in %i sequences in %s\n" % (file_hmm_cropping_warning, stem)
        )
        hmm_cropping_warning += file_hmm_cropping_warning

    a, b = save_nr_fasta(counts, output_fasta, min_abundance)
    return a, b, max(counts.values())


def main(fastq, controls, out_dir, min_abundance=100, tmp_dir=None, debug=False, cpu=0):
    """Implement the thapbi_pict prepare-reads command.

    If there are controls, they will be used to potentially increase
    the minimum abundance threshold used for the non-control files.
    """
    global hmm_cropping_warning
    hmm_cropping_warning = 0  # reset for each job

    assert isinstance(fastq, list)

    # Will keep control_min_abundance fixed,
    # will increase min_abundance using max from the controls
    control_min_abundance = min_abundance
    sample_min_abundance = min_abundance

    if not controls:
        control_file_pairs = []
    else:
        control_file_pairs = find_fastq_pairs(controls, debug=debug)

    fastq_file_pairs = find_fastq_pairs(fastq, debug=debug)
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

    if tmp_dir:
        # Up to the user to remove the files
        tmp_obj = None
        shared_tmp = tmp_dir
    else:
        tmp_obj = tempfile.TemporaryDirectory()
        shared_tmp = tmp_obj.name

    if debug:
        sys.stderr.write("DEBUG: Shared temp folder %s\n" % shared_tmp)

    for control, stem, raw_R1, raw_R2 in file_pairs:
        sys.stdout.flush()
        sys.stderr.flush()

        folder, stem = os.path.split(stem)
        if out_dir and out_dir != "-":
            folder = out_dir
        fasta_name = os.path.join(folder, "%s.fasta" % stem)

        if os.path.isfile(fasta_name):
            if control:
                (uniq_count, max_indiv_abundance) = abundance_values_in_fasta(
                    fasta_name
                )
                # TODO - Refactor this duplicated logging?
                sys.stderr.write(
                    "Control %s had %i unique ITS1 sequences "
                    "over control abundance threshold %i (max abundance %i), "
                    % (stem, uniq_count, control_min_abundance, max_indiv_abundance)
                )
                if min_abundance < max_indiv_abundance:
                    sys.stderr.write(
                        "increasing sample abundance threshold from %i\n"
                        % sample_min_abundance
                    )
                else:
                    sys.stderr.write(
                        "keeping sample abundance threshold at %i\n"
                        % sample_min_abundance
                    )
                sample_min_abundance = max(sample_min_abundance, max_indiv_abundance)
                continue
            else:
                sys.stderr.write(
                    "WARNING: Skipping %s as already exists\n" % fasta_name
                )
                continue

        min_abundance = control_min_abundance if control else sample_min_abundance

        sys.stderr.write(
            "Starting to prepare %s %s (min abundance set to %i)\n"
            % ("control" if control else "sample", fasta_name, min_abundance)
        )
        sys.stdout.flush()
        sys.stderr.flush()

        tmp = os.path.join(shared_tmp, stem)
        if not os.path.isdir(tmp):
            # If using tempfile.TemporaryDirectory() for shared_tmp
            # this will be deleted automatically, otherwise user must:
            os.mkdir(tmp)

        if True:
            if debug:
                sys.stderr.write("DEBUG: Temp folder of %s is %s\n" % (stem, tmp))

            # trimmomatic
            trim_R1 = os.path.join(tmp, "trimmomatic_R1.fastq")
            trim_R2 = os.path.join(tmp, "trimmomatic_R2.fastq")
            run_trimmomatic(raw_R1, raw_R2, trim_R1, trim_R2, debug=debug, cpu=cpu)
            for _ in (trim_R1, trim_R2):
                if not os.path.isfile(_):
                    sys.exit("ERROR: Expected file %r from trimmomatic\n" % _)

            # pear
            pear_prefix = os.path.join(tmp, "pear")
            merged_fastq = os.path.join(tmp, "pear.assembled.fastq")
            run_pear(trim_R1, trim_R2, pear_prefix, debug=debug, cpu=cpu)
            if not os.path.isfile(merged_fastq):
                sys.exit("ERROR: Expected file %r from pear\n" % merged_fastq)

            merged_fasta = os.path.join(tmp, "dedup_long.fasta")
            # Do not apply min_abundance threshold here as after ITS1
            # trimming pooling entries would increase their counts.
            count, _ = make_nr_fastq_to_fasta(
                merged_fastq, merged_fasta, min_abundance=0
            )
            if debug:
                sys.stderr.write(
                    "Merged paired FASTQ reads into %i unique sequences\n" % count
                )

            # Find the ITS1 region (if present) using hmmscan,
            # make NR, and name as MD5_abundance
            # Apply the min_abundance threshold here (at the final step)
            dedup = os.path.join(tmp, "dedup_its1.fasta")
            uniq_count, acc_uniq_count, max_indiv_abundance = make_nr_its1(
                merged_fasta, dedup, stem, min_abundance, debug
            )

            # File done
            shutil.move(dedup, fasta_name)

            # Update threshold and logging
            if control:
                sys.stderr.write(
                    "Control %s has %i unique ITS1 sequences "
                    "over control abundance threshold %i (max abundance %i), "
                    % (stem, acc_uniq_count, control_min_abundance, max_indiv_abundance)
                )
                if sample_min_abundance < max_indiv_abundance:
                    sys.stderr.write(
                        "increasing sample abundance threshold from %i\n"
                        % sample_min_abundance
                    )
                else:
                    sys.stderr.write(
                        "keeping sample abundance threshold at %i\n"
                        % sample_min_abundance
                    )
                sample_min_abundance = max(sample_min_abundance, max_indiv_abundance)
            else:
                sys.stderr.write(
                    "Wrote %s with %i unique sequences "
                    "over sample abundance theshold %i (max abundance %i)\n"
                    % (stem, acc_uniq_count, sample_min_abundance, max_indiv_abundance)
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
    return 0
