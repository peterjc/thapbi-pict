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
from .utils import expand_IUPAC_ambiguity_codes
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

    Returns the number of sequences accepted (above any minimum
    abundance specified).
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
    return accepted


def make_nr_fastq_to_fasta(
    input_fastq, output_fasta, left_primer, right_primer, min_abundance=0, debug=False
):
    """Trim and make non-redundant FASTA file from FASTQ inputs.

    The FASTQ read names are ignored and treated as abundance one!

    Applies the specified fixed triming, then makes a non-redundant
    FASTA file with the sequences named MD5_abundance.

    Returns the number of sequences before de-duplication (integer),
    number of unique sequences (integer), and the number of those
    which passed the minimum abundance threshold, and the maximum
    abundance of any one read (for use with controls for setting
    the threshold).
    """
    bad_start = 0
    bad_end = 0
    trim_left = len(left_primer)
    trim_right = len(right_primer)
    trim_starts = tuple(expand_IUPAC_ambiguity_codes(left_primer))
    trim_ends = tuple(expand_IUPAC_ambiguity_codes(right_primer))
    counts = dict()  # OrderedDict on older Python?
    with open(input_fastq) as handle:
        for _, seq, _ in FastqGeneralIterator(handle):
            seq = seq.upper()
            if not seq.startswith(trim_starts):
                bad_start += 1
            if not seq.endswith(trim_ends):
                bad_end += 1
            seq = seq[trim_left:-trim_right]
            try:
                counts[seq] += 1
            except KeyError:
                counts[seq] = 1
    if debug:
        sys.stderr.write(
            "DEBUG: %i unique from %i merged FASTQ reads, "
            "%i and %i didn't have expected primer-based start and end\n"
            % (len(counts), sum(counts.values()), bad_start, bad_end)
        )
    return (
        sum(counts.values()),
        len(counts),
        save_nr_fasta(counts, output_fasta, min_abundance),
        max(counts.values()),
    )


def filter_fasta_for_its1(input_fasta, output_fasta, stem, debug=False):
    """Filter for ITS1 regions.

    Assumes you have already applied trimming (and are not using
    the HMM results for this).

    Assumes the SWARM naming convention.

    Returns the number of unique ITS1 sequences (integer),
    """
    margin = 10
    cropping_warning = 0
    max_indiv_abundance = 0
    # This could be generalised if need something else, e.g.
    # >name;size=6; for VSEARCH.
    count = 0
    with open(output_fasta, "w") as out_handle:
        for title, full_seq, hmm_seq in filter_for_ITS1(input_fasta, debug=debug):
            if not hmm_seq:
                # Using HMM match as a presense/absense filter
                continue
            out_handle.write(">%s\n%s\n" % (title, full_seq))
            count += 1
            max_indiv_abundance = max(
                max_indiv_abundance, abundance_from_read_name(title.split(None, 1)[0])
            )

            left = full_seq.index(hmm_seq)
            right = len(full_seq) - left - len(hmm_seq)
            if margin < left or margin < right:
                cropping_warning += 1
                if debug:
                    sys.stderr.write(
                        "WARNING: %s has HMM cropping %i left, %i right "
                        "(on top of fixed trimming)\n"
                        % (title.split(None, 1)[0], left, right)
                    )

    if cropping_warning:
        sys.stderr.write(
            "WARNING: HMM cropping very different from fixed trimming "
            "in %i sequences in %s\n" % (cropping_warning, stem)
        )

    return count, max_indiv_abundance, cropping_warning


def main(
    fastq,
    controls,
    out_dir,
    left_primer,
    right_primer,
    min_abundance=100,
    debug=False,
    cpu=0,
):
    """Implement the thapbi_pict prepare-reads command.

    If there are controls, they will be used to potentially increase
    the minimum abundance threshold used for the non-control files.
    """
    global hmm_cropping_warning
    hmm_cropping_warning = 0  # reset for each job

    assert isinstance(fastq, list)

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
                    "Control %s had %i unique ITS1 sequences, "
                    "%i of most abundant, " % (stem, uniq_count, max_indiv_abundance)
                )
                if min_abundance < max_indiv_abundance:
                    sys.stderr.write(
                        "increasing abundance threshold from %i\n" % min_abundance
                    )
                else:
                    sys.stderr.write(
                        "keeping abundance threshold at %i\n" % min_abundance
                    )
                min_abundance = max(min_abundance, max_indiv_abundance)
                continue
            else:
                sys.stderr.write(
                    "WARNING: Skipping %s as already exists\n" % fasta_name
                )
                continue

        sys.stderr.write("Starting to prepare %s\n" % fasta_name)
        sys.stdout.flush()
        sys.stderr.flush()

        # Context manager should remove the temp dir:
        with tempfile.TemporaryDirectory() as tmp:
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

            merged_fasta = os.path.join(tmp, "dedup_trimmed.fasta")
            (
                count,
                uniq_count,
                acc_uniq_count,
                max_indiv_abundance,
            ) = make_nr_fastq_to_fasta(
                merged_fastq,
                merged_fasta,
                left_primer,
                right_primer,
                min_abundance=0 if control else min_abundance,
                debug=debug,
            )
            if control:
                assert uniq_count == acc_uniq_count
                if debug:
                    sys.stderr.write(
                        "Merged %i paired FASTQ reads into %i unique sequences "
                        "(max abundance %i)\n"
                        % (count, uniq_count, acc_uniq_count, max_indiv_abundance)
                    )
            else:
                if debug:
                    sys.stderr.write(
                        "Merged %i paired FASTQ reads into %i unique sequences, "
                        "%i above min abundance %i (max abundance %i)\n"
                        % (
                            count,
                            uniq_count,
                            acc_uniq_count,
                            min_abundance,
                            max_indiv_abundance,
                        )
                    )

            # Find the ITS1 region (if present) using hmmscan,
            dedup = os.path.join(tmp, "dedup_its1.fasta")
            uniq_count, max_indiv_abundance, cropping = filter_fasta_for_its1(
                merged_fasta, dedup, stem, debug=debug
            )
            hmm_cropping_warning += cropping
            if control:
                sys.stderr.write(
                    "Control %s has %i unique ITS1 sequences, "
                    "%i of most abundant, " % (stem, uniq_count, max_indiv_abundance)
                )
                if min_abundance < max_indiv_abundance:
                    sys.stderr.write(
                        "increasing abundance threshold from %i\n" % min_abundance
                    )
                else:
                    sys.stderr.write(
                        "keeping abundance threshold at %i\n" % min_abundance
                    )
                min_abundance = max(min_abundance, max_indiv_abundance)
            elif debug:
                sys.stderr.write(
                    "Filtered %s down to %i unique ITS1 sequences, "
                    "passing abundance threshold %i, "
                    "with top abundance %i\n"
                    % (stem, uniq_count, min_abundance, max_indiv_abundance)
                )

            # File done
            shutil.move(dedup, fasta_name)
            if control:
                sys.stderr.write(
                    "Wrote %s with %i unique control reads\n" % (stem, uniq_count)
                )
            else:
                sys.stderr.write(
                    "Wrote %s with %i unique reads over abundance %i\n"
                    % (stem, uniq_count, min_abundance)
                )

    if hmm_cropping_warning:
        sys.stderr.write(
            "WARNING: HMM cropping very different from fixed trimming "
            "in %i sequences over %i paired FASTQ files\n"
            % (hmm_cropping_warning, len(file_pairs))
        )

    sys.stdout.flush()
    sys.stderr.flush()
    return 0
