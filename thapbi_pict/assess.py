"""Assess classification of ITS1 reads.

This implements the ``thapbi_pict assess ...`` command.
"""

import os
import sys
import tempfile

from collections import Counter

from .utils import find_requested_files


def parse_species_tsv(tabular_file):
    """Parse file of species assignments/predictions by sequence."""
    with open(tabular_file) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            name, genus, species, etc = line.split("\t", 3)
            yield name, genus, species


def tally_files(expected_file, predicted_file):
    """Make dictionary tally confusion matrix of species assignements."""
    counter = Counter()
    # Sorting because currently not all the classifiers produce out in
    # the same order. The identify classifier respects the input FASTA
    # order (which is by decreasing abundance), while the swarm classifier
    # uses the cluster order. Currently the outputs are all small, so fine.
    for expt, pred in zip(
        sorted(parse_species_tsv(expected_file)),
        sorted(parse_species_tsv(predicted_file)),
    ):
        if not expt[0] == pred[0]:
            sys.exit(
                "Sequence name mismatch in %s vs %s, %s vs %s\n"
                % (expected_file, predicted_file, expt[0], pred[0])
            )
        # Might only have genus with species "", thus strip whitespace:
        expt_sp = ("%s %s" % (expt[1], expt[2])).strip()
        pred_sp = ("%s %s" % (pred[1], pred[2])).strip()
        counter[expt_sp, pred_sp] += 1
    return counter


def save_confusion_matrix(tally, filename, debug=False):
    """Output a multi-class confusion matrix as a tab-separated table."""
    species = set()
    for expt, pred in tally:
        species.add(expt)
        species.add(pred)
    species = sorted(species)
    with open(filename, "w") as handle:
        handle.write("#\t%s\n" % "\t".join(species))
        for pred in species:
            handle.write(
                "%s\t%s\n"
                % (pred, "\t".join(str(tally[pred, expt]) for expt in species))
            )
    if debug:
        sys.stderr.write(
            "DEBUG: Wrote %i x %i confusion matrix to %s\n"
            % (len(species), len(species), filename)
        )


def main(fasta, known, method, out_dir, debug=False):
    """Implement the thapbi_pict assess command."""
    assert isinstance(fasta, list)

    fasta_list = find_requested_files(fasta, ext=".fasta", debug=debug)

    count = 0
    global_tally = Counter()

    # Context manager should remove the temp dir:
    with tempfile.TemporaryDirectory() as shared_tmp:
        if debug:
            sys.stderr.write("DEBUG: Shared temp folder %s\n" % shared_tmp)
        for f in fasta_list:
            sys.stderr.write("Assessing %s vs %s for %s\n" % (method, known, f))
            stem = f[:-6]
            predicted_file = stem + ".%s-reads.tsv" % method
            expected_file = stem + ".%s-reads.tsv" % known
            # Not aborting here as typically in a folder while should have
            # full set of predictions, will only have expected results for
            # control subset.
            if not os.path.isfile(predicted_file):
                sys.stderr.write("WARNING: Missing %s\n" % predicted_file)
                continue
            expected_file = stem + ".%s-reads.tsv" % known
            if not os.path.isfile(expected_file):
                sys.stderr.write("WARNING: Missing %s\n" % expected_file)
                continue

            file_tally = tally_files(expected_file, predicted_file)
            count += 1
            # print(file_tally)
            global_tally.update(file_tally)

    sys.stderr.write("Assessed %s vs %s in %i files\n" % (method, known, count))
    save_confusion_matrix(global_tally, "/dev/stdout", debug=debug)

    sys.stdout.flush()
    sys.stderr.flush()

    if not count:
        sys.exit("ERROR: Could not find files to assess\n")

    return 0
