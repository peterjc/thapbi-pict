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


def class_list_from_tally(tally):
    """Sorted list of all class names used in a confusion table dict."""
    classes = set()
    for expt, pred in tally:
        classes.add(expt)
        classes.add(pred)
    return sorted(classes)


def save_confusion_matrix(tally, filename, debug=False):
    """Output a multi-class confusion matrix as a tab-separated table."""
    species = class_list_from_tally(tally)
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


def species_level(prediction):
    """Is this prediction atspecies level.

    Returns True for a binomial name (at least one space), false for genus
    only or no prediction.
    """
    return prediction and " " in prediction


def extract_binary_tally(class_name, tally):
    """Compress a multi-class confusion matrix to a single-species binary one.

    Returns a 4-tuple of values, True Positives (TP), False Positves (FP),
    False Negatives (FN), True Negatives (TN)
    """
    bt = Counter()
    for (expt, pred), count in tally.items():
        bt[expt == class_name, pred == class_name] += count
    return bt[True, True], bt[False, True], bt[True, False], bt[False, False]


def extract_global_tally(tally):
    """Compress multi-class confusion matrix to global binary one.

    Treats no prediction and genus level only prediction as negatives,
    trests species level prediction as positives.

    If the input data has no negative controls, all there will be no
    true negatives.

    Returns a 4-tuple of values, True Positives (TP), False Positves (FP),
    False Negatives (FN), True Negatives (TN).
    """
    tp = fp = fn = tn = 0
    for (expt, pred), count in tally.items():
        if species_level(expt):
            # Have species level expectation...
            if species_level(pred):
                # this is either TP or FP
                if expt == pred:
                    tp += count
                else:
                    fp += count
            else:
                # No species level prediction, FN
                fn += count
        elif species_level(pred):
            # Have no species level expectation, but a (false) species
            # level prediction was made - FP
            fp += count
        else:
            # Have no species level expectation, no species
            # level prediction was made - TN
            tn += count
    return tp, fp, fn, tn


def main(inputs, known, method, assess_output, confusion_output, debug=False):
    """Implement the thapbi_pict assess command."""
    assert isinstance(inputs, list)

    input_list = find_requested_files(inputs, ext=".%s-reads.tsv" % method, debug=debug)

    count = 0
    global_tally = Counter()

    # Context manager should remove the temp dir:
    with tempfile.TemporaryDirectory() as shared_tmp:
        if debug:
            sys.stderr.write("DEBUG: Shared temp folder %s\n" % shared_tmp)
        for predicted_file in input_list:
            stem = predicted_file.rsplit(".", 2)[0]
            assert predicted_file == stem + ".%s-reads.tsv" % method
            expected_file = stem + ".%s-reads.tsv" % known
            # Not aborting here as typically in a folder while should have
            # full set of predictions, will only have expected results for
            # control subset.
            if not os.path.isfile(expected_file):
                sys.stderr.write("WARNING: Missing %s\n" % expected_file)
                continue
            if debug:
                sys.stderr.write("Assessing %s vs %s for %s\n" % (method, known, stem))

            file_tally = tally_files(expected_file, predicted_file)
            count += 1
            # print(file_tally)
            global_tally.update(file_tally)

    sys.stderr.write("Assessed %s vs %s in %i files\n" % (method, known, count))

    if not count:
        sys.exit("ERROR: Could not find files to assess\n")

    if confusion_output == "-":
        save_confusion_matrix(global_tally, "/dev/stdout", debug=debug)
    elif confusion_output:
        save_confusion_matrix(global_tally, confusion_output, debug=debug)

    if assess_output == "-":
        if debug:
            sys.stderr.write("DEBUG: Output to stdout...\n")
        handle = sys.stdout
    else:
        handle = open(assess_output, "w")

    handle.write("#Species\tTP\tFP\tFN\tTN\tsensitivity\tspecificity\tprecision\tF1\n")
    sp_list = [_ for _ in class_list_from_tally(global_tally) if species_level(_)]
    for species in [None] + sp_list:
        if species is None:
            # Special case flag to report global values at end
            tp, fp, fn, tn = extract_global_tally(global_tally)
            species = "OVERALL"
        else:
            tp, fp, fn, tn = extract_binary_tally(species, global_tally)
        # sensitivity, recall, hit rate, or true positive rate (TPR):
        sensitivity = float(tp) / (tp + fn) if tp else 0.0
        # specificity, selectivity or true negative rate (TNR)
        specificity = float(tn) / (tn + fp) if tn else 0.0
        # precision or positive predictive value (PPV)
        precision = float(tp) / (tp + fp) if tp else 0.0
        # F1 score
        f1 = tp * 2.0 / (2 * tp + fp + fn) if tp else 0.0
        handle.write(
            "%s\t%i\t%i\t%i\t%i\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n"
            % (species, tp, fp, fn, tn, sensitivity, specificity, precision, f1)
        )
    if assess_output != "-":
        handle.close()

    sys.stdout.flush()
    sys.stderr.flush()

    return 0
