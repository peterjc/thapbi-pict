"""Assess classification of ITS1 reads.

This implements the ``thapbi_pict assess ...`` command.
"""

import sys

from collections import Counter

from .utils import find_paired_files
from .utils import parse_species_tsv
from .utils import parse_species_list_from_tsv


def sp_in_tsv(classifier_file):
    """Return semi-colon separated list of species in column 2."""
    species = set()
    for line in open(classifier_file):
        if line.startswith("#"):
            continue
        name, taxid, sp, etc = line.split("\t", 3)
        if sp:
            species.update(sp.split(";"))
    return ";".join(sorted(species))


def tally_files(expected_file, predicted_file, min_abundance=0):
    """Make dictionary tally confusion matrix of species assignements."""
    counter = Counter()
    # Sorting because currently not all the classifiers produce out in
    # the same order. The identify classifier respects the input FASTA
    # order (which is by decreasing abundance), while the swarm classifier
    # uses the cluster order. Currently the outputs are all small, so fine.
    for expt, pred in zip(
        sorted(parse_species_tsv(expected_file, min_abundance)),
        sorted(parse_species_tsv(predicted_file, min_abundance)),
    ):
        if not expt[0] == pred[0]:
            sys.exit(
                "Sequence name mismatch in %s vs %s, %s vs %s\n"
                % (expected_file, predicted_file, expt[0], pred[0])
            )
        # TODO: Look at taxid?
        # Should now have (possibly empty) string of genus-species;...
        counter[expt[2], pred[2]] += 1
    return counter


def class_list_from_tally_and_db_list(tally, db_sp_list):
    """Sorted list of all class names used in a confusion table dict."""
    classes = set()
    for expt, pred in tally:
        if expt:
            for sp in expt.split(";"):
                classes.add(sp)
                if sp not in db_sp_list:
                    sys.stderr.write(
                        "WARNING: Expected species %s was not a possible prediction.\n"
                        % sp
                    )
        if pred:
            for sp in pred.split(";"):
                classes.add(sp)
                if sp and sp not in db_sp_list:
                    sys.exit(
                        "ERROR: Species %s was not in the prediction file's header!"
                        % sp
                    )
    classes.update(db_sp_list)
    return sorted(classes)


def save_mapping(tally, filename, debug=False):
    """Output tally table of expected species to predicted sp."""
    with open(filename, "w") as handle:
        handle.write("#Count\tExpected\tPredicted\n")
        for expt, pred in sorted(tally):
            handle.write("%i\t%s\t%s\n" % (tally[expt, pred], expt, pred))
    if debug:
        sys.stderr.write(
            "DEBUG: Wrote %i entry mapping table (total %i) to %s\n"
            % (len(tally), sum(tally.values()), filename)
        )


def save_confusion_matrix(tally, db_sp_list, sp_list, filename, exp_total, debug=False):
    """Output a multi-class confusion matrix as a tab-separated table."""
    total = 0

    assert "" not in db_sp_list
    assert "" not in sp_list
    for sp in db_sp_list:
        assert sp in sp_list

    # leaving out empty cols where can't predict expt (thus dp_sp_list over sp_list)
    # leaving out possible predictions where column would be all zeros
    predicted = set()
    for _, pred in tally:
        if pred:
            for sp in pred.split(";"):
                assert sp in db_sp_list
                predicted.add(sp)
    cols = ["(TN)", "(FN)"] + sorted(predicted)
    del predicted

    # Will report one row per possible combination of expected species
    # None entry (if expected), then single sp (A, B, C), then multi-species (A;B;C etc)
    rows = (
        list(set("(None)" for (expt, pred) in tally if expt == ""))
        + sorted(set(expt for (expt, pred) in tally if expt in sp_list))
        + sorted(set(expt for (expt, pred) in tally if expt and expt not in sp_list))
    )

    assert len(sp_list) * sum(tally.values()) == exp_total

    values = Counter()
    for (expt, pred), count in tally.items():
        if expt:
            assert expt in rows
            e_list = expt.split(";")
        else:
            expt = "(None)"
            e_list = []
        if pred:
            p_list = pred.split(";")
        else:
            p_list = []
        for sp in sp_list:
            if sp in p_list:
                # The TP and FP
                values[expt, sp] += count
            elif sp in e_list:
                values[expt, "(FN)"] += count
            else:
                values[expt, "(TN)"] += count
    total = sum(values.values())

    with open(filename, "w") as handle:
        handle.write("#Expected vs predicted\t%s\n" % "\t".join(cols))
        for expt in rows:
            handle.write(
                "%s\t%s\n" % (expt, "\t".join(str(values[expt, pred]) for pred in cols))
            )

    if debug:
        sys.stderr.write(
            "DEBUG: Wrote %i x %i confusion matrix (total %i) to %s\n"
            % (len(rows), len(cols), total, filename)
        )
    assert total >= sum(tally.values())
    if total != exp_total:
        sys.exit(
            "ERROR: Expected %i but confusion matrix total %i" % (exp_total, total)
        )


def species_level(prediction):
    """Is this prediction at species level.

    Returns True for a binomial name (at least one space), false for genus
    only or no prediction.
    """
    return prediction and " " in prediction


def extract_binary_tally(class_name, tally):
    """Extact single-class TP, FP, FN, TN from multi-class confusion tally.

    Reduces the mutli-class expectation/prediction to binary - did they
    include the class of interest, or not?

    Returns a 4-tuple of values, True Positives (TP), False Positves (FP),
    False Negatives (FN), True Negatives (TN), which sum to the tally total.
    """
    bt = Counter()
    for (expt, pred), count in tally.items():
        bt[class_name in expt.split(";"), class_name in pred.split(";")] += count
    return bt[True, True], bt[False, True], bt[True, False], bt[False, False]


def extract_global_tally(tally, sp_list):
    """Process multi-label confusion matrix (tally dict) to TP, FP, FN, TN.

    If the input data has no negative controls, all there will be no
    true negatives (TN).

    Returns a 4-tuple of values, True Positives (TP), False Positves (FP),
    False Negatives (FN), True Negatives (TN).

    These values are analagous to the classical binary classifier approach,
    but are NOT the same. Even if applied to single class expected and
    predicted values, results differ:

    - Expect none, predict none - 1xTN
    - Expect none, predict A - 1xFP
    - Expect A, predict none - 1xFN
    - Expect A, predict A - 1xTP
    - Expect A, predict B - 1xFP (the B), 1xFN (missing A)
    - Expect A, predict A&B - 1xTP (the A), 1xFP (the B)
    - Expect A&B, predict A&B - 2xTP
    - Expect A&B, predict A - 1xTP, 1xFN (missing B)
    - Expect A&B, predict A&C - 1xTP (the A), 1xFP (the C), 1xFN (missing B)

    The TP, FP, FN, TN sum will exceed the tally total.  For each tally
    entry, rather than one of TP, FP, FN, TN being incremented (weighted
    by the tally count), several can be increased.

    If the input data has no negative controls, all there will be no TN.
    """
    all_sp = set()
    for (expt, pred) in tally:
        if expt:
            all_sp.update(expt.split(";"))
        if pred:
            all_sp.update(pred.split(";"))
    assert "" not in all_sp
    for sp in sp_list:
        assert sp in sp_list, sp

    x = Counter()
    for (expt, pred), count in tally.items():
        if expt == "" and pred == "":
            # TN special case
            x[False, False] += count * len(sp_list)
        else:
            for class_name in sp_list:
                x[class_name in expt.split(";"), class_name in pred.split(";")] += count

    return x[True, True], x[False, True], x[True, False], x[False, False]


assert extract_global_tally({("", ""): 1}, ["A"]) == (0, 0, 0, 1)
assert extract_global_tally({("", ""): 1}, ["A", "B", "C", "D"]) == (0, 0, 0, 4)
assert extract_global_tally({("", "A"): 1}, ["A"]) == (0, 1, 0, 0)
assert extract_global_tally({("", "A"): 1}, ["A", "B", "C", "D"]) == (0, 1, 0, 3)
assert extract_global_tally({("A", ""): 1}, ["A"]) == (0, 0, 1, 0)
assert extract_global_tally({("A", "A"): 1}, ["A"]) == (1, 0, 0, 0)
assert extract_global_tally({("A", "A"): 1}, ["A", "B", "C", "D"]) == (1, 0, 0, 3)
assert extract_global_tally({("A", "B"): 1}, ["A", "B"]) == (0, 1, 1, 0)
assert extract_global_tally({("A", "A;B"): 1}, ["A", "B"]) == (1, 1, 0, 0)
assert extract_global_tally({("A;B", "A;B"): 1}, ["A", "B"]) == (2, 0, 0, 0)
assert extract_global_tally({("A;B", "A"): 1}, ["A", "B"]) == (1, 0, 1, 0)
assert extract_global_tally({("A;B", "A;C"): 1}, ["A", "B", "C"]) == (1, 1, 1, 0)
assert extract_global_tally({("A;B", "A;C"): 1}, ["A", "B", "C", "D"]) == (1, 1, 1, 1)


def main(
    inputs,
    level,
    known,
    method,
    min_abundance,
    assess_output,
    map_output,
    confusion_output,
    debug=False,
):
    """Implement the thapbi_pict assess command."""
    assert isinstance(inputs, list)
    assert level in ["sequence", "sample"], level

    input_list = find_paired_files(
        inputs, ".%s.tsv" % method, ".%s.tsv" % known, debug=False
    )

    db_sp_list = None
    file_count = 0
    global_tally = Counter()

    for predicted_file, expected_file in input_list:
        if debug:
            sys.stderr.write("Assessing %s vs %s\n" % (predicted_file, expected_file))
        if db_sp_list is None:
            db_sp_list = parse_species_list_from_tsv(predicted_file)
        elif db_sp_list != parse_species_list_from_tsv(predicted_file):
            sys.exit("ERROR: Inconsistent species lists in predicted file headers")
        if debug:
            assert db_sp_list is not None, db_sp_list
            sys.stderr.write(
                "DEBUG: %s says DB had %i species\n" % (predicted_file, len(db_sp_list))
            )

        file_count += 1
        if level == "sequence":
            global_tally.update(
                tally_files(expected_file, predicted_file, min_abundance)
            )
        else:
            global_tally[sp_in_tsv(expected_file), sp_in_tsv(predicted_file)] += 1

    if db_sp_list is None:
        sys.exit("ERROR: Failed to load DB species list from headers")
    if debug and db_sp_list:
        sys.stderr.write(
            "Classifier DB had %i species: %s\n"
            % (len(db_sp_list), ", ".join(db_sp_list))
        )
    sp_list = class_list_from_tally_and_db_list(global_tally, db_sp_list)
    if "" in sp_list:
        sp_list.remove("")
    assert sp_list
    if debug:
        sys.stderr.write(
            "Classifier DB had %i species, including expected values have %i species\n"
            % (len(db_sp_list), len(sp_list))
        )
    number_of_classes_and_examples = len(sp_list) * sum(global_tally.values())

    sys.stderr.write(
        "Assessed %s vs %s in %i files (%i species; %i %s level predictions)\n"
        % (method, known, file_count, len(sp_list), sum(global_tally.values()), level)
    )

    assert file_count == len(input_list)

    if not file_count:
        sys.exit("ERROR: Could not find files to assess\n")

    if map_output == "-":
        save_mapping(global_tally, "/dev/stdout", debug=debug)
    elif map_output:
        save_mapping(global_tally, map_output, debug=debug)

    if confusion_output == "-":
        save_confusion_matrix(
            global_tally,
            db_sp_list,
            sp_list,
            "/dev/stdout",
            number_of_classes_and_examples,
            debug=debug,
        )
    elif confusion_output:
        save_confusion_matrix(
            global_tally,
            db_sp_list,
            sp_list,
            confusion_output,
            number_of_classes_and_examples,
            debug=debug,
        )

    if assess_output == "-":
        if debug:
            sys.stderr.write("DEBUG: Output to stdout...\n")
        handle = sys.stdout
    else:
        handle = open(assess_output, "w")

    handle.write(
        "#Species\tTP\tFP\tFN\tTN\t"
        "sensitivity\tspecificity\tprecision\tF1\tHamming-loss\n"
    )
    multi_class_total1 = multi_class_total2 = 0
    for sp in [None] + sp_list:
        if sp is None:
            # Special case flag to report global values at end
            tp, fp, fn, tn = extract_global_tally(global_tally, sp_list)
            sp = "OVERALL"
            multi_class_total1 = tp + fp + fn + tn
        else:
            assert species_level(sp)
            tp, fp, fn, tn = extract_binary_tally(sp, global_tally)
            multi_class_total2 += tp + fp + fn + tn
        # sensitivity, recall, hit rate, or true positive rate (TPR):
        sensitivity = float(tp) / (tp + fn) if tp else 0.0
        # specificity, selectivity or true negative rate (TNR)
        specificity = float(tn) / (tn + fp) if tn else 0.0
        # precision or positive predictive value (PPV)
        precision = float(tp) / (tp + fp) if tp else 0.0
        # F1 score
        f1 = tp * 2.0 / (2 * tp + fp + fn) if tp else 0.0
        # Hamming Loss = (total number of mis-predicted class entries
        #                 / number of class-level predictions)
        hamming_loss = float(fp + fn) / number_of_classes_and_examples
        handle.write(
            "%s\t%i\t%i\t%i\t%i\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.4f\n"
            % (
                sp,
                tp,
                fp,
                fn,
                tn,
                sensitivity,
                specificity,
                precision,
                f1,
                hamming_loss,
            )
        )

    if multi_class_total1 != multi_class_total2 and multi_class_total2:
        sys.exit(
            "ERROR: Overall TP+FP+FN+TP = %i, but sum for species was %i\n"
            % (multi_class_total1, multi_class_total2)
        )

    if assess_output != "-":
        handle.close()

    try:
        sys.stdout.flush()
    except BrokenPipeError:
        pass
    try:
        sys.stderr.flush()
    except BrokenPipeError:
        pass

    return 0
