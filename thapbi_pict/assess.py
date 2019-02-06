"""Assess classification of ITS1 reads.

This implements the ``thapbi_pict assess ...`` command.
"""

import sys
import tempfile

from collections import Counter

from .utils import find_paired_files


def parse_species_tsv(tabular_file):
    """Parse file of species assignments/predictions by sequence."""
    with open(tabular_file) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            name, taxid, genus, species, etc = line.split("\t", 4)
            yield name, taxid, genus, species


def untangle_species(taxid, genus, species):
    """Untangle classifier predictions which might have ; entries.

    Returns semi-colon separated string of species names (in the
    binomial form, genus-species).

    Does not currently use the taxid.
    """
    if ";" in taxid:
        assert ";" in species, "%s %s %s" % (taxid, genus, species)
    if ";" in species:
        assert ";" in taxid or taxid == 0, "%s %s %s" % (taxid, genus, species)
    if ";" in genus:
        assert ";" in taxid or taxid == 0, "%s %s %s" % (taxid, genus, species)

    if not species:
        return ""  # No species level predictions

    answer = set()
    if ";" in genus:
        assert species.count(";") == genus.count(";")
        for s, g in zip(species.split(";"), genus.split(";")):
            assert s.strip()
            answer.add("%s %s" % (g, s))
    else:
        for s in species.split(";"):
            assert s.strip()
            answer.add("%s %s" % (genus, s))
    return ";".join(sorted(answer))


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
        # TODO: Look at taxid?
        # This has dealt with ambiguous entries with ; separated entries
        # Should now have (possibly empty) string of genus-species;...
        expt_sp = untangle_species(*expt[1:])
        pred_sp = untangle_species(*pred[1:])
        counter[expt_sp, pred_sp] += 1
    return counter


def class_list_from_tally(tally):
    """Sorted list of all class names used in a confusion table dict."""
    classes = set()
    for expt, pred in tally:
        for sp in expt.split(";"):
            classes.add(sp)
        for sp in pred.split(";"):
            classes.add(sp)
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


def save_confusion_matrix(tally, filename, debug=False):
    """Output a multi-class confusion matrix as a tab-separated table."""
    total = 0
    species = class_list_from_tally(tally)
    extra_rows = sorted(set(expt for (expt, pred) in tally if expt not in species))
    with open(filename, "w") as handle:
        handle.write("#\t%s\n" % "\t".join(species))
        for expt in species:
            values = Counter()
            for (e, p), count in tally.items():
                if expt == e:
                    for pred in p.split(";"):
                        values[pred] += count
                    # if expt not in p.split(";"):
                    #     values[""] += count
            assert values[""] >= tally.get((expt, ""), 0), values
            handle.write(
                "%s\t%s\n" % (expt, "\t".join(str(values[pred]) for pred in species))
            )
            total += sum(values.values())
        for expt in extra_rows:
            # These are the multi-species expected entries
            values = Counter()
            for (e, p), count in tally.items():
                if expt == e:
                    for pred in p.split(";"):
                        values[pred] += count
            assert values[""] >= tally.get((expt, ""), 0), values
            handle.write(
                "%s\t%s\n" % (expt, "\t".join(str(values[pred]) for pred in species))
            )
            total += sum(values.values())
    if debug:
        sys.stderr.write(
            "DEBUG: Wrote %i x %i confusion matrix (total %i) to %s\n"
            % (len(species) + len(extra_rows), len(species), total, filename)
        )
    assert total >= sum(tally.values())


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


def extract_global_tally(tally):
    """Process multi-label confusion matrix (tally dict) to TP, FP, DN, TN.

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
    tp = fp = fn = tn = 0
    for (expt, pred), count in tally.items():
        expt_sp_list = expt.split(";") if expt else []
        pred_sp_list = pred.split(";") if pred else []

        if expt_sp_list:
            # Hopefully some TP...
            if pred_sp_list:
                # Have some combination of TP, FP, FN
                for sp in expt_sp_list:
                    if sp in pred_sp_list:
                        tp += count
                    else:
                        fn += count
                for sp in pred_sp_list:
                    if sp not in expt_sp_list:
                        fp += count
            else:
                # No predictions, these all FN
                fn += count * len(expt_sp_list)
        else:
            # Have no species level expectation,
            if pred_sp_list:
                # False predictions made - FP
                fp += count * len(pred_sp_list)
            else:
                # No species level prediction was made - TN
                tn += count

    return tp, fp, fn, tn


assert extract_global_tally({("", ""): 1}) == (0, 0, 0, 1)
assert extract_global_tally({("", "A"): 1}) == (0, 1, 0, 0)
assert extract_global_tally({("A", ""): 1}) == (0, 0, 1, 0)
assert extract_global_tally({("A", "A"): 1}) == (1, 0, 0, 0)
assert extract_global_tally({("A", "B"): 1}) == (0, 1, 1, 0)
assert extract_global_tally({("A", "A;B"): 1}) == (1, 1, 0, 0)
assert extract_global_tally({("A;B", "A;B"): 1}) == (2, 0, 0, 0)
assert extract_global_tally({("A;B", "A"): 1}) == (1, 0, 1, 0)
assert extract_global_tally({("A;B", "A;C"): 1}) == (1, 1, 1, 0)


def main(
    inputs, known, method, assess_output, map_output, confusion_output, debug=False
):
    """Implement the thapbi_pict assess command."""
    assert isinstance(inputs, list)

    input_list = find_paired_files(
        inputs, ".%s.tsv" % method, ".%s.tsv" % known, debug=False
    )

    file_count = 0
    global_tally = Counter()

    # Context manager should remove the temp dir:
    with tempfile.TemporaryDirectory() as shared_tmp:
        if debug:
            sys.stderr.write("DEBUG: Shared temp folder %s\n" % shared_tmp)
        for predicted_file, expected_file in input_list:
            if debug:
                sys.stderr.write(
                    "Assessing %s vs %s\n" % (predicted_file, expected_file)
                )
            file_count += 1
            global_tally.update(tally_files(expected_file, predicted_file))

    sys.stderr.write(
        "Assessed %s vs %s in %i files (%i sequence entries)\n"
        % (method, known, file_count, sum(global_tally.values()))
    )

    assert file_count == len(input_list)

    if not file_count:
        sys.exit("ERROR: Could not find files to assess\n")

    if map_output == "-":
        save_mapping(global_tally, "/dev/stdout", debug=debug)
    elif map_output:
        save_mapping(global_tally, map_output, debug=debug)

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
    sp_list = class_list_from_tally(global_tally)
    # Ensure "" special case is at the start of the list (will be first if sorted)
    if "" in sp_list:
        sp_list.remove("")
    sp_list = [""] + sp_list
    for sp in sp_list:
        if sp:
            assert species_level(sp)
    for sp in sp_list:
        if not sp:
            # Special case flag to report global values at end
            tp, fp, fn, tn = extract_global_tally(global_tally)
            sp = "OVERALL"
        else:
            tp, fp, fn, tn = extract_binary_tally(sp, global_tally)
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
            % (sp, tp, fp, fn, tn, sensitivity, specificity, precision, f1)
        )
    if assess_output != "-":
        handle.close()

    sys.stdout.flush()
    sys.stderr.flush()

    return 0
