# Copyright 2018-2020 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Assess classification of ITS1 reads at species level.

This implements the ``thapbi_pict assess ...`` command.
"""
import sys
from collections import Counter

from .utils import abundance_from_read_name
from .utils import find_paired_files
from .utils import parse_species_list_from_tsv
from .utils import parse_species_tsv
from .utils import species_level
from .utils import split_read_name_abundance


def sp_in_tsv(classifier_file, min_abundance):
    """Return semi-colon separated list of species in column 2.

    Will ignore genus level predictions.
    """
    species = set()
    for line in open(classifier_file):
        if line.startswith("#") or not line.strip():
            continue
        try:
            name, taxid, sp, etc = line.split("\t", 3)
        except ValueError:
            sys.exit(f"ERROR: Wrong field count in {classifier_file}\n")
        if not sp:
            continue
        # Using "*" as name (without abundance) as wildcard
        if (
            name != "*"
            and min_abundance > 1
            and abundance_from_read_name(name) < min_abundance
        ):
            continue
        species.update(sp.split(";"))
    return ";".join(sorted(_ for _ in species if species_level(_)))


def tally_files(expected_file, predicted_file, min_abundance=0):
    """Make dictionary tally confusion matrix of species assignements.

    Rather than the values simply being an integer count, they are
    the set of MD5 identifiers (take the length for the count).
    """
    counter = {}
    # Sorting because currently not all the classifiers produce out in
    # the same order. The identify classifier respects the input FASTA
    # order (which is by decreasing abundance), while the swarm classifier
    # uses the cluster order. Currently the outputs are all small, so fine.
    try:
        for expt, pred in zip(
            sorted(parse_species_tsv(expected_file, req_species_level=True)),
            sorted(parse_species_tsv(predicted_file, req_species_level=True)),
        ):
            if not expt[0] == pred[0]:
                sys.exit(
                    f"ERROR: Sequence name mismatch in {expected_file} vs"
                    f" {predicted_file}, {expt[0]} vs {pred[0]}\n"
                )
            md5, abundance = split_read_name_abundance(expt[0])
            if min_abundance > 1 and abundance < min_abundance:
                continue
            # TODO: Look at taxid?
            # Should now have (possibly empty) string of genus-species;...
            for sp in expt[2].split(";"):
                if sp:
                    assert species_level(sp), (
                        f"Expectation {expt[2]} is not all species level from"
                        f" {expected_file}"
                    )
            for sp in pred[2].split(";"):
                if sp:
                    assert species_level(sp), (
                        f"Prediction {pred[2]} is not all species level from"
                        f" {predicted_file}"
                    )
            try:
                counter[expt[2], pred[2]].add(md5)
            except KeyError:
                counter[expt[2], pred[2]] = {md5}
    except ValueError as e:
        # This is used for single sample controls,
        # where all reads are expected to be from species X.
        if str(e) != "Wildcard species name found":
            raise
        expt_sp_genus = None
        with open(expected_file) as handle:
            for line in handle:
                if line.startswith("*\t"):
                    _, _, expt_sp_genus, _ = line.split("\t", 3)
        assert expt_sp_genus, "Didn't find expected wildcard species line"
        for pred in parse_species_tsv(
            predicted_file, min_abundance, req_species_level=True
        ):
            md5, abundance = split_read_name_abundance(pred[0])
            if min_abundance > 1 and abundance < min_abundance:
                continue
            try:
                counter[expt_sp_genus, pred[2]].add(md5)
            except KeyError:
                counter[expt_sp_genus, pred[2]] = {md5}
    return counter


def class_list_from_tally_and_db_list(tally, db_sp_list):
    """Sorted list of all class names used in a confusion table dict."""
    classes = set()
    impossible = set()
    for expt, pred in tally:
        if expt:
            for sp in expt.split(";"):
                assert species_level(
                    sp
                ), f"{sp} from expected value {pred} is not species-level"
                classes.add(sp)
                if sp not in db_sp_list:
                    impossible.add(sp)
        if pred:
            for sp in pred.split(";"):
                assert species_level(
                    sp
                ), f"{sp} from prediction {pred} is not species-level"
                classes.add(sp)
                if sp and sp not in db_sp_list:
                    sys.exit(
                        f"ERROR: Species {sp} was not in the prediction file's header!"
                    )
    if impossible:
        sys.stderr.write(
            f"WARNING: {len(impossible)} expected species were not a"
            f" possible prediction: {';'.join(sorted(impossible))}\n"
        )

    classes.update(db_sp_list)
    return sorted(classes)


def save_mapping(tally, filename, debug=False):
    """Output tally table of expected species to predicted sp."""
    with open(filename, "w") as handle:
        handle.write("#sample-count\tExpected\tPredicted\n")
        for expt, pred in sorted(tally):
            handle.write(f"{tally[expt, pred]}\t{expt}\t{pred}\n")
    if debug:
        sys.stderr.write(
            f"DEBUG: Wrote {len(tally)} entry mapping table"
            f" (total {sum(tally.values())}) to {filename}\n"
        )


def save_confusion_matrix(tally, db_sp_list, sp_list, filename, exp_total, debug=False):
    """Output a multi-class confusion matrix as a tab-separated table."""
    total = 0

    assert "" not in db_sp_list
    assert "" not in sp_list
    for sp in db_sp_list:
        assert sp in sp_list, sp
    for sp in sp_list:
        assert species_level(sp), sp

    # Prediction lines will cover TP and FP
    predicted = set()
    for _, pred in tally:
        if pred:
            for sp in pred.split(";"):
                assert sp in db_sp_list, sp
                predicted.add(sp)
    cols = ["(TN)", "(FN)"] + sorted(predicted)
    del predicted

    # Will report one row per possible combination of expected species
    # None entry (if expected), then single sp (A, B, C), then multi-species (A;B;C etc)
    rows = (
        list({"(None)" for (expt, pred) in tally if expt == ""})
        + sorted({expt for (expt, pred) in tally if expt in sp_list})
        + sorted({expt for (expt, pred) in tally if expt and expt not in sp_list})
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
        handle.write("#Expected vs predicted\tsample count\t%s\n" % "\t".join(cols))
        for expt in rows:
            if expt == "(None)":
                level_count = sum(count for ((e, _), count) in tally.items() if not e)
            else:
                level_count = sum(
                    count for ((e, _), count) in tally.items() if e == expt
                )
            handle.write(
                "%s\t%i\t%s\n"
                % (
                    expt,
                    level_count,
                    "\t".join(str(values[expt, pred]) for pred in cols),
                )
            )

    if debug:
        sys.stderr.write(
            f"DEBUG: Wrote {len(rows)} x {len(cols)} confusion matrix"
            f" (total {total}) to {filename}\n"
        )
    assert total >= sum(tally.values())
    if total != exp_total:
        sys.exit(f"ERROR: Expected {exp_total} but confusion matrix total {total}")


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
    for sp in all_sp:
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
    known,
    method,
    min_abundance,
    assess_output,
    map_output,
    confusion_output,
    ignore_prefixes=None,
    debug=False,
):
    """Implement the (sample/species level) ``thapbi_pict assess`` command."""
    assert isinstance(inputs, list), inputs

    input_list = find_paired_files(
        inputs, f".{method}.tsv", f".{known}.tsv", ignore_prefixes, debug, strict=False
    )

    if not input_list:
        sys.exit(f"ERROR: Found no pairs *.{method}.tsv with *.{known}.tsv")

    db_sp_list = None
    file_count = 0
    global_tally = {}

    for predicted_file, expected_file in input_list:
        if debug:
            sys.stderr.write(f"Assessing {predicted_file} vs {expected_file}\n")
        if db_sp_list is None:
            # Note this will ignore genus level classifications:
            db_sp_list = parse_species_list_from_tsv(predicted_file)
        elif db_sp_list != parse_species_list_from_tsv(predicted_file):
            sys.exit("ERROR: Inconsistent species lists in predicted file headers")
        if debug:
            assert db_sp_list is not None, db_sp_list

        file_count += 1
        expt = sp_in_tsv(expected_file, min_abundance)
        pred = sp_in_tsv(predicted_file, min_abundance)
        global_tally[expt, pred] = global_tally.get((expt, pred), 0) + 1
    assert len(input_list) == sum(global_tally.values())

    if debug:
        sys.stderr.write(f"DEBUG: Assessing {sum(global_tally.values())} predictions\n")

    if db_sp_list is None:
        sys.exit("ERROR: Failed to load DB species list from headers")
    for sp in db_sp_list:
        assert species_level(sp), sp
    sp_list = class_list_from_tally_and_db_list(global_tally, db_sp_list)
    if "" in sp_list:
        sp_list.remove("")
    assert sp_list
    if debug:
        sys.stderr.write(
            f"Classifier DB had {len(db_sp_list)} species,"
            f" including expected values have {len(sp_list)} species\n"
        )
    for sp in sp_list:
        assert species_level(sp), sp

    number_of_classes_and_examples = len(sp_list) * sum(global_tally.values())

    sys.stderr.write(
        f"Assessed {method} vs {known} in {file_count} files"
        f" ({len(sp_list)} species)\n"
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
        "sensitivity\tspecificity\tprecision\tF1\t"
        "Hamming-loss\tAd-hoc-loss\n"
    )
    multi_class_total1 = multi_class_total2 = 0
    boring_species_count = 0
    boring_species_tn = 0
    for sp in [""] + sp_list + [None]:
        if sp == "":
            # Special case flag to report global values at top
            tp, fp, fn, tn = extract_global_tally(global_tally, sp_list)
            sp = "OVERALL"
            multi_class_total1 = tp + fp + fn + tn
        elif sp is None:
            # Special case flag to report boring sp summary at end
            sp = f"OTHER {boring_species_count} SPECIES IN DB"
            tp = fp = fn = 0
            tn = boring_species_tn
            if not boring_species_count:
                assert boring_species_tn == 0
                # No point printing this line
                continue
        else:
            assert species_level(sp)
            tp, fp, fn, tn = extract_binary_tally(sp, global_tally)
            multi_class_total2 += tp + fp + fn + tn
            if tp == fp == fn:
                # Boring output line, just DB padding
                boring_species_tn += tn
                boring_species_count += 1
                continue
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
        try:
            hamming_loss = float(fp + fn) / (tp + fp + fn + tn)
        except ZeroDivisionError:
            hamming_loss = 0.0
        # Ad-hoc Loss = (total number of mis-predicted class entries
        #                 / number of class-level predictions ignoring TN
        try:
            ad_hoc_loss = float(fp + fn) / (tp + fp + fn)
        except ZeroDivisionError:
            ad_hoc_loss = 0.0
        handle.write(
            "%s\t%i\t%i\t%i\t%i\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.4f\t%0.3f\n"
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
                ad_hoc_loss,
            )
        )

    if multi_class_total1 != number_of_classes_and_examples:
        sys.exit(
            f"ERROR: Overall TP+FP+FN+TP = {multi_class_total1},"
            f" but species times samples = {number_of_classes_and_examples}\n"
        )
    if multi_class_total1 != multi_class_total2 and multi_class_total2:
        sys.exit(
            f"ERROR: Overall TP+FP+FN+TP = {multi_class_total1},"
            f" but sum for species was {multi_class_total2}\n"
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
