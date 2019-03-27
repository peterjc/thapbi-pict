"""Summarise ITS1 classification results at sample (multi-plate) level.

This implements the ``thapbi_pict sample-summary ...`` command.
"""

import os
import sys

from collections import Counter

from .utils import find_requested_files
from .utils import parse_species_tsv
from .utils import abundance_from_read_name


def main(inputs, output, human_output, method, min_abundance=1, debug=False):
    """Implement the thapbi_pict sample-summary command.

    The expectation is that the inputs represent all the samples from
    a meaningful group, likely from multiple sequencing runs (plates).
    """
    assert isinstance(inputs, list)

    if not (output or human_output):
        sys.exit("No output file specified.\n")

    samples = set()
    counts = Counter()
    sp_to_taxid = {}
    tsv_files = find_requested_files(inputs, ".%s.tsv" % method, debug)
    if debug:
        sys.stderr.write(
            "Loading %i sample predictions using method %s\n" % (len(tsv_files), method)
        )
    for predicted_file in tsv_files:
        sample = os.path.basename(predicted_file).rsplit(".", 2)[0]
        samples.add(sample)  # Check for duplicates (aside from methods)
        for name, taxid_list, sp_list in parse_species_tsv(
            predicted_file, min_abundance
        ):
            taxid_list = taxid_list.split(";")
            sp_list = sp_list.split(";")
            assert len(taxid_list) == len(sp_list), predicted_file
            unambig = len(sp_list) == 1
            for sp, taxid in zip(sp_list, taxid_list):
                if sp in sp_to_taxid:
                    assert sp_to_taxid[sp] == taxid, "Clash for taxid %s" % taxid
                else:
                    sp_to_taxid[sp] = taxid
                counts[sample, sp, unambig] += abundance_from_read_name(name)

    if debug:
        sys.stderr.write("Loaded predictions for %i samples\n" % len(samples))

    if output == "-":
        if debug:
            sys.stderr.write("DEBUG: Output to stdout...\n")
        handle = sys.stdout
    elif output:
        handle = open(output, "w")
    else:
        handle = None

    if human_output == "-":
        if debug:
            sys.stderr.write("DEBUG: Output human report to stdout...\n")
        human = sys.stdout
    elif human_output:
        human = open(human_output, "w")
    else:
        human = None

    if handle:
        handle.write("#Sample\tTaxID\tSpecies\tUnambiguous\tSeq-count\n")
    for sample in samples:
        all_sp = set()
        unambig_sp = set()
        for sp in sp_to_taxid:
            for unambig in [True, False]:
                count = counts[sample, sp, unambig]
                if count:
                    if handle:
                        handle.write(
                            "%s\t%s\t%s\t%s\t%i\n"
                            % (sample, sp_to_taxid[sp], sp, unambig, count)
                        )
                    all_sp.add(sp)
                    if unambig:
                        unambig_sp.add(sp)
        if not all_sp and handle:
            # Match unclassified count output: TaxID zero, blank species, True
            handle.write("%s\t%s\t%s\t%s\t%i\n" % (sample, "0", "", True, 0))
        if human:
            human.write("%s\n\n" % sample)
            for sp in sorted(all_sp):
                if sp not in unambig_sp:
                    sp = "(%s)" % sp
                if not sp:
                    sp = "Unknown"
                human.write(" - %s\n" % sp)
            if not all_sp:
                human.write(" - No data\n")
            human.write("\n")

    if output != "-" and handle:
        handle.close()
    if human_output != "-" and human:
        human.close()

    try:
        sys.stdout.flush()
    except BrokenPipeError:
        pass
    try:
        sys.stderr.flush()
    except BrokenPipeError:
        pass

    return 0
