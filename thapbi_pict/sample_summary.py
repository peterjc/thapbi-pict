"""Summarise ITS1 classification results at sample (multi-plate) level.

This implements the ``thapbi_pict sample-summary ...`` command.
"""

import os
import sys

from collections import Counter

from .utils import find_requested_files
from .utils import parse_species_tsv
from .utils import abundance_from_read_name


def main(inputs, output, method, min_abundance=1, debug=False):
    """Implement the thapbi_pict sample-summary command.

    The expectation is that the inputs represent all the samples from
    a meaningful group, likely from multiple sequencing runs (plates).
    """
    assert isinstance(inputs, list)

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
    else:
        handle = open(output, "w")

    handle.write("#Sample\tTaxID\tSpecies\tUnambiguous\tSeq-count\n")
    for sample in samples:
        output_from_sample = False
        for sp in sp_to_taxid:
            for unambig in [True, False]:
                count = counts[sample, sp, unambig]
                if count:
                    handle.write(
                        "%s\t%s\t%s\t%s\t%i\n"
                        % (sample, sp_to_taxid[sp], sp, unambig, count)
                    )
                    output_from_sample = True
        if not output_from_sample:
            # Match unclassified count output: TaxID zero, blank species, True
            handle.write("%s\t%s\t%s\t%s\t%i\n" % (sample, "0", "", True, 0))

    if output != "-":
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
