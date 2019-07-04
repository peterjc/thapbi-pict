# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

"""Summarise ITS1 classification results at sample (multi-plate) level.

This implements the ``thapbi_pict sample-summary ...`` command.
"""

import os
import sys

from collections import Counter

from .utils import find_requested_files
from .utils import load_metadata
from .utils import parse_species_tsv
from .utils import abundance_from_read_name
from .utils import sample_sort


def main(
    inputs,
    output,
    human_output,
    method,
    min_abundance=1,
    metadata_file=None,
    metadata_cols=None,
    metadata_fieldnames=None,
    metadata_index=None,
    debug=False,
):
    """Implement the thapbi_pict sample-summary command.

    The expectation is that the inputs represent all the samples from
    a meaningful group, likely from multiple sequencing runs (plates).
    """
    assert isinstance(inputs, list)

    if not (output or human_output):
        sys.exit("ERROR: No output file specified.\n")

    samples = set()
    counts = Counter()
    sp_to_taxid = {}
    tsv_files = find_requested_files(inputs, ".%s.tsv" % method, debug)
    if debug:
        sys.stderr.write(
            "Loading %i sample predictions using method %s\n" % (len(tsv_files), method)
        )
    if not tsv_files:
        sys.exit("ERROR: No input files found\n")
    for predicted_file in tsv_files:
        sample = os.path.basename(predicted_file).rsplit(".", 2)[0]
        if sample in samples:
            sys.exit("ERROR: Duplicate sample name: %s\n" % sample)
        samples.add(sample)
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

    samples = sample_sort(samples)
    (
        metadata_rows,
        metadata_samples,
        meta_names,
        meta_default,
        missing_meta,
    ) = load_metadata(
        metadata_file,
        metadata_cols,
        metadata_fieldnames,
        metadata_index,
        sequenced_samples=samples,
        metadata_sort=True,
        debug=debug,
    )

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
    if human:
        human.write(
            "NOTE: Species listed with (uncertain/ambiguous) in brackets are where "
            "sequences matched multiple species equally well. For example, "
            "Phytophthora andina, P. infestans, and P. ipomoeae, share an identical "
            "marker.\n\n"
        )

    # Note already sorted on metadata values, discarded the order in the table
    batches = list(zip(metadata_rows, metadata_samples))
    if missing_meta:
        batches.append([meta_default, missing_meta])
    for metadata, sample_batch in batches:
        if human and meta_names:
            # Write the metadata header
            try:
                human.write("-" * 60 + "\n\n")
                if metadata:
                    for name, value in zip(meta_names, metadata):
                        if value:
                            human.write("%s: %s\n" % (name, value))
                    human.write("\n")
                else:
                    human.write("Missing metadata\n\n")
                if not sample_batch:
                    human.write("Has not been sequenced.\n\n")
                # elif len(sample_batch) == 1:
                #     human.write("Has been sequenced once:\n\n")
                # else:
                #     human.write("Has been sequenced %i times:\n\n")
            except BrokenPipeError:
                human = None
        # Now do the samples in this batch
        for sample in sample_batch:
            if sample not in samples:
                sys.stderr.write("WARNING: Missing %s\n" % sample)
            all_sp = set()
            unambig_sp = set()
            for sp in sp_to_taxid:
                for unambig in [True, False]:
                    count = counts[sample, sp, unambig]
                    if count:
                        if handle:
                            try:
                                handle.write(
                                    "%s\t%s\t%s\t%s\t%i\n"
                                    % (sample, sp_to_taxid[sp], sp, unambig, count)
                                )
                            except BrokenPipeError:
                                # Stop trying to write to stdout (eg piped to head)
                                handle = None
                        all_sp.add(sp)
                        if unambig:
                            unambig_sp.add(sp)
            if not all_sp and handle:
                try:
                    # Match unclassified count output: TaxID zero, blank species, True
                    handle.write("%s\t%s\t%s\t%s\t%i\n" % (sample, "0", "", True, 0))
                except BrokenPipeError:
                    handle = None

            if human:
                try:
                    if meta_names:
                        human.write("Sequencing sample: %s\n\n" % sample)
                    else:
                        human.write("%s\n\n" % sample)
                    for sp in sorted(all_sp):
                        if sp not in unambig_sp:
                            sp = "%s (uncertain/ambiguous)" % sp
                        if not sp:
                            sp = "Unknown"
                        human.write(" - %s\n" % sp)
                    if not all_sp:
                        human.write(" - No data\n")
                    human.write("\n")
                except BrokenPipeError:
                    # Stop trying to write to stdout (e.g. piped to head)
                    human = None

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
