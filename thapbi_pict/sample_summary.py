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
    excel,
    human_output,
    method,
    min_abundance=1,
    metadata_file=None,
    metadata_cols=None,
    metadata_fieldnames=None,
    metadata_index=None,
    require_metadata=False,
    ignore_prefixes=None,
    debug=False,
):
    """Implement the ``thapbi_pict sample-summary`` command.

    The expectation is that the inputs represent all the samples from
    a meaningful group, likely from multiple sequencing runs (plates).
    """
    assert isinstance(inputs, list)

    if not (output or human_output):
        sys.exit("ERROR: No output file specified.\n")

    samples = set()
    tsv_files = find_requested_files(inputs, f".{method}.tsv", ignore_prefixes, debug)
    if debug:
        sys.stderr.write(
            f"Loading {len(tsv_files):d} sample predictions using method {method}\n"
        )
    if not tsv_files:
        sys.exit("ERROR: No input files found\n")

    genus_predictions = set()
    sample_genus_counts = {}
    species_predictions = set()  # includes A;B;C ambiguous entries
    sample_species_counts = {}
    for predicted_file in tsv_files:
        sample = os.path.basename(predicted_file).rsplit(".", 2)[0]
        if sample in samples:
            sys.exit(f"ERROR: Duplicate sample name: {sample}\n")
        samples.add(sample)
        sample_species_counts[sample] = Counter()
        sample_genus_counts[sample] = Counter()
        for name, _taxid_list, sp_list in parse_species_tsv(
            predicted_file, min_abundance
        ):
            species_predictions.add(sp_list)  # as string with any ; included
            sample_species_counts[sample][sp_list] += abundance_from_read_name(name)
            genus_list = {sp.split(" ", 1)[0] for sp in sp_list.split(";")}
            genus_predictions.update(genus_list)
            if len(genus_list) > 1:
                sys.stderr.write(
                    f"WARNING: Conflicting genus for {name} from {sample}\n"
                )
            for genus in genus_list:
                sample_genus_counts[sample][genus] += abundance_from_read_name(name)
    species_predictions = sorted(species_predictions)  # turn into a list
    genus_predictions = sorted(genus_predictions)  # turn into a list
    species_columns = [_ for _ in species_predictions if _ not in genus_predictions]

    if debug:
        sys.stderr.write(
            f" {len(samples):d} samples"
            f" with predictions for {len(genus_predictions):d} genera\n"
        )

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
        require_metadata=require_metadata,
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

    if excel:
        import xlsxwriter

        workbook = xlsxwriter.Workbook(excel)
        worksheet = workbook.add_worksheet("Sequence vs samples")
        red_conditional_format = workbook.add_format(
            # Maraschino red
            {"bg_color": "#FF2600", "font_color": "#000000"}
        )
        vertical_text_format = workbook.add_format(
            # Vertical text, reading up the page
            {"rotation": 90}
        )
    else:
        workbook = None
        worksheet = None
    current_row = 0

    if human_output == "-":
        if debug:
            sys.stderr.write("DEBUG: Output human report to stdout...\n")
        human = sys.stdout
    elif human_output:
        human = open(human_output, "w")
    else:
        human = None

    if handle:
        if meta_names:
            handle.write(
                "#%s\tSequencing sample\tSeq-count\t%s\t%s\n"
                % (
                    "\t".join(meta_names),
                    "\t".join(_ if _ else "Unknown" for _ in genus_predictions),
                    "\t".join(species_columns),
                )
            )
        else:
            handle.write(
                "#Sequencing sample\tSeq-count\t%s\t%s\n"
                % (
                    "\t".join(_ if _ else "Unknown" for _ in genus_predictions),
                    "\t".join(species_columns),
                )
            )
    if worksheet:
        col_offset = len(meta_names)
        # Set first row to be tall, with vertical text
        worksheet.set_row(0, 150, vertical_text_format)
        # If there are lots of species, set narrow column widths
        cols = len(genus_predictions) + len(species_columns)
        if cols > 50:
            # Set column width to 2
            worksheet.set_column(col_offset + 1, col_offset + 1 + cols, 2)
        elif cols > 20:
            # Set column width to 4
            worksheet.set_column(col_offset + 1, col_offset + 1 + cols, 4)
        del cols
        for offset, name in enumerate(meta_names):
            worksheet.write_string(current_row, offset, name)
        col_offset = len(meta_names)
        worksheet.write_string(current_row, col_offset, "Sequencing sample")
        worksheet.write_string(current_row, col_offset + 1, "Seq-count")
        for offset, genus in enumerate(genus_predictions):
            worksheet.write_string(
                current_row, col_offset + 2 + offset, genus if genus else "Unknown"
            )
        for offset, sp in enumerate(species_columns):
            worksheet.write_string(
                current_row, col_offset + 2 + len(genus_predictions) + offset, sp
            )
        worksheet.freeze_panes(current_row + 1, col_offset + 2)

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
        assert not require_metadata
        batches.append([meta_default, missing_meta])
    for metadata, sample_batch in batches:
        if human and meta_names:
            # Write the metadata header
            try:
                human.write("-" * 60 + "\n\n")
                if metadata:
                    for name, value in zip(meta_names, metadata):
                        if value:
                            human.write(f"{name}: {value}\n")
                    human.write("\n")
                else:
                    human.write("Missing metadata\n\n")
                    assert not require_metadata
                if not sample_batch:
                    human.write("Has not been sequenced.\n\n")
            except BrokenPipeError:
                human = None
        # Now do the samples in this batch
        for sample in sample_batch:
            if sample not in samples:
                sys.stderr.write(f"WARNING: Missing {sample}\n")
            else:
                if handle:
                    try:
                        if metadata:
                            handle.write("\t".join(metadata) + "\t")
                        handle.write(
                            "%s\t%i\t%s\t%s\n"
                            % (
                                sample,
                                sum(sample_species_counts[sample].values()),
                                "\t".join(
                                    str(sample_genus_counts[sample][genus])
                                    for genus in genus_predictions
                                ),
                                "\t".join(
                                    str(sample_species_counts[sample][sp])
                                    for sp in species_predictions
                                    if sp not in genus_predictions
                                ),
                            )
                        )
                    except BrokenPipeError:
                        # Stop trying to write to stdout (eg piped to head)
                        handle = None
                if worksheet:
                    current_row += 1
                    assert len(meta_names) == len(metadata)
                    for offset, value in enumerate(metadata):
                        worksheet.write_string(current_row, offset, value)
                    assert col_offset == len(meta_names)
                    worksheet.write_string(current_row, col_offset, sample)
                    worksheet.write_number(
                        current_row,
                        col_offset + 1,
                        sum(sample_species_counts[sample].values()),
                    )
                    for offset, genus in enumerate(genus_predictions):
                        worksheet.write_number(
                            current_row,
                            col_offset + 2 + offset,
                            sample_genus_counts[sample][genus],
                        )
                    for offset, sp in enumerate(species_columns):
                        worksheet.write_number(
                            current_row,
                            col_offset + 2 + len(genus_predictions) + offset,
                            sample_species_counts[sample][sp],
                        )

            if human:
                all_sp = set()
                unambig_sp = set()
                for sp_list, count in sample_species_counts[sample].items():
                    if count:
                        sp_list = sp_list.split(";")
                        all_sp.update(sp_list)
                        if len(sp_list) == 1:
                            unambig_sp.add(sp_list[0])
                try:
                    if meta_names:
                        human.write(f"Sequencing sample: {sample}\n\n")
                    else:
                        human.write(f"{sample}\n\n")
                    for sp in sorted(all_sp):
                        if sp not in unambig_sp:
                            sp = f"{sp} (uncertain/ambiguous)"
                        if not sp:
                            sp = "Unknown"
                        human.write(f" - {sp}\n")
                    if not all_sp:
                        human.write(" - No data\n")
                    human.write("\n")
                except BrokenPipeError:
                    # Stop trying to write to stdout (e.g. piped to head)
                    human = None

    if worksheet:
        worksheet.conditional_format(
            1,
            col_offset + 1,
            current_row,
            col_offset + 1 + len(genus_predictions) + len(species_columns),
            {
                "type": "cell",
                "criteria": "greater than",
                "value": 0,
                "format": red_conditional_format,
            },
        )
    if workbook:
        workbook.close()

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
