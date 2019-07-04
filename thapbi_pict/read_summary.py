# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

"""Summarise ITS1 classification results at folder (plate) level.

This implements the ``thapbi_pict plate-summary ...`` command.
"""

import os
import sys

from collections import Counter

from Bio.SeqIO.FastaIO import SimpleFastaParser

from .utils import find_requested_files
from .utils import load_metadata
from .utils import parse_species_tsv
from .utils import sample_sort
from .utils import split_read_name_abundance


def color_bands(meta_groups, sample_color_bands, debug=False):
    """Return a list for formats, one for each sample."""
    default = [None] * len(meta_groups)

    min_groups = 3
    max_groups = 0.8 * len(meta_groups)
    if min_groups > max_groups:
        if debug:
            sys.stderr.write("DEBUG: Not enough columns to bother with coloring\n")
        return default

    if len(set(meta_groups)) == 1:
        # All the same not helpful for banding
        if debug:
            sys.stderr.write(
                "DEBUG: All samples had same metadata in color grouping field: %r\n"
                % meta_groups[0]
            )
        return default

    if max_groups < len(set(meta_groups)):
        if debug:
            sys.stderr.write(
                "DEBUG: Too many coloring groups, trying first word only\n"
            )
        # Attempting heuristic, taking the first word/field will work on schemes like
        # SITE_DATE_NUMBER or SPECIES-SAMPLE etc.
        meta_groups = [
            _.replace("-", " ").replace("_", " ").split(None, 1)[0] if _ else ""
            for _ in meta_groups
        ]
        if len(set(meta_groups)) == 1:
            # That didn't work.
            if debug:
                sys.stderr.write(
                    "DEBUG: Too many coloring groups, but first word only was unique\n"
                )
            return default

    if len(set(meta_groups)) < min_groups or max_groups < len(set(meta_groups)):
        # (Almost) all same or (almost) unique not helpful
        if debug:
            sys.stderr.write(
                "DEBUG: %i groups not suitable for coloring\n" % len(set(meta_groups))
            )
        return default

    bands = []
    debug_msg = []
    for s, value in enumerate(meta_groups):
        if s == 0:
            bands.append(0)
            debug_msg.append(value)
        elif value == meta_groups[s - 1]:
            # Same
            bands.append(bands[-1])
        else:
            # Different
            if value in meta_groups[:s]:
                # Metadata values are not grouped, can't use for banding
                # (might be able to use with a color key?)
                sys.stderr.write(
                    "WARNING: Metadata not grouped nicely for coloring: "
                    "%s and then %s (again).\n" % (", ".join(debug_msg), value)
                )
                # if debug:
                #     sys.stderr.write("DEBUG: %r\n" % meta_groups)
                return default
            bands.append(max(bands) + 1)
            debug_msg.append(value)
    assert len(set(bands)) == max(bands) + 1, bands
    return [sample_color_bands[_ % len(sample_color_bands)] for _ in bands]


def main(
    inputs,
    output,
    method,
    min_abundance=1,
    excel=None,
    metadata_file=None,
    metadata_cols=None,
    metadata_groups=None,
    metadata_fieldnames=None,
    metadata_index=None,
    debug=False,
):
    """Implement the thapbi_pict plate-summary command.

    The expectation is that the inputs represent all the samples
    from one (96 well) plate, or some other meaningful batch.
    """
    assert isinstance(inputs, list)

    if not output:
        sys.exit("ERROR: No output file specified.\n")

    # Which column are we using for group coloring?
    if metadata_groups and not metadata_cols:
        sys.exit("ERROR: Using -g / --metagroups requires -c / --metacols")
    if metadata_cols:
        # Tiny code duplication from load_metadata:
        try:
            value_cols = [int(_) - 1 for _ in metadata_cols.split(",")]
        except ValueError:
            sys.exit(
                "ERROR: Output metadata columns should be a comma separated list "
                "of positive integers, not %r." % metadata_cols
            )
    if metadata_groups:
        try:
            group_col = int(metadata_groups) - 1
        except ValueError:
            sys.exit(
                "ERROR: Invalid metadata group column, should be positive or 0, not %r."
                % metadata_groups
            )
        if group_col not in value_cols:
            sys.exit(
                "ERROR: Metadata group column %i not included in reported metadata.\n"
                % metadata_groups
            )
        group_col = value_cols.index(group_col)  # i.e. which of requested columns
        del value_cols
    else:
        group_col = 0  # Default is the first requested column

    samples = set()
    md5_abundance = Counter()
    abundance_by_samples = {}
    md5_species = {}
    md5_to_seq = {}

    if debug:
        sys.stderr.write("Loading FASTA sequences and abundances\n")
    for fasta_file in find_requested_files(
        [_ for _ in inputs if not _.endswith(".tsv")], ".fasta", debug
    ):
        sample = os.path.basename(fasta_file).rsplit(".", 1)[0]
        samples.add(sample)
        with open(fasta_file) as handle:
            for title, seq in SimpleFastaParser(handle):
                md5, abundance = split_read_name_abundance(title.split(None, 1)[0])
                if min_abundance > 1 and abundance < min_abundance:
                    continue
                abundance_by_samples[md5, sample] = abundance
                md5_abundance[md5] += abundance
                md5_to_seq[md5] = seq
                md5_species[md5] = set()
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
    # Turn row-centric metadata into a dictionary keyed on sequenced sample name
    metadata = {}
    new = []
    # Already sorted rows of the metadata values, discarded the order in the table
    for row, r_samples in zip(metadata_rows, metadata_samples):
        for sample in r_samples:
            if sample in samples:
                # print(sample, row)
                metadata[sample] = row
                assert sample not in new, sample
                new.append(sample)
    for sample in missing_meta:
        # print("Missing metadata for %s" % sample)
        assert sample not in new, sample
        new.append(sample)
    assert set(samples) == set(new)
    assert len(samples) == len(new)
    samples = new
    del metadata_rows, metadata_samples, new

    methods = method.split(",")
    for method in methods:
        if debug:
            sys.stderr.write("Loading predictions for %s\n" % method)
        tsv_files = find_requested_files(
            [_ for _ in inputs if not _.endswith(".fasta")], ".%s.tsv" % method, debug
        )
        if len(samples) != len(tsv_files):
            sys.exit(
                "ERROR: Identified %i samples from FASTA files, but %i TSV files for %s"
                % (len(samples), len(tsv_files), method)
            )
        for predicted_file in tsv_files:
            sample = os.path.basename(predicted_file).rsplit(".", 2)[0]
            assert sample in samples, predicted_file
            # TODO: Look at taxid here?
            for name, _, sp in parse_species_tsv(predicted_file, min_abundance):
                md5, abundance = split_read_name_abundance(name)
                if min_abundance > 1 and abundance < min_abundance:
                    continue
                assert abundance_by_samples[md5, sample] == abundance, name
                # Combining over all methods!
                if sp:
                    md5_species[md5].update(sp.split(";"))

    if output == "-":
        if debug:
            sys.stderr.write("DEBUG: Output to stdout...\n")
        handle = sys.stdout
    else:
        handle = open(output, "w")

    if excel:
        import xlsxwriter

        workbook = xlsxwriter.Workbook(excel)
        worksheet = workbook.add_worksheet("Sequence vs samples")
        cell_rightalign_format = workbook.add_format({"align": "right"})
        red_conditional_format = workbook.add_format(
            # Maraschino red
            {"bg_color": "#FF2600", "font_color": "#000000"}
        )
        sample_color_bands = [
            # Simple rolling rainbow pastel pallet
            workbook.add_format({"bg_color": c, "font_color": "#000000"})
            for c in [
                "#FFCCDA",  # pink
                "#F7D6B7",  # orange
                "#FFFFCC",  # yellow
                "#CCFFDD",  # green
                "#CCF7FF",  # blue
            ]
        ]

        # If there are lots of samples, set narrow column widths
        if len(samples) > 50:
            # Set column width to 2
            worksheet.set_column(5, 5 + len(samples), 2)
        elif len(samples) > 20:
            # Set column width to 4
            worksheet.set_column(5, 5 + len(samples), 4)
    else:
        workbook = None
        worksheet = None

    current_row = 0
    first_data_row = 0
    sample_formats = [None] * len(samples)
    if metadata:
        # Insert extra header rows at start for sample meta-data
        # Make a single metadata call for each sample
        meta = [metadata.get(sample, meta_default) for sample in samples]
        for i, name in enumerate(meta_names):
            handle.write("#\t\t\t\t%s\t%s\n" % (name, "\t".join(_[i] for _ in meta)))
        if worksheet:
            sample_formats = color_bands(
                [_[group_col] for _ in meta], sample_color_bands, debug=debug
            )
            for i, name in enumerate(meta_names):
                worksheet.write_string(i, 4, name, cell_rightalign_format)
                for s, sample in enumerate(samples):
                    worksheet.write_string(
                        i,
                        5 + s,
                        metadata.get(sample, meta_default)[i],
                        sample_formats[s],
                    )
            current_row += len(meta_names)
    handle.write(
        "#ITS1-MD5\t%s-predictions\tSequence\tSample-count\tTotal-abundance\t%s\n"
        % (",".join(methods), "\t".join(samples))
    )
    handle.write(
        "TOTAL\t-\t-\t%i\t%i\t%s\n"
        % (
            sum(
                1
                for md5 in md5_to_seq
                for sample in samples
                if (md5, sample) in abundance_by_samples
            ),
            sum(md5_abundance.values()),
            "\t".join(
                str(
                    sum(
                        abundance_by_samples.get((md5, sample), 0) for md5 in md5_to_seq
                    )
                )
                for sample in samples
            ),
        )
    )
    if worksheet:
        worksheet.write_string(current_row, 0, "ITS1-MD5")
        worksheet.write_string(current_row, 1, "%s-predictions" % ",".join(methods))
        worksheet.write_string(current_row, 2, "Sequence")
        worksheet.write_string(current_row, 3, "Sample-count")
        worksheet.write_string(current_row, 4, "Total-abundance")
        for s, sample in enumerate(samples):
            worksheet.write_string(current_row, 5 + s, sample, sample_formats[s])
        current_row += 1
        first_data_row = current_row
        worksheet.write_string(current_row, 0, "TOTAL")
        worksheet.write_string(current_row, 1, "-")
        worksheet.write_string(current_row, 2, "-")
        worksheet.write_number(
            current_row,
            3,
            sum(
                1
                for md5 in md5_to_seq
                for sample in samples
                if (md5, sample) in abundance_by_samples
            ),
        )
        worksheet.write_number(current_row, 4, sum(md5_abundance.values()))
        for s, sample in enumerate(samples):
            worksheet.write_number(
                current_row,
                5 + s,
                sum(abundance_by_samples.get((md5, sample), 0) for md5 in md5_to_seq),
                sample_formats[s],
            )
        current_row += 1
        worksheet.freeze_panes(current_row, 5)  # keep total line in view plus headers

    # Build the first few columns as a list of lists, which we can sort
    data = [
        [
            md5,
            ";".join(sorted(md5_species[md5])),
            md5_to_seq[md5],
            sum(1 for _ in samples if (md5, _) in abundance_by_samples),
            total_abundance,
        ]
        for md5, total_abundance in md5_abundance.items()
    ]
    # Sort on species prediction (with blank last, sorting as tilde);
    # number of samples (decreasing), total abundance (decreasing), md5
    data.sort(key=lambda row: (row[1] if row[1] else "~", -row[3], -row[4], row[0]))
    for md5, sp, seq, md5_in_xxx_samples, total_abundance in data:
        sample_counts = [abundance_by_samples.get((md5, _), 0) for _ in samples]
        handle.write(
            "%s\t%s\t%s\t%i\t%i\t%s\n"
            % (
                md5,
                sp,
                seq,
                md5_in_xxx_samples,
                total_abundance,
                "\t".join(str(_) for _ in sample_counts),
            )
        )
        if worksheet:
            worksheet.write_string(current_row, 0, md5)
            worksheet.write_string(current_row, 1, sp)
            worksheet.write_string(current_row, 2, seq)
            worksheet.write_number(current_row, 3, md5_in_xxx_samples)
            worksheet.write_number(current_row, 4, total_abundance)
            for s, count in enumerate(sample_counts):
                worksheet.write_number(current_row, 5 + s, count, sample_formats[s])
            current_row += 1
    del data

    if worksheet:
        worksheet.conditional_format(
            first_data_row,
            5,
            current_row,
            5 + len(samples),
            {
                "type": "cell",
                "criteria": "greater than",
                "value": 0,
                "format": red_conditional_format,
            },
        )

    if output != "-":
        handle.close()
    if workbook:
        workbook.close()

    try:
        sys.stdout.flush()
    except BrokenPipeError:
        pass
    try:
        sys.stderr.flush()
    except BrokenPipeError:
        pass

    return 0
