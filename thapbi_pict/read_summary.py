# Copyright 2019-2020 by Peter Cock, The James Hutton Institute.
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

from .utils import color_bands
from .utils import find_requested_files
from .utils import load_metadata
from .utils import parse_species_tsv
from .utils import sample_sort
from .utils import split_read_name_abundance


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
    ignore_prefixes=None,
    debug=False,
):
    """Implement the ``thapbi_pict plate-summary`` command.

    The expectation is that the inputs represent all the samples
    from one (96 well) plate, or some other meaningful batch.
    """
    assert isinstance(inputs, list)

    if not output:
        sys.exit("ERROR: No output file specified.\n")

    samples = set()
    md5_abundance = Counter()
    abundance_by_samples = {}
    md5_species = {}
    md5_to_seq = {}

    if debug:
        sys.stderr.write("Loading FASTA sequences and abundances\n")
    for fasta_file in find_requested_files(
        [_ for _ in inputs if not _.endswith(".tsv")],
        ".fasta",
        ignore_prefixes,
        debug=debug,
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
        group_col,
    ) = load_metadata(
        metadata_file,
        metadata_cols,
        metadata_groups,
        metadata_fieldnames,
        metadata_index,
        sequenced_samples=samples,
        metadata_sort=True,
        ignore_prefixes=ignore_prefixes,
        debug=debug,
    )
    # Turn row-centric metadata into a dictionary keyed on sequenced sample name
    metadata = {}
    new = []
    # Already sorted rows of the metadata values, discarded the order in the table
    for row, r_samples in zip(metadata_rows, metadata_samples):
        for sample in r_samples:
            if sample in samples:
                metadata[sample] = row
                assert sample not in new, sample
                new.append(sample)
    for sample in missing_meta:
        assert sample not in new, sample
        new.append(sample)
    assert set(samples) == set(new)
    assert len(samples) == len(new)
    samples = new
    del metadata_rows, metadata_samples, new

    methods = method.split(",")
    for method in methods:
        if debug:
            sys.stderr.write(f"Loading predictions for {method}\n")
        tsv_files = find_requested_files(
            [_ for _ in inputs if not _.endswith(".fasta")],
            f".{method}.tsv",
            ignore_prefixes,
            debug,
        )
        if len(samples) != len(tsv_files):
            sys.exit(
                f"ERROR: Identified {len(samples):d} samples from FASTA files,"
                f" but {len(tsv_files):d} TSV files for {method}"
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

    LEADING_COLS = 6
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
            worksheet.set_column(LEADING_COLS, LEADING_COLS + len(samples), 2)
        elif len(samples) > 20:
            # Set column width to 4
            worksheet.set_column(LEADING_COLS, LEADING_COLS + len(samples), 4)
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
            handle.write(
                "#%s\t%s\t%s\n"
                % ("\t" * LEADING_COLS, name, "\t".join(_[i] for _ in meta))
            )
        if worksheet:
            sample_formats = color_bands(
                [_[group_col] for _ in meta], sample_color_bands, debug=debug
            )
            for i, name in enumerate(meta_names):
                worksheet.write_string(
                    i, LEADING_COLS - 1, name, cell_rightalign_format
                )
                for s, sample in enumerate(samples):
                    worksheet.write_string(
                        i,
                        LEADING_COLS + s,
                        metadata.get(sample, meta_default)[i],
                        sample_formats[s],
                    )
            current_row += len(meta_names)
    handle.write(
        "#ITS1-MD5\t%s-predictions\tSequence\tSample-count"
        "\tMax-sample-abundance\tTotal-abundance\t%s\n"
        % (",".join(methods), "\t".join(samples))
    )
    handle.write(
        "TOTAL\t-\t-\t%i\t%i\t%i\t%s\n"
        % (
            max(
                (
                    abundance_by_samples.get((md5, sample), 0)
                    for md5 in md5_to_seq
                    for sample in samples
                ),
                default=0,
            ),
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
        worksheet.write_string(current_row, 1, ",".join(methods) + "-predictions")
        worksheet.write_string(current_row, 2, "Sequence")
        worksheet.write_string(current_row, 3, "Sample-count")
        worksheet.write_string(current_row, 4, "Max-sample-abundance")
        worksheet.write_string(current_row, 5, "Total-abundance")
        for s, sample in enumerate(samples):
            worksheet.write_string(
                current_row, LEADING_COLS + s, sample, sample_formats[s]
            )
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
        worksheet.write_number(
            current_row,
            4,
            max(
                (
                    abundance_by_samples.get((md5, sample), 0)
                    for md5 in md5_to_seq
                    for sample in samples
                ),
                default=0,
            ),
        )
        worksheet.write_number(current_row, 5, sum(md5_abundance.values()))
        for s, sample in enumerate(samples):
            worksheet.write_number(
                current_row,
                LEADING_COLS + s,
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
            "%s\t%s\t%s\t%i\t%i\t%i\t%s\n"
            % (
                md5,
                sp,
                seq,
                md5_in_xxx_samples,
                max(sample_counts),
                total_abundance,
                "\t".join(str(_) for _ in sample_counts),
            )
        )
        if worksheet:
            worksheet.write_string(current_row, 0, md5)
            worksheet.write_string(current_row, 1, sp)
            worksheet.write_string(current_row, 2, seq)
            worksheet.write_number(current_row, 3, md5_in_xxx_samples)
            worksheet.write_number(current_row, 4, max(sample_counts))
            worksheet.write_number(current_row, 5, total_abundance)
            for s, count in enumerate(sample_counts):
                worksheet.write_number(
                    current_row, LEADING_COLS + s, count, sample_formats[s]
                )
            current_row += 1
    del data

    if worksheet:
        worksheet.conditional_format(
            first_data_row,
            LEADING_COLS,
            current_row,
            LEADING_COLS + len(samples),
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
