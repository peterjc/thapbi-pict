# Copyright 2019-2020 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Summarise classification results at sample and read level.

This implements the ``thapbi_pict summary ...`` command.

The code uses the term metadata to refer to the user-provided
information about each sample (via a plain text TSV table),
and statistics for the internally tracked information about
each sample like the number of raw reads in the original FASTQ
files (via header lines in the intermediate FASTA files).
"""
import os
import sys
from collections import Counter

import xlsxwriter
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .prepare import load_fasta_header
from .utils import abundance_from_read_name
from .utils import color_bands
from .utils import file_to_sample_name
from .utils import find_paired_files
from .utils import load_metadata
from .utils import parse_species_tsv
from .utils import sample_sort
from .utils import split_read_name_abundance


def load_fasta_headers(sample_to_filename, fields, default=""):
    """Load requested fields from FASTA headers.

    Argument sample_to_filename is a dict of sample names
    as keys, filenames as paths. Arguments fields is a list
    of field names. Default is a single value marker.

    Returns a dict with sample names as keys, and lists
    of the requested fields as values.

    If all the values are missing and/or match the default,
    raises a KeyError for that sample.
    """
    answer = {}
    blanks = [default] * len(fields)
    for sample, filename in sample_to_filename.items():
        headers = load_fasta_header(filename)
        answer[sample] = [headers.get(_, default) for _ in fields]
        if answer[sample] == blanks:
            raise KeyError(sample)
    return answer


def sample_summary(
    tsv_files,
    metadata,
    meta_names,
    group_col,
    sample_stats,
    stats_fields,
    output,
    excel,
    human_output,
    method,
    min_abundance=1,
    debug=False,
):
    """Implement the ``thapbi_pict sample-summary`` command.

    The expectation is that the inputs represent all the samples from
    a meaningful group, likely from multiple sequencing runs (plates).
    """
    if not (output or human_output):
        sys.exit("ERROR: No output file specified.\n")

    missing_meta = set()
    genus_predictions = set()
    sample_genus_counts = {}
    species_predictions = set()  # includes A;B;C ambiguous entries
    sample_species_counts = {}
    for sample, predicted_file in tsv_files.items():
        if sample not in metadata:
            missing_meta.add(sample)
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
            f"DEBUG: {len(tsv_files)} samples with predictions for"
            f" {len(genus_predictions)} genera\n"
        )

    # Open files and write headers
    # ============================

    # TSV header
    # ----------
    handle = open(output, "w")
    handle.write(
        "#%sSequencing sample\t%sSeq-count\t%s\t%s\n"
        % (
            "\t".join(meta_names) + "\t" if meta_names else "",
            "\t".join(stats_fields) + "\t" if stats_fields else "",
            "\t".join(_ if _ else "Unknown" for _ in genus_predictions),
            "\t".join(species_columns),
        )
    )

    # Excel setup
    # -----------
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
    sample_formats = color_bands(
        [metadata[_][group_col] for _ in metadata], sample_color_bands, debug=debug,
    )

    # Excel header
    # ------------
    current_row = 0
    col_offset = len(meta_names) + 1 + len(stats_fields)
    # Set first row to be tall, with vertical text
    worksheet.set_row(0, 150, vertical_text_format)
    # If there are lots of species, set narrow column widths
    cols = len(genus_predictions) + len(species_columns)
    if cols > 50:
        # Set column width to 2
        worksheet.set_column(col_offset + 1, col_offset + cols, 2)
    elif cols > 20:
        # Set column width to 4
        worksheet.set_column(col_offset + 1, col_offset + cols, 4)
    del cols
    for offset, name in enumerate(meta_names):
        worksheet.write_string(current_row, offset, name)
    col_offset = len(meta_names)
    worksheet.write_string(current_row, col_offset, "Sequencing sample")
    col_offset += 1
    for offset, name in enumerate(stats_fields):
        worksheet.write_string(current_row, col_offset + offset, name)
    col_offset += len(stats_fields)
    worksheet.write_string(current_row, col_offset, "Seq-count")  # offset reference!
    for offset, genus in enumerate(genus_predictions):
        worksheet.write_string(
            current_row, col_offset + 1 + offset, genus if genus else "Unknown"
        )
    for offset, sp in enumerate(species_columns):
        worksheet.write_string(
            current_row, col_offset + 1 + len(genus_predictions) + offset, sp
        )
    worksheet.freeze_panes(current_row + 1, col_offset + 1)

    # Human header
    # -------------
    human = open(human_output, "w")
    human.write(
        "NOTE: Species listed with (uncertain/ambiguous) in brackets are where "
        "sequences matched multiple species equally well. For example, "
        "Phytophthora andina, P. infestans, and P. ipomoeae, share an identical "
        "marker.\n\n"
    )

    # Main body
    # =========
    # Note already sorted on metadata values, discarded the order in the table
    batches = []
    current_batch = []
    current_meta = None
    for sample, meta in metadata.items():
        if meta == current_meta:
            current_batch.append(sample)
            continue
        if current_batch:
            batches.append((current_meta, current_batch))
        current_batch = [sample]
        current_meta = meta
    if current_batch:
        batches.append((current_meta, current_batch))
    del current_batch, current_meta

    if missing_meta:
        batches.append(([""] * len(meta_names), sample_sort(missing_meta)))

    for metadata, sample_batch in batches:
        if meta_names:
            # Write the human readable metadata header
            human.write("-" * 60 + "\n\n")
            if metadata:
                for name, value in zip(meta_names, metadata):
                    if value:
                        human.write(f"{name}: {value}\n")
                human.write("\n")
            else:
                human.write("Missing metadata\n\n")
            if not sample_batch:
                human.write("Has not been sequenced.\n\n")
        # Now do the samples in this batch
        for sample in sample_batch:
            # TSV
            # ---
            if metadata:
                handle.write("\t".join(metadata) + "\t")
            handle.write(sample + "\t")
            if stats_fields:
                handle.write("\t".join(str(_) for _ in sample_stats[sample]) + "\t")
            handle.write(
                "%i\t%s\t%s\n"
                % (
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

            # Excel
            # -----
            try:
                cell_format = sample_formats[current_row]  # sample number
            except IndexError:
                cell_format = None
            current_row += 1
            assert len(meta_names) == len(metadata)
            assert col_offset == len(meta_names) + 1 + len(stats_fields)
            for offset, value in enumerate(metadata):
                worksheet.write_string(current_row, offset, value, cell_format)
            worksheet.write_string(current_row, len(meta_names), sample, cell_format)
            for offset, _ in enumerate(stats_fields):
                worksheet.write_number(
                    current_row,
                    len(meta_names) + 1 + offset,
                    sample_stats[sample][offset],
                    cell_format,
                )
            worksheet.write_number(
                current_row,
                col_offset,
                sum(sample_species_counts[sample].values()),
                cell_format,
            )
            for offset, genus in enumerate(genus_predictions):
                worksheet.write_number(
                    current_row,
                    col_offset + 1 + offset,
                    sample_genus_counts[sample][genus],
                    cell_format,
                )
            for offset, sp in enumerate(species_columns):
                worksheet.write_number(
                    current_row,
                    col_offset + 1 + len(genus_predictions) + offset,
                    sample_species_counts[sample][sp],
                    cell_format,
                )

            # Human report
            # ------------
            all_sp = set()
            unambig_sp = set()
            if sample in sample_species_counts:
                for sp_list, count in sample_species_counts[sample].items():
                    if count:
                        sp_list = sp_list.split(";")
                        all_sp.update(sp_list)
                        if len(sp_list) == 1:
                            unambig_sp.add(sp_list[0])
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
    workbook.close()
    handle.close()
    human.close()


def read_summary(
    fasta_files,
    tsv_files,
    metadata,
    meta_names,
    group_col,
    sample_stats,
    stats_fields,
    output,
    method,
    min_abundance=1,
    excel=None,
    debug=False,
):
    """Implement the ``thapbi_pict plate-summary`` command.

    The expectation is that the inputs represent all the samples
    from one (96 well) plate, or some other meaningful batch.
    """
    if not output:
        sys.exit("ERROR: No output file specified.\n")

    missing_meta = set()
    md5_abundance = Counter()
    abundance_by_samples = {}
    md5_species = {}
    md5_to_seq = {}

    if debug:
        sys.stderr.write("Loading FASTA sequences and abundances\n")
    for sample, fasta_file in fasta_files.items():
        if sample not in metadata:
            missing_meta.add(sample)
        with open(fasta_file) as handle:
            for title, seq in SimpleFastaParser(handle):
                md5, abundance = split_read_name_abundance(title.split(None, 1)[0])
                if min_abundance > 1 and abundance < min_abundance:
                    continue
                abundance_by_samples[md5, sample] = abundance
                md5_abundance[md5] += abundance
                md5_to_seq[md5] = seq
                md5_species[md5] = set()

    if debug:
        sys.stderr.write(f"Loading predictions for {method}\n")
    for sample, predicted_file in tsv_files.items():
        # TODO: Look at taxid here?
        for name, _, sp in parse_species_tsv(predicted_file, min_abundance):
            md5, abundance = split_read_name_abundance(name)
            if min_abundance > 1 and abundance < min_abundance:
                continue
            assert abundance_by_samples[md5, sample] == abundance, name
            if sp:
                md5_species[md5].update(sp.split(";"))

    if missing_meta:
        for sample in sample_sort(missing_meta):
            metadata[sample] = [""] * len(meta_names)

    # Excel setup
    # -----------
    LEADING_COLS = 6
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
    if len(metadata) > 50:
        # Set column width to 2
        worksheet.set_column(LEADING_COLS, LEADING_COLS + len(metadata), 2)
    elif len(metadata) > 20:
        # Set column width to 4
        worksheet.set_column(LEADING_COLS, LEADING_COLS + len(metadata), 4)

    # TSV setup
    # ---------
    handle = open(output, "w")

    # Metadata rows (one column per sample)
    # -------------
    current_row = 0
    first_data_row = 0
    sample_formats = [None] * len(metadata)
    if meta_names:
        # Insert extra header rows at start for sample meta-data
        # Make a single metadata call for each sample
        meta = [metadata[sample] for sample in metadata]
        for i, name in enumerate(meta_names):
            handle.write(
                "#%s%s\t%s\n"
                % ("\t" * (LEADING_COLS - 1), name, "\t".join(_[i] for _ in meta))
            )
        sample_formats = color_bands(
            [metadata[_][group_col] for _ in metadata], sample_color_bands, debug=debug,
        )
        for i, name in enumerate(meta_names):
            worksheet.write_string(
                current_row + i, LEADING_COLS - 1, name, cell_rightalign_format
            )
            for s, value in enumerate(meta):
                worksheet.write_string(
                    current_row + i, LEADING_COLS + s, value[i], sample_formats[s],
                )
        current_row += len(meta_names)
    if stats_fields:
        # Insert extra header rows at start for sample meta-data
        # Make a single metadata call for each sample
        meta = [sample_stats[sample] for sample in metadata]
        for i, name in enumerate(stats_fields):
            handle.write(
                "#%s%s\t%s\n"
                % ("\t" * (LEADING_COLS - 1), name, "\t".join(str(_[i]) for _ in meta))
            )
        for i, name in enumerate(stats_fields):
            worksheet.write_string(
                current_row + i, LEADING_COLS - 1, name, cell_rightalign_format
            )
            for s, value in enumerate(meta):
                worksheet.write_number(
                    current_row + i, LEADING_COLS + s, value[i], sample_formats[s],
                )
        current_row += len(stats_fields)

    # TSV main header
    # ---------------
    handle.write(
        "#ITS1-MD5\t%s-predictions\tSequence\tSample-count"
        "\tMax-sample-abundance\tTotal-abundance\t%s\n" % (method, "\t".join(metadata))
    )
    handle.write(
        "TOTAL\t-\t-\t%i\t%i\t%i\t%s\n"
        % (
            sum(
                1
                for md5 in md5_to_seq
                for sample in metadata
                if (md5, sample) in abundance_by_samples
            ),
            max(
                (
                    abundance_by_samples.get((md5, sample), 0)
                    for md5 in md5_to_seq
                    for sample in metadata
                ),
                default=0,
            ),
            sum(md5_abundance.values()),
            "\t".join(
                str(
                    sum(
                        abundance_by_samples.get((md5, sample), 0) for md5 in md5_to_seq
                    )
                )
                for sample in metadata
            ),
        )
    )

    # Excel main header
    # -----------------
    worksheet.write_string(current_row, 0, "ITS1-MD5")
    worksheet.write_string(current_row, 1, method + "-predictions")
    worksheet.write_string(current_row, 2, "Sequence")
    worksheet.write_string(current_row, 3, "Sample-count")
    worksheet.write_string(current_row, 4, "Max-sample-abundance")
    worksheet.write_string(current_row, 5, "Total-abundance")
    for s, sample in enumerate(metadata):
        worksheet.write_string(current_row, LEADING_COLS + s, sample, sample_formats[s])
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
            for sample in metadata
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
                for sample in metadata
            ),
            default=0,
        ),
    )
    worksheet.write_number(current_row, 5, sum(md5_abundance.values()))
    for s, sample in enumerate(metadata):
        worksheet.write_number(
            current_row,
            LEADING_COLS + s,
            sum(abundance_by_samples.get((md5, sample), 0) for md5 in md5_to_seq),
            sample_formats[s],
        )
    current_row += 1
    worksheet.freeze_panes(current_row, 5)  # keep total line in view plus headers

    # Main body
    # ---------
    # Build the first few columns as a list of lists, which we can sort
    data = [
        [
            md5,
            ";".join(sorted(md5_species[md5])),
            md5_to_seq[md5],
            sum(1 for _ in metadata if (md5, _) in abundance_by_samples),
            total_abundance,
        ]
        for md5, total_abundance in md5_abundance.items()
    ]
    # Sort on species prediction (with blank last, sorting as tilde);
    # number of samples (decreasing), total abundance (decreasing), md5
    data.sort(key=lambda row: (row[1] if row[1] else "~", -row[3], -row[4], row[0]))
    for md5, sp, seq, md5_in_xxx_samples, total_abundance in data:
        sample_counts = [abundance_by_samples.get((md5, _), 0) for _ in metadata]
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

    worksheet.conditional_format(
        first_data_row,
        LEADING_COLS,
        current_row,
        LEADING_COLS + len(metadata),
        {
            "type": "cell",
            "criteria": "greater than",
            "value": 0,
            "format": red_conditional_format,
        },
    )

    handle.close()
    workbook.close()


def main(
    inputs,
    out_dir,
    report,
    method,
    min_abundance=1,
    metadata_file=None,
    metadata_cols=None,
    metadata_groups=None,
    metadata_fieldnames=None,
    metadata_index=None,
    require_metadata=False,
    ignore_prefixes=None,
    debug=False,
):
    """Implement the ``thapbi_pict summary`` command.

    The expectation is that the inputs represent all the samples from
    a meaningful group, likely from multiple sequencing runs (plates).
    """
    # TODO - refactor the old separate reporting code
    assert isinstance(inputs, list)

    (metadata, meta_names, group_col,) = load_metadata(
        metadata_file,
        metadata_cols,
        metadata_groups,
        metadata_fieldnames,
        metadata_index,
        metadata_sort=True,
        ignore_prefixes=ignore_prefixes,
        debug=debug,
    )

    fasta_files = {}
    tsv_files = {}
    for fasta_file, tsv_file in find_paired_files(
        inputs, ".fasta", f".{method}.tsv", ignore_prefixes, debug, strict=True
    ):
        sample = file_to_sample_name(fasta_file)
        if require_metadata and sample not in metadata:
            continue
        fasta_files[sample] = fasta_file
        tsv_files[sample] = tsv_file

    if debug:
        sys.stderr.write(
            f"Have metadata for {len(metadata)} samples,"
            f" found {len(fasta_files)} FASTA files,"
            f" and {len(tsv_files)} TSV for method {method}\n"
        )
    if not tsv_files:
        sys.exit("ERROR: No input FASTA and TSV files found")

    if report:
        stem = os.path.join(out_dir, report)
    else:
        # Include version number here?
        stem = os.path.join(out_dir, "thapbi-pict")

    # Not loading the post-abundance-threshold count, or the threshold.
    # Count should match the Seq-count column, but will not if running
    # report with higher abundance threshold - simpler to exclude them:
    stats_fields = ("Raw FASTQ", "Trimmomatic", "Flash", "Cutadapt")
    try:
        sample_stats = load_fasta_headers(
            fasta_files, ("raw_fastq", "trimmomatic", "flash", "cutadapt"), -1
        )
    except KeyError as err:
        sys.stderr.write(
            f"WARNING: Missing header information in FASTA file(s): {err}\n"
        )
        sample_stats = {}
        stats_fields = []

    bad_fields = []
    for i in range(len(stats_fields)):
        if -1 in (_[i] for _ in sample_stats.values()):
            bad_fields.append(i)
    if bad_fields:
        for sample, value in sample_stats.items():
            sample_stats[sample] = [
                _ for i, _ in enumerate(value) if i not in bad_fields
            ]
        stats_fields = [_ for i, _ in enumerate(stats_fields) if i not in bad_fields]
        sys.stderr.write(
            f"WARNING: Dropping {len(bad_fields)} missing stats column(s)\n"
        )
    del bad_fields

    sample_summary(
        tsv_files,
        metadata,
        meta_names,
        group_col,
        sample_stats,
        stats_fields,
        output=f"{stem}.samples.{method}.tsv",
        excel=f"{stem}.samples.{method}.xlsx",
        human_output=f"{stem}.samples.{method}.txt",
        method=method,
        min_abundance=min_abundance,
        debug=debug,
    )
    sys.stderr.write(f"Wrote {stem}.samples.{method}.*\n")

    read_summary(
        fasta_files,
        tsv_files,
        metadata,
        meta_names,
        group_col,
        sample_stats,
        stats_fields,
        output=f"{stem}.reads.{method}.tsv",
        excel=f"{stem}.reads.{method}.xlsx",
        method=method,
        min_abundance=min_abundance,
        debug=debug,
    )
    sys.stderr.write(f"Wrote {stem}.reads.{method}.*\n")

    return 0
