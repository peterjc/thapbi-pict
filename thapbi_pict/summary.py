# Copyright 2019-2022 by Peter Cock, The James Hutton Institute.
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
from collections import defaultdict

import xlsxwriter
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .prepare import load_fasta_header
from .utils import color_bands
from .utils import file_to_sample_name
from .utils import find_requested_files
from .utils import load_metadata
from .utils import parse_species_tsv
from .utils import split_read_name_abundance

MISSING_META = ""
MISSING_DATA = "-"


def load_fasta_headers(sample_to_filenames, fields, default=""):
    """Load requested fields from FASTA headers.

    Argument sample_to_filenames is a dict of sample names as keys, with lists
    of filenames as values. Arguments fields is a list of field names. Default
    is a single value marker.

    WARNING: Currently looks at the first file for each sample only.

    Returns a dict with sample names as keys, and a list of the requested
    fields as values.

    If all the values are missing and/or match the default, raises a KeyError
    for that sample.
    """
    answer = {}
    blanks = [default] * len(fields)
    for sample, filenames in sample_to_filenames.items():
        if filenames:
            filename = filenames[0]  # TODO: Ignores multiple FASTA files!
            headers = load_fasta_header(filename)
            values = [headers.get(_, default) for _ in fields]
            if values == blanks:
                raise KeyError(filename)
            answer[sample] = values
    return answer


def _sp_display(species):
    """Format species classification for reports, see also _sp_sort_key."""
    if " " in species:
        return species
    elif species:
        return species + " (unknown species)"
    else:
        return "Unknown"


def _sp_sort_key(species):
    """Sort unknowns after knowns, see also _sp_display.

    Want this order:

      * Genus1 species1
      * Genus1 species2
      * Genus1 (unknown species) -- aka "Genus1"
      * Genus2 species1
      * ...
      * GenusN species1
      * ...
      * GenusN (unknown species) -- aka "GenusN"
      * Unknown - aka ""

    See also function _sp_display which handles the genus only and unknown.
    """
    assert species != "Unknown" and not species.endswith(" (unknown species)"), species
    if " " in species:
        return species
    elif species:
        return species + " {unknown species}"
    else:
        return "{Unknown}"


def sample_summary(
    sample_species_counts,
    meta_to_stem,
    stem_to_meta,
    meta_names,
    group_col,
    sample_stats,
    stats_fields,
    show_unsequenced,
    output,
    excel,
    method,
    min_abundance=1,
    debug=False,
):
    """Create samples (rows) vs species (cols) report.

    The expectation is that the inputs represent all the samples from
    a meaningful group, likely from multiple sequencing runs (plates).
    """
    if not output:
        sys.exit("ERROR: No output file specified.\n")

    species_predictions = set()  # includes A;B;C ambiguous entries
    for sample in sample_species_counts:
        for sp_list, _ in sample_species_counts[sample].items():
            species_predictions.add(sp_list)  # as string with any ; included
            genus_list = {sp_list.split(" ", 1)[0] for sp in sp_list.split(";")}
            if len(genus_list) > 1:
                sys.stderr.write(
                    f"WARNING: Conflicting genus from {sample}: {genus_list}\n"
                )

    species_predictions = sorted(species_predictions, key=_sp_sort_key)

    missing_stats = [MISSING_DATA] * len(stats_fields)

    # Open files and write headers
    # ============================

    # TSV header
    # ----------
    handle = open(output, "w")
    handle.write(
        "#%sSequencing sample\tClassification summary\t%sAccepted\tUnique\t%s\n"
        % (
            "\t".join(meta_names) + "\t" if meta_names else "",
            "\t".join(stats_fields) + "\t" if stats_fields else "",
            "\t".join(_sp_display(_) for _ in species_predictions),
        )
    )

    # Excel setup
    # -----------
    workbook = xlsxwriter.Workbook(excel)
    worksheet = workbook.add_worksheet("Samples vs species")
    red_conditional_format = workbook.add_format(
        # Maraschino red
        {"bg_color": "#FF2600", "font_color": "#000000"}
    )
    grey_conditional_format = workbook.add_format(
        {"bg_color": "#D3D3D3", "font_color": "#000000"}
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
    if meta_names:
        group_value = []
        for metadata, sample_batch in meta_to_stem.items():
            if sample_batch:
                group_value += [metadata[group_col]] * len(sample_batch)
            elif show_unsequenced:
                group_value += [metadata[group_col]]
        sample_formats = color_bands(
            group_value,
            sample_color_bands,
            debug=debug,
        )
        del group_value
    else:
        sample_formats = []

    # Excel header
    # ------------
    current_row = 0
    col_offset = len(meta_names) + 2 + len(stats_fields)
    # Set first row to be tall, with vertical text
    worksheet.set_row(0, 150, vertical_text_format)
    # If there are lots of species, set narrow column widths
    cols = len(species_predictions)
    if cols > 50:
        # Set column width to 2
        worksheet.set_column(col_offset + 2, col_offset + 1 + cols, 2)
    elif cols > 20:
        # Set column width to 4
        worksheet.set_column(col_offset + 2, col_offset + 1 + cols, 4)
    del cols
    for offset, name in enumerate(meta_names):
        worksheet.write_string(current_row, offset, name)
    col_offset = len(meta_names)
    worksheet.write_string(current_row, col_offset, "Sequencing sample")
    col_offset += 1
    worksheet.write_string(current_row, col_offset, "Classification summary")
    col_offset += 1
    for offset, name in enumerate(stats_fields):
        worksheet.write_string(current_row, col_offset + offset, name)
    col_offset += len(stats_fields)
    worksheet.write_string(current_row, col_offset, "Accepted")  # offset reference!
    worksheet.write_string(current_row, col_offset + 1, "Unique")  # offset reference!
    for offset, sp in enumerate(species_predictions):
        worksheet.write_string(current_row, col_offset + 2 + offset, _sp_display(sp))
    worksheet.freeze_panes(current_row + 1, col_offset + 2)

    # Main body
    # =========
    # Note already sorted on metadata values, discarded the order in the table
    for metadata, sample_batch in meta_to_stem.items():
        if not sample_batch and not show_unsequenced:
            # Nothing sequenced for this metadata entry,don't report it
            continue
        if not sample_batch and show_unsequenced:
            # Missing data in TSV:
            blanks = len(stats_fields) + 4 + len(species_predictions)
            # Using "-" for missing data, could use "NA" or "?"
            handle.write("\t".join(metadata) + ("\t" + MISSING_DATA) * blanks + "\n")
            # Missing data in Excel:
            try:
                cell_format = sample_formats[current_row]  # sample number
            except IndexError:
                cell_format = None
            current_row += 1
            assert len(meta_names) == len(metadata)
            assert col_offset == len(meta_names) + 2 + len(stats_fields)
            for offset, value in enumerate(metadata):
                worksheet.write_string(current_row, offset, value, cell_format)
            for offset in range(blanks):
                # Using formula "=NA()" or value "#N/A" work with our
                # original simpler conditional formatting color rule
                worksheet.write_string(
                    current_row, len(meta_names) + offset, MISSING_DATA, cell_format
                )

        # Now do the samples in this batch
        for sample in sample_batch:
            all_sp = set()
            unambig_sp = set()
            if sample in sample_species_counts:
                for sp_list, count in sample_species_counts[sample].items():
                    if count:
                        sp_list = sp_list.split(";")
                        all_sp.update(sp_list)
                        if len(sp_list) == 1:
                            unambig_sp.add(sp_list[0])
            human_sp_list = [
                _sp_display(sp) if sp in unambig_sp else sp + "(*)"
                for sp in sorted(all_sp, key=_sp_sort_key)
            ]
            if not all_sp:
                human_sp_list = ["-"]  # To match other fields, not using N/A

            # TSV
            # ---
            if metadata:
                handle.write("\t".join(metadata) + "\t")
            handle.write(sample + "\t" + ", ".join(human_sp_list) + "\t")
            if stats_fields:
                handle.write(
                    "\t".join(str(_) for _ in sample_stats.get(sample, missing_stats))
                    + "\t"
                )
            if sample in sample_species_counts:
                handle.write(
                    "%i\t%i\t%s\n"
                    % (
                        sum(sample_species_counts[sample].values()),
                        sum(1 for _ in sample_species_counts[sample].values() if _),
                        "\t".join(
                            str(sample_species_counts[sample][sp])
                            for sp in species_predictions
                        ),
                    )
                )
            else:
                # e.g. unsequenced sample, use "-" for missing data
                handle.write(
                    MISSING_DATA
                    + ("\t" + MISSING_DATA) * (1 + len(species_predictions))
                    + "\n"
                )

            # Excel
            # -----
            try:
                cell_format = sample_formats[current_row]  # sample number
            except IndexError:
                cell_format = None
            current_row += 1
            assert len(meta_names) == len(metadata)
            assert col_offset == len(meta_names) + 2 + len(stats_fields)
            for offset, value in enumerate(metadata):
                worksheet.write_string(current_row, offset, value, cell_format)
            worksheet.write_string(current_row, len(meta_names), sample, cell_format)
            worksheet.write_string(
                current_row, len(meta_names) + 1, ", ".join(human_sp_list), cell_format
            )
            for offset, _ in enumerate(stats_fields):
                # Currently all the statistics are integers, or strings
                v = sample_stats.get(sample, missing_stats)[offset]
                if isinstance(v, int):
                    worksheet.write_number(
                        current_row,
                        len(meta_names) + 2 + offset,
                        sample_stats[sample][offset],
                        cell_format,
                    )
                else:
                    worksheet.write_string(
                        current_row,
                        len(meta_names) + 2 + offset,
                        v,
                        cell_format,
                    )
            if sample in sample_species_counts:
                worksheet.write_number(
                    current_row,
                    col_offset,
                    sum(sample_species_counts[sample].values()),
                    cell_format,
                )
                worksheet.write_number(
                    current_row,
                    col_offset + 1,
                    sum(1 for _ in sample_species_counts[sample].values() if _),
                    cell_format,
                )
                for offset, sp in enumerate(species_predictions):
                    worksheet.write_number(
                        current_row,
                        col_offset + 2 + offset,
                        sample_species_counts[sample][sp],
                        cell_format,
                    )
            else:
                for offset in range(2 + len(species_predictions)):
                    worksheet.write_string(
                        current_row,
                        col_offset + offset,
                        MISSING_DATA,
                        cell_format,
                    )

    # Defined first, but takes priority over later conditional rules:
    worksheet.conditional_format(
        1,
        col_offset - 2 - len(stats_fields),  # go back to sample name
        current_row,
        col_offset + 2 + len(species_predictions),
        {
            "type": "cell",
            "criteria": "equal to",
            "value": f'"{MISSING_DATA}"',
            "format": grey_conditional_format,
        },
    )
    worksheet.conditional_format(
        1,
        col_offset + 2,
        current_row,
        col_offset + 2 + len(species_predictions),
        {
            "type": "cell",
            "criteria": "greater than",
            "value": 0,
            "format": red_conditional_format,
        },
    )
    workbook.close()
    handle.close()


def read_summary(
    markers,
    marker_md5_to_seq,
    marker_md5_species,
    marker_md5_abundance,
    abundance_by_samples,
    stem_to_meta,
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
    """Create reads (rows) vs species (cols) report.

    The expectation is that the inputs represent all the samples
    from one (96 well) plate, or some other meaningful batch.
    """
    if not output:
        sys.exit("ERROR: No output file specified.\n")

    # Excel setup
    # -----------
    LEADING_COLS = 7
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
    if len(stem_to_meta) > 50:
        # Set column width to 2
        worksheet.set_column(LEADING_COLS, LEADING_COLS + len(stem_to_meta), 2)
    elif len(stem_to_meta) > 20:
        # Set column width to 4
        worksheet.set_column(LEADING_COLS, LEADING_COLS + len(stem_to_meta), 4)

    # TSV setup
    # ---------
    handle = open(output, "w")

    # Metadata rows (one column per sample)
    # -------------
    current_row = 0
    first_data_row = 0
    sample_formats = [None] * len(stem_to_meta)
    if meta_names:
        # Insert extra header rows at start for sample meta-data
        # Make a single metadata call for each sample
        meta = [stem_to_meta[sample] for sample in stem_to_meta]
        for i, name in enumerate(meta_names):
            handle.write(
                "#%s%s\t%s\n"
                % ("\t" * (LEADING_COLS - 1), name, "\t".join(_[i] for _ in meta))
            )
        sample_formats = color_bands(
            [stem_to_meta[_][group_col] for _ in stem_to_meta],
            sample_color_bands,
            debug=debug,
        )
        for i, name in enumerate(meta_names):
            worksheet.write_string(
                current_row + i, LEADING_COLS - 1, name, cell_rightalign_format
            )
            for s, value in enumerate(meta):
                worksheet.write_string(
                    current_row + i,
                    LEADING_COLS + s,
                    value[i],
                    sample_formats[s],
                )
        current_row += len(meta_names)
    if stats_fields:
        # Insert extra header rows at start for sample meta-data
        # Make a single metadata call for each sample
        meta = [
            sample_stats.get(sample, [MISSING_DATA] * len(stats_fields))
            for sample in stem_to_meta
        ]
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
                if isinstance(value[i], int):
                    worksheet.write_number(
                        current_row + i,
                        LEADING_COLS + s,
                        value[i],
                        sample_formats[s],
                    )
                else:
                    worksheet.write_string(
                        current_row + i,
                        LEADING_COLS + s,
                        value[i],
                        sample_formats[s],
                    )
        current_row += len(stats_fields)

    # TSV main header
    # ---------------
    handle.write(
        "#Marker\tMD5\t%s-predictions\tMarker-sequence\tSample-count"
        "\tMax-sample-abundance\tTotal-abundance\t%s\n"
        % (method, "\t".join(stem_to_meta))
    )
    handle.write(
        "TOTAL or MAX\t-\t-\t-\t%i\t%i\t%i\t%s\n"
        % (
            max(
                sum(
                    1
                    for sample in stem_to_meta
                    if (marker, md5, sample) in abundance_by_samples
                )
                for (marker, md5) in marker_md5_to_seq
            )
            if marker_md5_to_seq
            else 0,
            max(
                (
                    abundance_by_samples.get((marker, md5, sample), 0)
                    for (marker, md5) in marker_md5_to_seq
                    for sample in stem_to_meta
                ),
                default=0,
            ),
            sum(marker_md5_abundance.values()),
            "\t".join(
                str(
                    sum(
                        abundance_by_samples.get((marker, md5, sample), 0)
                        for (marker, md5) in marker_md5_to_seq
                    )
                )
                for sample in stem_to_meta
            ),
        )
    )

    # Excel main header
    # -----------------
    worksheet.write_string(current_row, 0, "Marker")
    worksheet.write_string(current_row, 1, "MD5")
    worksheet.write_string(current_row, 2, method + "-predictions")
    worksheet.write_string(current_row, 3, "Marker-Sequence")
    worksheet.write_string(current_row, 4, "Sample-count")
    worksheet.write_string(current_row, 5, "Max-sample-abundance")
    worksheet.write_string(current_row, 6, "Total-abundance")
    for s, sample in enumerate(stem_to_meta):
        worksheet.write_string(current_row, LEADING_COLS + s, sample, sample_formats[s])
    current_row += 1
    first_data_row = current_row
    worksheet.write_string(current_row, 0, "TOTAL or MAX")
    worksheet.write_string(current_row, 1, "-")
    worksheet.write_string(current_row, 2, "-")
    worksheet.write_string(current_row, 3, "-")
    worksheet.write_number(
        current_row,
        4,
        max(
            sum(
                1
                for sample in stem_to_meta
                if (marker, md5, sample) in abundance_by_samples
            )
            for (marker, md5) in marker_md5_to_seq
        )
        if marker_md5_to_seq
        else 0,
    )
    worksheet.write_number(
        current_row,
        5,
        max(
            (
                abundance_by_samples.get((marker, md5, sample), 0)
                for (marker, md5) in marker_md5_to_seq
                for sample in stem_to_meta
            ),
            default=0,
        ),
    )
    worksheet.write_number(current_row, 6, sum(marker_md5_abundance.values()))
    for s, sample in enumerate(stem_to_meta):
        worksheet.write_number(
            current_row,
            LEADING_COLS + s,
            sum(
                abundance_by_samples.get((marker, md5, sample), 0)
                for (marker, md5) in marker_md5_to_seq
            ),
            sample_formats[s],
        )
    current_row += 1
    # keep total line in view plus headers:
    worksheet.freeze_panes(current_row, LEADING_COLS)

    # Main body
    # ---------
    # Build the first few columns as a list of lists, which we can sort
    data = [
        [
            marker,
            md5,
            ";".join(sorted(marker_md5_species[marker, md5])),
            marker_md5_to_seq[marker, md5],
            sum(1 for _ in stem_to_meta if (marker, md5, _) in abundance_by_samples),
            total_abundance,
        ]
        for (marker, md5), total_abundance in marker_md5_abundance.items()
    ]
    # Sort on marker, species prediction (with blank last, sorting as tilde);
    # number of samples (decreasing), total abundance (decreasing), md5
    data.sort(
        key=lambda row: (row[0], row[2] if row[2] else "~", -row[4], -row[5], row[1])
    )
    for marker, md5, sp, seq, md5_in_xxx_samples, total_abundance in data:
        sample_counts = [
            abundance_by_samples.get((marker, md5, _), 0) for _ in stem_to_meta
        ]
        assert sample_counts, f"{marker} {md5} sample counts: {sample_counts!r}"
        handle.write(
            "%s\t%s\t%s\t%s\t%i\t%i\t%i\t%s\n"
            % (
                marker,
                md5,
                sp,
                seq,
                md5_in_xxx_samples,
                max(sample_counts),
                total_abundance,
                "\t".join(str(_) for _ in sample_counts),
            )
        )
        worksheet.write_string(current_row, 0, marker)
        worksheet.write_string(current_row, 1, md5)
        worksheet.write_string(current_row, 2, sp)
        worksheet.write_string(current_row, 3, seq)
        worksheet.write_number(current_row, 4, md5_in_xxx_samples)
        worksheet.write_number(current_row, 5, max(sample_counts))
        worksheet.write_number(current_row, 6, total_abundance)
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
        LEADING_COLS + len(stem_to_meta),
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
    report_stem,
    method,
    min_abundance=1,
    metadata_file=None,
    metadata_encoding=None,
    metadata_cols=None,
    metadata_groups=None,
    metadata_fieldnames=None,
    metadata_index=None,
    require_metadata=False,
    show_unsequenced=True,
    ignore_prefixes=None,
    debug=False,
):
    """Implement the ``thapbi_pict summary`` command.

    The expectation is that the inputs represent all the samples from
    a meaningful group, likely from multiple sequencing runs (plates).
    """
    # TODO - refactor the old separate reporting code
    assert isinstance(inputs, list)

    if report_stem.endswith(os.path.sep) or os.path.isdir(report_stem):
        sys.exit("ERROR: Summary requires an output filename stem, not a directory")

    (stem_to_meta, meta_to_stem, meta_names, group_col,) = load_metadata(
        metadata_file,
        metadata_encoding,
        metadata_cols,
        metadata_groups,
        metadata_fieldnames,
        metadata_index,
        ignore_prefixes=ignore_prefixes,
        debug=debug,
    )

    meta_default = tuple([MISSING_META] * len(meta_names))
    fasta_files = defaultdict(list)
    tsv_files = defaultdict(list)
    markers = set()

    for fasta_file in find_requested_files(inputs, ".fasta", ignore_prefixes, debug):
        sample = file_to_sample_name(fasta_file)
        if sample not in stem_to_meta:
            if require_metadata:
                if debug:
                    sys.stderr.write(
                        f"DEBUG: Missing required metadata for {fasta_file}\n"
                    )
                continue
            stem_to_meta[sample] = meta_default
            if meta_default in meta_to_stem:
                meta_to_stem[meta_default].append(sample)
            else:
                meta_to_stem[meta_default] = [sample]
            if debug and metadata_file:
                sys.stderr.write(f"DEBUG: Missing metadata for {sample}\n")
        fasta_files[sample].append(fasta_file)
    for tsv_file in find_requested_files(
        inputs, f".{method}.tsv", ignore_prefixes, debug
    ):
        sample = file_to_sample_name(tsv_file)
        if sample.endswith(".all_reads") or sample in fasta_files:
            tsv_files[sample].append(tsv_file)
        elif debug:
            sys.stderr.write(f"DEBUG: Ignoring {tsv_file}\n")

    if debug:
        sys.stderr.write(
            f"DEBUG: Have metadata for {len(stem_to_meta)} samples, found"
            f" {sum(len(_) for _ in fasta_files.values())} FASTA files, and"
            f" {sum(len(_) for _ in tsv_files.values())} TSV for method {method}\n"
        )
    if not tsv_files:
        sys.exit("ERROR: No input FASTA and/or TSV files found")

    marker_md5_abundance = Counter()
    abundance_by_samples = {}
    marker_md5_species = {}  # sp values are sets (e.g. ambiguous matches)
    marker_md5_to_seq = {}
    sample_species_counts = {}  # 2nd key is sp list with semi-colons

    if debug:
        sys.stderr.write(f"Loading predictions for {method}\n")
    for sample, predicted_files in tsv_files.items():
        # Get MD5 to species mapping from the TSV only
        assert sample.endswith(".all_reads") or sample in stem_to_meta, sample
        for predicted_file in predicted_files:
            for marker, name, _taxid, sp_list in parse_species_tsv(
                predicted_file, min_abundance
            ):
                assert marker, predicted_file
                markers.add(marker)
                md5, abundance = split_read_name_abundance(name)
                if min_abundance > 1 and abundance < min_abundance:
                    continue
                if (marker, md5) in marker_md5_species:
                    marker_md5_species[marker, md5].update(sp_list.split(";"))
                else:
                    marker_md5_species[marker, md5] = set(sp_list.split(";"))

    if len(markers) == 1:
        # Not loading the post-abundance-threshold count,
        # Count should match the Seq-count column, but will not if running
        # report with higher abundance threshold - simpler to exclude.
        # For the threshold we have to update this if the report is stricter...
        stats_fields = (
            "Raw FASTQ",
            "Flash",
            "Cutadapt",
            "Threshold pool",
            "Threshold",
            "Max non-spike",
            "Max spike-in",
            "Singletons",
        )
    else:
        # TODO: How best to show the cutadapt and threshold values (per marker)?
        stats_fields = ("Raw FASTQ", "Flash")
    stats_headers = [_.lower().replace(" ", "_") for _ in stats_fields]
    blank_stats = [-1] * len(stats_fields)
    sample_stats = {}

    if not markers:
        # Corner case of zero sequences in all markers?
        assert not marker_md5_species

    if debug:
        sys.stderr.write("Loading FASTA sequences and abundances\n")
    for sample, fasta_list in fasta_files.items():
        assert sample in stem_to_meta, sample
        assert not sample.endswith(".all_fasta")
        sample_species_counts[sample] = Counter()
        for fasta_file in fasta_list:
            fasta_header = load_fasta_header(fasta_file)
            assert "marker" in fasta_header, f"{fasta_file} gave {fasta_header}"
            marker = fasta_header["marker"]

            sample_stats[sample] = values = [
                fasta_header.get(_, -1) for _ in stats_headers
            ]
            if values == blank_stats:
                sys.stderr.write(
                    f"WARNING: Missing header in FASTA file(s): {fasta_file}\n"
                )
                stats_fields = []
                sample_stats = {}
            del values

            with open(fasta_file) as handle:
                for title, seq in SimpleFastaParser(handle):
                    # Note: marker might not be in markers if FASTA has no seqs
                    assert marker in markers, fasta_file
                    md5, abundance = split_read_name_abundance(title.split(None, 1)[0])
                    if min_abundance > 1 and abundance < min_abundance:
                        continue
                    if (marker, md5) not in marker_md5_species:
                        if len(tsv_files) == 1:
                            sys.exit(
                                f"ERROR: {marker} {md5} from {fasta_file} "
                                f"not in {list(tsv_files.values())[0]}"
                            )
                        else:
                            sys.exit(
                                f"ERROR: {marker} {md5} from {fasta_file} "
                                f"not in {len(tsv_files)} TSV files"
                            )
                    abundance_by_samples[marker, md5, sample] = abundance
                    marker_md5_abundance[marker, md5] += abundance
                    marker_md5_to_seq[marker, md5] = seq
                    sample_species_counts[sample][
                        ";".join(sorted(marker_md5_species[marker, md5]))
                    ] += abundance

    if "Threshold" in stats_fields:
        # Apply any over-ride min_abundance for the reports
        i = stats_fields.index("Threshold")
        for _ in sample_stats:
            if sample_stats[_][i] < min_abundance:
                sample_stats[_][i] = min_abundance

    bad_fields = []
    for i, field in enumerate(stats_fields):
        if -1 in (_[i] for _ in sample_stats.values()):
            bad_fields.append(i)
            if debug:
                sys.stderr.write(f"WARNING: Dropping {field} column as missing stats\n")
    if bad_fields:
        for key, value in sample_stats.items():
            sample_stats[key] = [_ for i, _ in enumerate(value) if i not in bad_fields]
        stats_fields = [_ for i, _ in enumerate(stats_fields) if i not in bad_fields]
    del bad_fields

    if "Threshold pool" in stats_fields:
        # Try to remove any common folder prefix like raw_data/
        i = stats_fields.index("Threshold pool")
        paths = {_[i] for _ in sample_stats.values()}
        common = os.path.commonpath(paths)
        if len(paths) > 1 and common:
            if debug:
                sys.stderr.write(
                    f"DEBUG: Dropping threshold pool common prefix {common}\n"
                )
            for values in sample_stats.values():
                values[i] = values[i][len(common) + 1 :]

    if debug:
        sys.stderr.write(
            f"DEBUG: Have the following stats fields:{';'.join(stats_fields)}\n"
        )

    sample_summary(
        sample_species_counts,
        meta_to_stem,
        stem_to_meta,
        meta_names,
        group_col,
        sample_stats,
        stats_fields,
        show_unsequenced=show_unsequenced,
        output=f"{report_stem}.samples.{method}.tsv",
        excel=f"{report_stem}.samples.{method}.xlsx",
        method=method,
        min_abundance=min_abundance,
        debug=debug,
    )
    sys.stderr.write(f"Wrote {report_stem}.samples.{method}.*\n")

    read_summary(
        sorted(markers),
        marker_md5_to_seq,
        marker_md5_species,
        marker_md5_abundance,
        abundance_by_samples,
        stem_to_meta,
        meta_names,
        group_col,
        sample_stats,
        stats_fields,
        output=f"{report_stem}.reads.{method}.tsv",
        excel=f"{report_stem}.reads.{method}.xlsx",
        method=method,
        min_abundance=min_abundance,
        debug=debug,
    )
    sys.stderr.write(f"Wrote {report_stem}.reads.{method}.*\n")

    return 0
