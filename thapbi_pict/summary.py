# Copyright 2019-2024 by Peter Cock, The James Hutton Institute.
# Revisions copyright 2024 by Peter Cock, University of Strathclyde.
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

from __future__ import annotations

import os
import sys
from collections import Counter

import xlsxwriter

from .utils import export_sample_biom
from .utils import find_requested_files
from .utils import load_metadata
from .utils import md5seq
from .utils import parse_sample_tsv

MISSING_META = ""
MISSING_DATA = "-"


def _sp_display(species: str) -> str:
    """Format species classification for reports, see also _sp_sort_key."""
    if " " in species:
        return species
    elif species:
        return species + " (unknown species)"
    else:
        return "Unknown"


def _sp_sort_key(species: str) -> str:
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


def color_bands(
    meta_groups, sample_color_bands, default_fmt=None, debug: bool = False
) -> list:
    """Return a list for formats, one for each sample."""
    default = [default_fmt] * len(meta_groups)

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
                "DEBUG: All samples had same metadata in color grouping field:"
                f" {meta_groups[0]!r}\n"
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
                f"DEBUG: {len(set(meta_groups))} groups not suitable for coloring:\n"
            )
            sys.stderr.write(";".join(sorted(set(meta_groups))) + "\n")
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
                    "WARNING: Metadata not grouped nicely for coloring:"
                    f" {', '.join(debug_msg)} and then {value} (again).\n"
                )
                return default
            bands.append(max(bands) + 1)
            debug_msg.append(value)
    assert len(set(bands)) == max(bands) + 1, bands
    return [sample_color_bands[_ % len(sample_color_bands)] for _ in bands]


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
        for sp_list in sample_species_counts[sample].keys():
            species_predictions.add(sp_list)  # as string with any ; included
            genus_list = {sp_list.split(" ", 1)[0] for sp in sp_list.split(";")}
            if len(genus_list) > 1:
                sys.stderr.write(
                    f"WARNING: Conflicting genus from {sample}: {genus_list}\n"
                )

    species_predictions = sorted(species_predictions, key=_sp_sort_key)

    # Open files and write headers
    # ============================

    # TSV header
    # ----------
    handle = open(output, "w")
    handle.write(
        "#%sSequencing sample\tClassification summary\t%s%s\n"
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
    thousands_format = workbook.add_format({"num_format": "#,###,###,##0"})
    sample_color_bands = [
        # Simple rolling rainbow pastel pallet
        workbook.add_format(
            {
                "bg_color": c,
                "font_color": "#000000",
                # Integers with thousands separators (regional settings respected):
                "num_format": "#,###,###,##0",
            }
        )
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
            thousands_format,
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
        worksheet.set_column(col_offset, col_offset + cols - 1, 2)
    elif cols > 20:
        # Set column width to 4
        worksheet.set_column(col_offset, col_offset + cols - 1, 4)
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
    for offset, sp in enumerate(species_predictions):
        worksheet.write_string(current_row, col_offset + offset, _sp_display(sp))
    assert col_offset == len(meta_names) + 2 + len(stats_fields)

    # Main body
    # =========
    # Note already sorted on metadata values, discarded the order in the metadata table

    # Columns are:
    # - len(meta_names) == len(metadata)
    # - Sequencing sample
    # - Classification summary
    # - len(stats_fields)
    # - len(species_predictions)
    assert col_offset == len(meta_names) + 2 + len(stats_fields)

    for metadata, sample_batch in meta_to_stem.items():
        if not sample_batch and not show_unsequenced:
            # Nothing sequenced for this metadata entry,don't report it
            continue
        if not sample_batch and show_unsequenced:
            # Missing data in TSV:
            blanks = 2 + len(stats_fields) + len(species_predictions)
            # Using "-" for missing data, could use "NA" or "?"
            handle.write("\t".join(metadata) + ("\t" + MISSING_DATA) * blanks + "\n")
            # Missing data in Excel:
            try:
                cell_format = sample_formats[current_row]  # sample number
            except IndexError:
                cell_format = thousands_format
            current_row += 1
            assert len(meta_names) == len(metadata)
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
            if stats_fields and sample_stats:
                handle.write(
                    "\t".join(
                        str(sample_stats.get(sample, {}).get(field, MISSING_DATA))
                        for field in stats_fields
                    )
                    + "\t"
                )
            if sample in sample_species_counts:
                handle.write(
                    "\t".join(
                        str(sample_species_counts[sample][sp])
                        for sp in species_predictions
                    )
                    + "\n"
                )
            else:
                # e.g. unsequenced sample, use "-" for missing data
                handle.write(
                    "\t".join([MISSING_DATA] * len(species_predictions)) + "\n"
                )

            # Excel
            # -----
            try:
                cell_format = sample_formats[current_row]  # sample number
            except IndexError:
                cell_format = thousands_format
            current_row += 1
            assert len(meta_names) == len(metadata)
            assert col_offset == len(meta_names) + 2 + len(stats_fields)
            for offset, value in enumerate(metadata):
                worksheet.write_string(current_row, offset, value, cell_format)
            worksheet.write_string(current_row, len(meta_names), sample, cell_format)
            worksheet.write_string(
                current_row, len(meta_names) + 1, ", ".join(human_sp_list), cell_format
            )
            for offset, field in enumerate(stats_fields):
                # Currently all the statistics are integers, or strings
                v = sample_stats.get(sample, {}).get(field, MISSING_DATA)
                if isinstance(v, int):
                    worksheet.write_number(
                        current_row,
                        len(meta_names) + 2 + offset,
                        v,
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
                for offset, sp in enumerate(species_predictions):
                    worksheet.write_number(
                        current_row,
                        col_offset + offset,
                        sample_species_counts[sample][sp],
                        cell_format,
                    )
            else:
                for offset in range(len(species_predictions)):
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
        col_offset + len(species_predictions) - 1,
        {
            "type": "cell",
            "criteria": "equal to",
            "value": f'"{MISSING_DATA}"',
            "format": grey_conditional_format,
        },
    )
    worksheet.conditional_format(
        1,
        col_offset,
        current_row,
        col_offset + len(species_predictions) - 1,
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
    biom=None,  # filename
    debug=False,
) -> None:
    """Create reads (rows) vs species (cols) report.

    The expectation is that the inputs represent all the samples
    from one (96 well) plate, or some other meaningful batch.
    """
    if biom and marker_md5_to_seq and stem_to_meta:
        if export_sample_biom(
            biom,
            marker_md5_to_seq,
            # The only sequence metadata is our classification,
            # Use field name f"{method}-predictions" to match TSV/Excel?
            # TODO - include list of taxids?
            {
                key: {"genus-species": value}
                for key, value in marker_md5_species.items()
            },
            # User-supplied sample data, plus stats from the pipeline:
            {
                sample: dict(
                    list(zip(meta_names, stem_to_meta[sample]))
                    + list(zip(stats_fields, sample_stats[sample]))
                )
                for sample in sample_stats
            },
            abundance_by_samples,
        ):
            if debug:
                sys.stderr.write(f"DEBUG: Wrote {biom}\n")
        else:
            sys.exit("ERROR: Missing optional Python library for BIOM output")

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
    thousands_format = workbook.add_format({"num_format": "#,###,###,##0"})
    sample_color_bands = [
        # Simple rolling rainbow pastel pallet
        workbook.add_format(
            {
                "bg_color": c,
                "font_color": "#000000",
                # Integers with thousands separators (regional settings respected):
                "num_format": "#,###,###,##0",
            }
        )
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
    sample_formats = [thousands_format] * len(stem_to_meta)
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
            thousands_format,
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
            [
                sample_stats.get(sample, {}).get(field, MISSING_DATA)
                for field in stats_fields
            ]
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
        "MAX or TOTAL\t-\t-\t-\t%i\t%i\t%i\t%s\n"
        % (
            # max Sample-count
            (
                max(
                    sum(
                        1
                        for sample in stem_to_meta
                        if abundance_by_samples.get((marker, md5, sample), 0)
                    )
                    for (marker, md5) in marker_md5_to_seq
                )
                if marker_md5_to_seq
                else 0
            ),
            # max of Max-sample-abundance
            max(
                (
                    abundance_by_samples.get((marker, md5, sample), 0)
                    for (marker, md5) in marker_md5_to_seq
                    for sample in stem_to_meta
                ),
                default=0,
            ),
            # sum of Total-abundance
            sum(marker_md5_abundance.values()),
            "\t".join(
                str(
                    max(
                        (
                            abundance_by_samples.get((marker, md5, sample), 0)
                            for (marker, md5) in marker_md5_to_seq
                        ),
                        default=0,
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
    worksheet.write_string(current_row, 0, "MAX or TOTAL")
    worksheet.write_string(current_row, 1, "-")
    worksheet.write_string(current_row, 2, "-")
    worksheet.write_string(current_row, 3, "-")
    worksheet.write_number(
        current_row,
        4,
        # max Sample-count:
        (
            max(
                sum(
                    1
                    for sample in stem_to_meta
                    if abundance_by_samples.get((marker, md5, sample), 0)
                )
                for (marker, md5) in marker_md5_to_seq
            )
            if marker_md5_to_seq
            else 0
        ),
        thousands_format,
    )
    worksheet.write_number(
        current_row,
        5,
        # max Max-sample-abundance
        max(
            (
                abundance_by_samples.get((marker, md5, sample), 0)
                for (marker, md5) in marker_md5_to_seq
                for sample in stem_to_meta
            ),
            default=0,
        ),
        thousands_format,
    )
    worksheet.write_number(
        # sum of Total-abundance
        current_row,
        6,
        sum(marker_md5_abundance.values()),
        thousands_format,
    )
    for s, sample in enumerate(stem_to_meta):
        worksheet.write_number(
            current_row,
            LEADING_COLS + s,
            max(
                (
                    abundance_by_samples.get((marker, md5, sample), 0)
                    for (marker, md5) in marker_md5_to_seq
                ),
                default=0,
            ),
            sample_formats[s],
        )
    current_row += 1

    # Main body
    # ---------
    # Build the first few columns as a list of lists, which we can sort
    data = [
        [
            marker,
            md5,
            marker_md5_species[marker, md5],
            marker_md5_to_seq[marker, md5],
            sum(
                1
                for sample in stem_to_meta
                if abundance_by_samples.get((marker, md5, sample), 0)
            ),
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
            abundance_by_samples.get((marker, md5, sample), 0)
            for sample in stem_to_meta
        ]
        if not sample_counts:
            # Must have lost it due to abundance filter?
            assert min_abundance, f"{marker} {md5} gave zero sample counts!"
            continue
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
        worksheet.write_number(current_row, 4, md5_in_xxx_samples, thousands_format)
        worksheet.write_number(current_row, 5, max(sample_counts), thousands_format)
        worksheet.write_number(current_row, 6, total_abundance, thousands_format)
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
    report_stem: str,
    method: str,
    min_abundance: int = 1,
    metadata_file: str | None = None,
    metadata_encoding: str | None = None,
    metadata_cols: str | None = None,
    metadata_groups: str | None = None,
    metadata_fieldnames: str | None = None,
    metadata_index: str | None = None,
    require_metadata: bool = False,
    show_unsequenced: bool = True,
    ignore_prefixes: tuple[str] | None = None,
    biom: bool = False,
    debug: bool = False,
) -> int:
    """Implement the ``thapbi_pict summary`` command.

    The expectation is that the inputs represent all the samples from
    a meaningful group, likely from multiple sequencing runs (plates).
    """
    # TODO - refactor the old separate reporting code
    assert isinstance(inputs, list)

    if report_stem.endswith(os.path.sep) or os.path.isdir(report_stem):
        sys.exit("ERROR: Summary requires an output filename stem, not a directory")

    (
        stem_to_meta,
        meta_to_stem,
        meta_names,
        group_col,
    ) = load_metadata(
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
    markers: set[str] = set()

    classifications_tsv = [
        _
        for _ in find_requested_files(inputs, f".{method}.tsv", ignore_prefixes, debug)
        # These are what we'll call the output files:
        if not _.endswith((f".reads.{method}.tsv", f".samples.{method}.tsv"))
    ]
    if not classifications_tsv:
        sys.exit(f"ERROR: No classification files matching *.{method}.tsv")

    if debug:
        sys.stderr.write(
            f"DEBUG: Have metadata for {len(stem_to_meta)} samples, found "
            f"{len(classifications_tsv)} classifier files\n"
        )

    marker_md5_abundance: dict[tuple[str, str], int] = Counter()
    abundance_by_samples: dict[tuple[str, str, str], int] = {}
    marker_md5_species: dict[
        tuple[str, str], str
    ] = {}  # sp values are ";" separated sorted (e.g. ambiguous matches)
    marker_md5_to_seq: dict[tuple[str, str], str] = {}
    sample_species_counts: dict[
        str, dict[str, int]
    ] = {}  # 2nd key is sp list with semi-colons

    # Not loading the post-abundance-threshold count,
    # Count should match the Seq-count column, but will not if running
    # report with higher abundance threshold - simpler to exclude.
    # For the threshold we have to update this if the report is stricter...
    stats_fields: tuple[str, ...] = (
        "Raw FASTQ",
        "Flash",
        "Cutadapt",
        "Threshold pool",
        "Threshold",
        "Control",
        "Max non-spike",
        "Max spike-in",
        "Singletons",
    )
    blank_stat = -1
    # nested dict, keys are sample, then field name:
    sample_stats: dict[str, dict[str, str | int]] = {}

    if not markers:
        # Corner case of zero sequences in all markers?
        assert not marker_md5_species

    for filename in classifications_tsv:
        if debug:
            sys.stderr.write(f"DEBUG: Loading samples sequences etc from {filename}\n")
        seqs, seq_meta, sample_headers, counts = parse_sample_tsv(
            filename, min_abundance=min_abundance, debug=debug
        )
        if require_metadata:
            # Drop unwanted samples & thus unwanted md5
            counts = {
                (marker, md5, sample): a
                for (marker, md5, sample), a in counts.items()
                if sample in stem_to_meta
            }
            sample_headers = {
                sample: fasta_header
                for (sample, fasta_header) in sample_headers.items()
                if sample in stem_to_meta
            }
            # Can now potentially drop some sequences:
            totals: dict[tuple[str, str], int] = Counter()
            for (idn, marker, sample), a in counts.items():
                assert sample in stem_to_meta  # we filtered counts earlier
                if a:
                    totals[idn, marker] += a
            seqs = {
                (idn, marker): seq
                for (idn, marker), seq in seqs.items()
                if totals[idn, marker]
            }
            seq_meta = {
                (idn, marker): value
                for (idn, marker), value in seq_meta.items()
                if totals[idn, marker]
            }
            del totals
        for sample, fasta_header in sample_headers.items():
            if sample not in stem_to_meta:
                if require_metadata:
                    sys.exit(f"ERROR: Failed to drop {sample} with no metadata\n")
                stem_to_meta[sample] = meta_default
                if meta_default in meta_to_stem:
                    meta_to_stem[meta_default].append(sample)
                else:
                    meta_to_stem[meta_default] = [sample]
                if debug and metadata_file:
                    sys.stderr.write(f"DEBUG: Missing metadata for {sample}\n")
            if sample not in sample_species_counts:
                sample_species_counts[sample] = Counter()
            if sample not in sample_stats:
                sample_stats[sample] = {}
            for field in stats_fields:
                value = fasta_header.get(field, blank_stat)
                if field in ("Sample", "Control", "Threshold pool"):
                    # Text fields; should be consistent between markers
                    if field in sample_stats[sample]:
                        assert sample_stats[sample][field] == value
                    else:
                        sample_stats[sample][field] = value
                else:
                    # Integer fields
                    value = int(value)
                    if value == blank_stat:
                        pass
                    elif field in ("Cutadapt", "Singletons"):
                        # Sum over the markers
                        sample_stats[sample][field] = (
                            int(sample_stats[sample].get(field, 0)) + value
                        )
                    elif field in ("Threshold"):
                        # Report value if shared between markers
                        if field not in sample_stats[sample]:
                            sample_stats[sample][field] = value
                        elif sample_stats[sample][field] != value:
                            sample_stats[sample][field] = "By marker"
                    elif field.startswith("Max "):
                        # Take largest value over the markers
                        if field not in sample_stats[sample]:
                            sample_stats[sample][field] = value
                        else:
                            sample_stats[sample][field] = max(
                                value, sample_stats[sample][field]
                            )
                    else:
                        # e.g. Threshold pool
                        # Should be consistent between markers
                        if field in sample_stats[sample]:
                            assert sample_stats[sample][field] == value
                        else:
                            sample_stats[sample][field] = value
                del value
            if (
                not sample_stats[sample]
                or set(sample_stats[sample].values()) == blank_stat
            ):
                sys.stderr.write(f"WARNING: Missing stats header for {sample}\n")

        if debug:
            sys.stderr.write(f"DEBUG: Loaded sample headers from {filename}\n")

        for (marker, md5), seq in seqs.items():
            markers.add(marker)
            assert md5 == md5seq(seq), (marker, md5, filename)
            marker_md5_to_seq[marker, md5] = seq
            marker_md5_species[marker, md5] = ";".join(
                sorted(set(seq_meta[marker, md5]["genus-species"].split(";")))
            )
            for sample in sample_headers:
                if require_metadata and sample not in stem_to_meta:
                    sys.exit(f"ERROR: Should have dropped {sample} with no metadata")
                try:
                    abundance = counts.pop((marker, md5, sample))  # empty the dict
                except KeyError:
                    abundance = 0
                assert (marker, md5, sample) not in abundance_by_samples
                abundance_by_samples[marker, md5, sample] = abundance
                marker_md5_abundance[marker, md5] += abundance
                # assert sample in sample_species_counts, sorted(
                #    sample_species_counts.keys()
                # )
                sample_species_counts[sample][marker_md5_species[marker, md5]] += (
                    abundance
                )
        del seqs, sample_headers, counts

    if debug:
        sys.stderr.write(
            f"DEBUG: Loaded sample counts for {len(sample_species_counts)} sequences\n"
        )

    if "Threshold" in stats_fields:
        # Apply any over-ride min_abundance for the reports
        for sample in sample_stats:
            if (
                sample_stats[sample]["Threshold"] != "By marker"
                and int(sample_stats[sample]["Threshold"]) < min_abundance
            ):
                sample_stats[sample]["Threshold"] = min_abundance

    bad_fields = []
    for field in stats_fields:
        if blank_stat in {_.get(field, blank_stat) for _ in sample_stats.values()}:
            bad_fields.append(field)
            sys.stderr.write(f"WARNING: Dropping {field} column as missing stats\n")
    if bad_fields:
        for sample in sample_stats:
            for field in bad_fields:
                if field in sample_stats[sample]:
                    del sample_stats[sample][field]
        stats_fields = tuple(field for field in stats_fields if field not in bad_fields)
    del bad_fields

    if "Threshold pool" in stats_fields:
        # Try to remove any common folder prefix like raw_data/
        i: int = stats_fields.index("Threshold pool")
        paths = {_["Threshold pool"] for _ in sample_stats.values()}
        common = os.path.commonpath([str(_) for _ in paths])  # make type explicit!
        if len(paths) > 1 and common:
            if debug:
                sys.stderr.write(
                    f"DEBUG: Dropping threshold pool common prefix {common}\n"
                )
            for values in sample_stats.values():
                assert isinstance(values, list)
                values[i] = values[i][len(common) + 1 :]

    # Both of these are potentially altered by the min_abundance argument,
    # simplest to calculate them here rather than in the separate reports:
    stats_fields = tuple(list(stats_fields) + ["Accepted", "Unique"])
    tmp_accepted_by_sample: dict[str, list[int]] = {}
    for (marker, _md5, sample), a in abundance_by_samples.items():
        if a >= min_abundance and marker in markers:
            try:
                tmp_accepted_by_sample[sample].append(a)
            except KeyError:
                tmp_accepted_by_sample[sample] = [a]
    for sample in sample_stats:
        asv_counts = tmp_accepted_by_sample.get(sample, [])
        sample_stats[sample]["Accepted"] = sum(asv_counts)
        sample_stats[sample]["Unique"] = len(asv_counts)
        del asv_counts
    del tmp_accepted_by_sample

    if require_metadata:
        assert set(sample_stats) == set(sample_species_counts) == set(stem_to_meta)
    else:
        assert set(sample_stats) == set(sample_species_counts), (
            f"{sorted(set(sample_stats))} vs {sorted(set(sample_species_counts))}"
        )
        assert len(sample_stats) <= len(stem_to_meta)

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

    del sample_species_counts, meta_to_stem

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
        biom=f"{report_stem}.reads.{method}.biom" if biom else None,
        method=method,
        min_abundance=min_abundance,
        debug=debug,
    )
    sys.stderr.write(f"Wrote {report_stem}.reads.{method}.*\n")

    return 0
