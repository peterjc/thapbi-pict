# Copyright 2019-2021 by Peter Cock, The James Hutton Institute.
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

    # Open files and write headers
    # ============================

    # TSV header
    # ----------
    handle = open(output, "w")
    handle.write(
        "#%sSequencing sample\tClassification summary\t%sRead count\t%s\n"
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
    worksheet.write_string(current_row, col_offset, "Classification summary")
    col_offset += 1
    for offset, name in enumerate(stats_fields):
        worksheet.write_string(current_row, col_offset + offset, name)
    col_offset += len(stats_fields)
    worksheet.write_string(current_row, col_offset, "Read count")  # offset reference!
    for offset, sp in enumerate(species_predictions):
        worksheet.write_string(current_row, col_offset + 1 + offset, _sp_display(sp))
    worksheet.freeze_panes(current_row + 1, col_offset + 1)

    # Human header
    # -------------
    human = open(human_output, "w")
    human.write(
        "NOTE: Species listed with (*) are where sequences matched multiple "
        "species equally well. For example, Phytophthora andina, P. infestans, "
        "and P. ipomoeae, share an identical ITS1 marker.\n\n"
    )

    # Main body
    # =========
    # Note already sorted on metadata values, discarded the order in the table
    for metadata, sample_batch in meta_to_stem.items():
        if not sample_batch and not show_unsequenced:
            # Nothing sequenced for this metadata entry,don't report it
            continue
        # Write the human readable metadata header
        if meta_names:
            human.write("-" * 60 + "\n\n")
            if any(metadata):
                for name, value in zip(meta_names, metadata):
                    if value:
                        human.write(f"{name}: {value}\n")
                human.write("\n")
            else:
                # Could report, but redundant with "Sequencing sample: ..."
                # human.write("Missing metadata\n\n")
                pass
        if not sample_batch and show_unsequenced:
            human.write("Has not been sequenced.\n\n")
            # Missing data in TSV:
            blanks = len(stats_fields) + 3 + len(species_predictions)
            # Using "-" for missing data, could use "NA" or "?"
            handle.write("\t".join(metadata) + ("\t-") * blanks + "\n")
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
                    current_row, len(meta_names) + offset, "-", cell_format
                )

        # Now do the samples in this batch
        for sample in sample_batch:
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
            human_sp_list = [
                _sp_display(sp) if sp in unambig_sp else sp + "(*)"
                for sp in sorted(all_sp, key=_sp_sort_key)
            ]
            for sp in human_sp_list:
                human.write(f" - {sp}\n")
            if not all_sp:
                human.write(" - No data\n")
                human_sp_list = ["-"]  # To match other fields, not using N/A
            human.write("\n")
            del all_sp, unambig_sp

            # TSV
            # ---
            if metadata:
                handle.write("\t".join(metadata) + "\t")
            handle.write(sample + "\t" + ", ".join(human_sp_list) + "\t")
            if stats_fields:
                handle.write("\t".join(str(_) for _ in sample_stats[sample]) + "\t")
            handle.write(
                "%i\t%s\n"
                % (
                    sum(sample_species_counts[sample].values()),
                    "\t".join(
                        str(sample_species_counts[sample][sp])
                        for sp in species_predictions
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
            assert col_offset == len(meta_names) + 2 + len(stats_fields)
            for offset, value in enumerate(metadata):
                worksheet.write_string(current_row, offset, value, cell_format)
            worksheet.write_string(current_row, len(meta_names), sample, cell_format)
            worksheet.write_string(
                current_row, len(meta_names) + 1, ", ".join(human_sp_list), cell_format
            )
            for offset, _ in enumerate(stats_fields):
                worksheet.write_number(
                    current_row,
                    len(meta_names) + 2 + offset,
                    sample_stats[sample][offset],
                    cell_format,
                )
            worksheet.write_number(
                current_row,
                col_offset,
                sum(sample_species_counts[sample].values()),
                cell_format,
            )
            for offset, sp in enumerate(species_predictions):
                worksheet.write_number(
                    current_row,
                    col_offset + 1 + offset,
                    sample_species_counts[sample][sp],
                    cell_format,
                )

    # Defined first, but takes priority over later conditional rules:
    worksheet.conditional_format(
        1,
        col_offset - 2 - len(stats_fields),  # go back to sample name
        current_row,
        col_offset + 1 + len(species_predictions),
        {
            "type": "cell",
            "criteria": "equal to",
            "value": '"-"',
            "format": grey_conditional_format,
        },
    )
    worksheet.conditional_format(
        1,
        col_offset + 1,
        current_row,
        col_offset + 1 + len(species_predictions),
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
    md5_to_seq,
    md5_species,
    md5_abundance,
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
    """Implement the ``thapbi_pict plate-summary`` command.

    The expectation is that the inputs represent all the samples
    from one (96 well) plate, or some other meaningful batch.
    """
    if not output:
        sys.exit("ERROR: No output file specified.\n")

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
        meta = [sample_stats[sample] for sample in stem_to_meta]
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
                    current_row + i,
                    LEADING_COLS + s,
                    value[i],
                    sample_formats[s],
                )
        current_row += len(stats_fields)

    # TSV main header
    # ---------------
    handle.write(
        "#ITS1-MD5\t%s-predictions\tSequence\tSample-count"
        "\tMax-sample-abundance\tTotal-abundance\t%s\n"
        % (method, "\t".join(stem_to_meta))
    )
    handle.write(
        "TOTAL or MAX\t-\t-\t%i\t%i\t%i\t%s\n"
        % (
            max(
                sum(
                    1
                    for sample in stem_to_meta
                    if (md5, sample) in abundance_by_samples
                )
                for md5 in md5_to_seq
            ),
            max(
                (
                    abundance_by_samples.get((md5, sample), 0)
                    for md5 in md5_to_seq
                    for sample in stem_to_meta
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
                for sample in stem_to_meta
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
    for s, sample in enumerate(stem_to_meta):
        worksheet.write_string(current_row, LEADING_COLS + s, sample, sample_formats[s])
    current_row += 1
    first_data_row = current_row
    worksheet.write_string(current_row, 0, "TOTAL or MAX")
    worksheet.write_string(current_row, 1, "-")
    worksheet.write_string(current_row, 2, "-")
    worksheet.write_number(
        current_row,
        3,
        max(
            sum(1 for sample in stem_to_meta if (md5, sample) in abundance_by_samples)
            for md5 in md5_to_seq
        ),
    )
    worksheet.write_number(
        current_row,
        4,
        max(
            (
                abundance_by_samples.get((md5, sample), 0)
                for md5 in md5_to_seq
                for sample in stem_to_meta
            ),
            default=0,
        ),
    )
    worksheet.write_number(current_row, 5, sum(md5_abundance.values()))
    for s, sample in enumerate(stem_to_meta):
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
            sum(1 for _ in stem_to_meta if (md5, _) in abundance_by_samples),
            total_abundance,
        ]
        for md5, total_abundance in md5_abundance.items()
    ]
    # Sort on species prediction (with blank last, sorting as tilde);
    # number of samples (decreasing), total abundance (decreasing), md5
    data.sort(key=lambda row: (row[1] if row[1] else "~", -row[3], -row[4], row[0]))
    for md5, sp, seq, md5_in_xxx_samples, total_abundance in data:
        sample_counts = [abundance_by_samples.get((md5, _), 0) for _ in stem_to_meta]
        assert sample_counts, f"{md5} sample counts: {sample_counts!r}"
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

    (stem_to_meta, meta_to_stem, meta_names, group_col,) = load_metadata(
        metadata_file,
        metadata_cols,
        metadata_groups,
        metadata_fieldnames,
        metadata_index,
        ignore_prefixes=ignore_prefixes,
        debug=debug,
    )

    meta_default = tuple([""] * len(meta_names))
    fasta_files = {}
    tsv_files = {}
    for fasta_file, tsv_file in find_paired_files(
        inputs, ".fasta", f".{method}.tsv", ignore_prefixes, debug, strict=True
    ):
        sample = file_to_sample_name(fasta_file)
        if sample not in stem_to_meta:
            if require_metadata:
                continue
            else:
                stem_to_meta[sample] = meta_default
                if meta_default in meta_to_stem:
                    meta_to_stem[meta_default].append(sample)
                else:
                    meta_to_stem[meta_default] = [sample]
                if debug:
                    sys.stderr.write(f"DEBUG: Missing metadata for {sample}\n")
        fasta_files[sample] = fasta_file
        tsv_files[sample] = tsv_file
    assert (
        len(fasta_files) == len(tsv_files) == len(stem_to_meta)
    ), f"{len(fasta_files)} FASTA, {len(tsv_files)} TSV, {len(stem_to_meta)} meta"

    if debug:
        sys.stderr.write(
            f"DEBUG: Have metadata for {len(stem_to_meta)} samples,"
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
    stats_fields = ("Raw FASTQ", "Flash", "Cutadapt", "Threshold")
    try:
        sample_stats = load_fasta_headers(
            fasta_files,
            ("raw_fastq", "flash", "cutadapt", "threshold"),
            -1,
        )
    except KeyError as err:
        sys.stderr.write(
            f"WARNING: Missing header information in FASTA file(s): {err}\n"
        )
        sample_stats = {}
        stats_fields = []

    bad_fields = []
    for i, field in enumerate(stats_fields):
        if -1 in (_[i] for _ in sample_stats.values()):
            bad_fields.append(i)
            sys.stderr.write(f"WARNING: Dropping {field} column as missing stats\n")
        elif field == "Threshold" and len({_[i] for _ in sample_stats.values()}) == 1:
            bad_fields.append(i)
            if debug:
                sys.stderr.write("DEBUG: Dropping abundance threshold as all same\n")
    if bad_fields:
        for sample, value in sample_stats.items():
            sample_stats[sample] = [
                _ for i, _ in enumerate(value) if i not in bad_fields
            ]
        stats_fields = [_ for i, _ in enumerate(stats_fields) if i not in bad_fields]
    del bad_fields

    md5_abundance = Counter()
    abundance_by_samples = {}
    md5_species = {}  # sp values are sets (e.g. ambiguous matches)
    md5_to_seq = {}
    sample_species_counts = {}  # 2nd key is sp list with semi-colons

    if debug:
        sys.stderr.write("Loading FASTA sequences and abundances\n")
    for sample, fasta_file in fasta_files.items():
        assert sample in stem_to_meta, sample
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
        sample_species_counts[sample] = Counter()
        assert sample in stem_to_meta, sample
        # TODO: Look at taxid here?
        for name, _, sp_list in parse_species_tsv(predicted_file, min_abundance):
            md5, abundance = split_read_name_abundance(name)
            if min_abundance > 1 and abundance < min_abundance:
                continue
            assert abundance_by_samples[md5, sample] == abundance, name
            if sp_list:
                md5_species[md5].update(sp_list.split(";"))

            sample_species_counts[sample][sp_list] += abundance_from_read_name(name)

    sample_summary(
        sample_species_counts,
        meta_to_stem,
        stem_to_meta,
        meta_names,
        group_col,
        sample_stats,
        stats_fields,
        show_unsequenced=show_unsequenced,
        output=f"{stem}.samples.{method}.tsv",
        excel=f"{stem}.samples.{method}.xlsx",
        human_output=f"{stem}.samples.{method}.txt",
        method=method,
        min_abundance=min_abundance,
        debug=debug,
    )
    sys.stderr.write(f"Wrote {stem}.samples.{method}.*\n")

    read_summary(
        md5_to_seq,
        md5_species,
        md5_abundance,
        abundance_by_samples,
        stem_to_meta,
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
