# Copyright 2019-2020 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Summarise classification results at sample and read level.

This implements the ``thapbi_pict summary ...`` command.
"""
import os
import sys
from collections import Counter

from Bio.SeqIO.FastaIO import SimpleFastaParser

from .utils import abundance_from_read_name
from .utils import color_bands
from .utils import find_requested_files
from .utils import load_metadata
from .utils import parse_species_tsv
from .utils import sample_sort
from .utils import split_read_name_abundance


def file_to_sample_name(filename):
    """Given filename (without and directory name), return sample name only.

    i.e. XXX.fasta --> and XXX.method.tsv --> XXX
    """
    if filename.endswith(".fasta"):
        return os.path.basename(filename).rsplit(".", 1)[0]
    elif filename.endswith(".tsv"):
        return os.path.basename(filename).rsplit(".", 2)[0]
    else:
        raise ValueError(f"Invalid file_to_sample_name arg: {filename}")


def sample_summary(
    tsv_files,
    metadata,
    meta_names,
    group_col,
    output,
    excel,
    human_output,
    method,
    min_abundance=1,
    ignore_prefixes=None,
    debug=False,
):
    """Implement the ``thapbi_pict sample-summary`` command.

    The expectation is that the inputs represent all the samples from
    a meaningful group, likely from multiple sequencing runs (plates).
    """
    if not (output or human_output):
        sys.exit("ERROR: No output file specified.\n")

    samples = set()
    missing_meta = set()
    genus_predictions = set()
    sample_genus_counts = {}
    species_predictions = set()  # includes A;B;C ambiguous entries
    sample_species_counts = {}
    for sample, predicted_file in tsv_files.items():
        if sample not in metadata:
            missing_meta.add(sample)
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
            f"DEBUG: {len(samples):d} samples with predictions for"
            f" {len(genus_predictions):d} genera\n"
        )

    samples = sample_sort(samples)

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
                if not sample_batch:
                    human.write("Has not been sequenced.\n\n")
            except BrokenPipeError:
                human = None
        # Now do the samples in this batch
        for sample in sample_batch:
            if sample not in samples:
                if not sample.startswith(ignore_prefixes):
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
                    try:
                        cell_format = sample_formats[current_row]  # sample number
                    except IndexError:
                        cell_format = None
                    current_row += 1
                    assert len(meta_names) == len(metadata)
                    for offset, value in enumerate(metadata):
                        worksheet.write_string(current_row, offset, value, cell_format)
                    assert col_offset == len(meta_names)
                    worksheet.write_string(current_row, col_offset, sample, cell_format)
                    worksheet.write_number(
                        current_row,
                        col_offset + 1,
                        sum(sample_species_counts[sample].values()),
                        cell_format,
                    )
                    for offset, genus in enumerate(genus_predictions):
                        worksheet.write_number(
                            current_row,
                            col_offset + 2 + offset,
                            sample_genus_counts[sample][genus],
                            cell_format,
                        )
                    for offset, sp in enumerate(species_columns):
                        worksheet.write_number(
                            current_row,
                            col_offset + 2 + len(genus_predictions) + offset,
                            sample_species_counts[sample][sp],
                            cell_format,
                        )

            if human:
                all_sp = set()
                unambig_sp = set()
                if sample in sample_species_counts:
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


def read_summary(
    fasta_files,
    tsv_files,
    metadata,
    meta_names,
    group_col,
    output,
    method,
    min_abundance=1,
    excel=None,
    ignore_prefixes=None,
    debug=False,
):
    """Implement the ``thapbi_pict plate-summary`` command.

    The expectation is that the inputs represent all the samples
    from one (96 well) plate, or some other meaningful batch.
    """
    if not output:
        sys.exit("ERROR: No output file specified.\n")

    samples = set()
    md5_abundance = Counter()
    abundance_by_samples = {}
    md5_species = {}
    md5_to_seq = {}

    if debug:
        sys.stderr.write("Loading FASTA sequences and abundances\n")
    samples = set()
    for sample, fasta_file in fasta_files.items():
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
    if len(samples) < len(metadata):
        missing = set(metadata).difference(samples)
        sys.stderr.write(
            f"WARNING: {len(missing)} samples in metadata have not been sequenced,"
            f" {sample_sort(missing)[0]} etc.\n"
        )
    samples = sample_sort(samples)

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
    meta_default = [""] * len(meta_names)
    sample_formats = [None] * len(samples)
    if metadata:
        # Insert extra header rows at start for sample meta-data
        # Make a single metadata call for each sample
        meta = [metadata.get(sample, meta_default) for sample in samples]
        for i, name in enumerate(meta_names):
            handle.write(
                "#%s%s\t%s\n"
                % ("\t" * (LEADING_COLS - 1), name, "\t".join(_[i] for _ in meta))
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
        "\tMax-sample-abundance\tTotal-abundance\t%s\n" % (method, "\t".join(samples))
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
        worksheet.write_string(current_row, 1, method + "-predictions")
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
    for filename in find_requested_files(
        [_ for _ in inputs if not _.endswith(".tsv")],
        ".fasta",
        ignore_prefixes,
        debug=debug,
    ):
        sample = file_to_sample_name(filename)
        if require_metadata and sample not in metadata:
            continue
        elif sample in fasta_files:
            sys.exit(f"ERROR: Multiple FASTA files using {sample} naming")
        fasta_files[sample] = filename

    tsv_files = {}
    for filename in find_requested_files(
        [_ for _ in inputs if not _.endswith(".fasta")],
        f".{method}.tsv",
        ignore_prefixes,
        debug,
    ):
        sample = file_to_sample_name(filename)
        if require_metadata and sample not in metadata:
            continue
        elif sample in tsv_files:
            sys.exit(f"ERROR: Multiple TSV files using {sample} naming")
        elif sample not in fasta_files:
            sys.exit(f"ERROR: {filename} without {sample} FASTA file")
        tsv_files[sample] = filename

    if debug:
        sys.stderr.write(
            f"Found {len(fasta_files)} FASTA files, and"
            f" {len(tsv_files)} sample predictions using method {method}\n"
        )
    if set(fasta_files) != set(tsv_files):
        sys.exit("ERROR: FASTA vs TSV sample name mismatch")
    if not tsv_files:
        sys.exit("ERROR: No input FASTA and TSV files found")

    if report:
        stem = os.path.join(out_dir, report)
    else:
        # Include version number here?
        stem = os.path.join(out_dir, "thapbi-pict")

    return_code = sample_summary(
        tsv_files,
        metadata,
        meta_names,
        group_col,
        output=f"{stem}.samples.{method}.tsv",
        excel=f"{stem}.samples.{method}.xlsx",
        human_output=f"{stem}.samples.{method}.txt",
        method=method,
        min_abundance=min_abundance,
        ignore_prefixes=ignore_prefixes,
        debug=debug,
    )
    if return_code:
        return return_code
    sys.stderr.write(f"Wrote {stem}.samples.{method}.*\n")

    return_code = read_summary(
        fasta_files,
        tsv_files,
        metadata,
        meta_names,
        group_col,
        output=f"{stem}.reads.{method}.tsv",
        excel=f"{stem}.reads.{method}.xlsx",
        method=method,
        min_abundance=min_abundance,
        ignore_prefixes=ignore_prefixes,
        debug=debug,
    )
    if return_code:
        return return_code
    sys.stderr.write(f"Wrote {stem}.reads.{method}.*\n")

    return return_code
