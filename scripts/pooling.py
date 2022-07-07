#!/usr/bin/env python3
"""Pool THAPBI PICT sample report using metadata."""
import argparse
import sys

import numpy as np
import xlsxwriter

# from collections import Counter

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.3")
    sys.exit(0)

# Parse Command Line
usage = """\
The input file should be a THAPBI PICT sample summary report as
a plain text tab separated variable (TSV) file. Give a filename
stem as the output argument (will make TSV and Excel files). The
column argument is which metadata to retain and group by. e.g.

$ python pooling.py -i input.tsv -o pooled -c 1,2,5
"""

parser = argparse.ArgumentParser(
    prog="pooling.py",
    description="Pooling the THAPBI PICT sample report (e.g. by replicate).",
    epilog=usage,
)
parser.add_argument(
    "-i",
    "--input",
    dest="input",
    default="/dev/stdin",
    metavar="FILE",
    help="Input TSV filename, default stdin.",
)
parser.add_argument(
    "-c",
    "--columns",
    type=str,
    default="1",
    metavar="COLUMNS",
    help="Comma separated list (e.g, '1,3,5') of metadata columns from the "
    "input table to retain and pool on.",
)
parser.add_argument(
    "-p",
    "--pending",
    type=str,
    default="0",
    metavar="COLUMN",
    help="Columns from the input table containing Y, Yes, True (in any case) "
    "where further sequencing is pending. Such entries will get a row of ??? "
    "(in addition to a row of values if any were non-zero).",
)
parser.add_argument(
    "-g",
    "--first-genus",
    type=str,
    default="",
    metavar="GENUS-LIST",
    help="Optional semi-colon separated list of genera to list first "
    "(overriding existing column order). e.g. Unknown;Phytophthora",
)
parser.add_argument(
    "-G",
    "--last-genus",
    type=str,
    default="",
    metavar="GENUS-LIST",
    help="Optional semi-colon separated list of genera to list last "
    "(overriding existing column order). e.g. Synthetic;Unknown",
)
parser.add_argument(
    "-b",
    "--boolean",
    default=False,
    action="store_true",
    help="Replace accepted read counts with boolean (Y and N entries).",
)
parser.add_argument(
    "-z",
    "--hide-zeros",
    default=False,
    action="store_true",
    help="Hide zero count columns (e.g. from using filtered input data).",
)
parser.add_argument(
    "--pcr",
    default=False,
    action="store_true",
    help="Replace sequence sample count with inferred PCR status, recommend "
    "using with the -p / --pending setting otherwise no Illumina data is "
    "taken to imply PCR was negative.",
)
parser.add_argument(
    "-o",
    "--output",
    dest="output",
    default="-",
    metavar="STEM",
    help="Output filename stem (defaults to '-' meaning TSV only to stdout)",
)

if len(sys.argv) == 1:
    sys.exit("ERROR: Invalid command line, try -h or --help.")
options = parser.parse_args()


def pool(
    input_filename,
    output_stem,
    columns_str,
    column_pending,
    show_boolean,
    hide_zeros,
    first_genus,
    last_genus,
    pcr_status,
):
    """Pool samples to make a more consise report."""
    try:
        value_cols = [int(_) - 1 for _ in columns_str.split(",")]
    except ValueError:
        sys.exit(
            "ERROR: Columns should be a comma separated list"
            f" of positive integers, not {columns_str!r}."
        )
    if min(value_cols) < 0:
        sys.exit("ERROR: Invalid column, should all be positive.")
    try:
        column_pending = int(column_pending) - 1
    except ValueError:
        sys.exit(
            "ERROR: Pending column should be a positive integer (or zero for none)."
        )
    if column_pending < 0:
        column_pending = None

    meta_samples = {}
    meta_species = {}
    meta_pending = set()

    with open(input_filename) as handle:
        line = handle.readline().rstrip("\n")
        if not line.startswith("#") or "\t" not in line:
            sys.exit("ERROR: Invalid TSV input file.")
        header = [_.rstrip() for _ in line.split("\t")]
        try:
            sample_col = header.index("Sequencing sample")
        except ValueError:
            sys.exit("ERROR: Header does not match THAPBI PICT sample report.")
        try:
            count_col = header.index("Accepted")
            assert count_col + 1 == header.index("Unique")
            first_sp_col = count_col + 2
        except ValueError:
            try:
                # Prior to v0.12.0 this was called "Read count",
                count_col = header.index("Read count")
                first_sp_col = count_col + 1
            except ValueError:
                sys.exit("ERROR: Header does not match THAPBI PICT sample report.")
        if max(value_cols) >= min(sample_col, count_col):
            sys.exit(
                f"ERROR: Requested column {max(value_cols)+1} not in metadata range."
            )
        if column_pending is not None and column_pending >= min(sample_col, count_col):
            sys.exit("ERROR: Pending column not in metadata range.")
        sp_headers = header[first_sp_col:]
        sp_null = ["-"] * len(sp_headers)
        meta_headers = [header[_] for _ in value_cols]

        for line in handle:
            parts = [_.rstrip() for _ in line.rstrip("\n").split("\t")]
            if len(parts) != len(header):
                sys.exit("ERROR: Inconsistent field counts")
            meta = tuple(parts[_] for _ in value_cols)
            if column_pending is not None and parts[column_pending].upper().strip() in (
                "Y",
                "YES",
                "TRUE",
            ):
                meta_pending.add(meta)
            samples = [_ for _ in parts[sample_col].split(";") if _ != "-"]

            sp_counts = parts[first_sp_col:]
            if sp_counts == sp_null:
                # Was not sequenced
                if meta in meta_samples:
                    meta_samples[meta].update(samples)
                    # Don't record replace with a None value
                else:
                    meta_samples[meta] = set(samples)
                    meta_species[meta] = None
            else:
                assert "-" not in sp_counts, sp_counts
                sp_counts = np.array([int(_) for _ in sp_counts], int)
                if meta in meta_samples:
                    meta_samples[meta].update(samples)
                    if meta_species[meta] is None:
                        # Replace the None value
                        meta_species[meta] = sp_counts
                    else:
                        meta_species[meta] += sp_counts
                else:
                    meta_samples[meta] = set(samples)
                    meta_species[meta] = sp_counts

    total_counts = sum(v for v in meta_species.values() if v is not None)
    assert len(total_counts) == len(sp_headers), total_counts.shape
    if hide_zeros:
        sp_headers = [v for v, t in zip(sp_headers, total_counts) if t]
        for meta in meta_species:
            if meta_species[meta] is not None:
                meta_species[meta] = [
                    v for v, t in zip(meta_species[meta], total_counts) if t
                ]
        total_counts = [t for t in total_counts if t]

    if first_genus:
        # Need to resort output columns: sp_headers, meta_species, total_counts
        first_genus = [_.strip() for _ in first_genus.replace(",", ";").split(";")]
        last_genus = [
            _.strip()
            for _ in last_genus.replace(",", ";").split(";")
            if _ not in first_genus
        ]
        g = [_.split(" ")[0] for _ in sp_headers]
        a = [i for i, v in enumerate(g) if v in first_genus]
        b = [i for i, v in enumerate(g) if v not in (first_genus + last_genus)]
        c = [i for i, v in enumerate(g) if v in last_genus]
        new_order = a + b + c
        del a, b, c, g
        sp_headers = [sp_headers[i] for i in new_order]
        total_counts = [total_counts[i] for i in new_order]
        for meta, old_counts in list(meta_species.items()):
            if old_counts is not None:
                meta_species[meta] = [old_counts[i] for i in new_order]
        del old_counts, new_order

    if show_boolean:
        # Convert counts to Y/N
        def display(count):
            return "Y" if count else "N"

    else:
        # Use counts
        display = str

    if output_stem == "-":
        handle = sys.stdout
    else:
        handle = open(output_stem + ".tsv", "w")

        # Excel setup
        # -----------
        workbook = xlsxwriter.Workbook(output_stem + ".xlsx")
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
        central_text_format = workbook.add_format({"align": "center"})

        # Set first row to be tall, with vertical text
        current_row = 0
        col_offset = len(meta_headers)

        worksheet.set_row(0, 150, vertical_text_format)
        for offset, name in enumerate(meta_headers):
            worksheet.write_string(current_row, offset, name)
        worksheet.write_string(
            current_row, col_offset, "PCR status" if pcr_status else "Samples-sequenced"
        )
        for offset, name in enumerate(sp_headers):
            worksheet.write_string(current_row, col_offset + 1 + offset, name)
        worksheet.freeze_panes(current_row + 1, col_offset + 1)

    # TSV setup
    # ---------
    handle.write(
        "\t".join(meta_headers)
        + ("\tPCR status\t" if pcr_status else "\tSamples-sequenced\t")
        + "\t".join(sp_headers)
        + "\n"
    )

    for meta, sp_counts in meta_species.items():
        # Cases: pending (True/False), data (positive/zero, missing)
        # If pending and missing, one line of output only!
        #
        #           | Pending      | Not Pending
        # ----------*--------------*----------------
        # Counts    | Counts & ??? | One counts line
        # No counts | One ??? line | One --- line

        if sp_counts is not None:
            # Line of values for unsequenced data
            # Might have "???" pending line as well!
            if pcr_status:
                sample_status = "Positive"
            else:
                sample_status = str(len(meta_samples[meta]))
            handle.write(
                "\t".join(meta)
                + "\t"
                + sample_status
                + "\t"
                + "\t".join(display(_) for _ in sp_counts)
                + "\n"
            )
            if output_stem != "-":
                current_row += 1
                for offset, value in enumerate(meta):
                    worksheet.write_string(current_row, offset, value)
                if pcr_status:
                    worksheet.write_string(current_row, col_offset, sample_status)
                else:
                    # Special cased to write a number not a string:
                    assert sample_status == str(len(meta_samples[meta]))
                    worksheet.write_number(
                        current_row, col_offset, len(meta_samples[meta])
                    )
                if show_boolean:
                    for offset, value in enumerate(sp_counts):
                        worksheet.write_string(
                            current_row,
                            col_offset + 1 + offset,
                            display(value),
                            central_text_format,
                        )
                else:
                    for offset, value in enumerate(sp_counts):
                        worksheet.write_number(
                            current_row, col_offset + 1 + offset, value
                        )
        if meta in meta_pending:
            # Line of "?" for unsequenced but pending data
            if pcr_status:
                sample_status = "Positive (NS)"
            else:
                sample_status = str(len(meta_samples[meta]))
            handle.write(
                "\t".join(meta) + "\t" + sample_status + "\t?" * len(sp_headers) + "\n"
            )
            if output_stem != "-":
                current_row += 1
                for offset, value in enumerate(meta):
                    worksheet.write_string(current_row, offset, value)
                if pcr_status:
                    worksheet.write_string(current_row, col_offset, sample_status)
                else:
                    # Special cased to write a number not a string:
                    assert sample_status == "0"
                    worksheet.write_number(current_row, col_offset, 0)
                for offset in range(len(sp_headers)):
                    worksheet.write_string(
                        current_row, col_offset + 1 + offset, "?", central_text_format
                    )
        elif sp_counts is None:
            # Line of "-" for unsequenced data
            if pcr_status:
                sample_status = "Negative"
            else:
                assert len(meta_samples[meta]) == 0
                sample_status = "0"
            handle.write(
                "\t".join(meta) + "\t" + sample_status + "\t-" * len(sp_headers) + "\n"
            )
            if output_stem != "-":
                current_row += 1
                for offset, value in enumerate(meta):
                    worksheet.write_string(current_row, offset, value)
                if pcr_status:
                    worksheet.write_string(current_row, col_offset, sample_status)
                else:
                    # Special cased to write a number not a string:
                    assert sample_status == "0"
                    worksheet.write_number(current_row, col_offset, 0)
                for offset in range(len(sp_headers)):
                    worksheet.write_string(
                        current_row, col_offset + 1 + offset, "-", central_text_format
                    )

    if output_stem != "-":
        handle.close()
        # Excel polish...
        if pcr_status:
            # Can set this just right
            worksheet.set_column(col_offset, col_offset, 9.7)
        # Defined first, but takes priority over later conditional rules:
        worksheet.conditional_format(
            1,
            col_offset + 1,
            current_row,
            col_offset + 1 + len(sp_headers),
            {
                "type": "cell",
                "criteria": "equal to",
                "value": '"-"',
                "format": grey_conditional_format,
            },
        )
        if show_boolean:
            # With just single characters only need narrow columns:
            worksheet.set_column(col_offset + 1, col_offset + 1 + len(sp_headers), 2.7)
            worksheet.conditional_format(
                1,
                col_offset + 1,
                current_row,
                col_offset + 1 + len(sp_headers),
                {
                    "type": "cell",
                    "criteria": "equal to",
                    "value": '"Y"',
                    "format": red_conditional_format,
                },
            )
        else:
            worksheet.conditional_format(
                1,
                col_offset + 1,
                current_row,
                col_offset + 1 + len(sp_headers),
                {
                    "type": "cell",
                    "criteria": "greater than",
                    "value": 0,
                    "format": red_conditional_format,
                },
            )
        workbook.close()


pool(
    options.input,
    options.output,
    options.columns,
    options.pending,
    options.boolean,
    options.hide_zeros,
    options.first_genus,
    options.last_genus,
    options.pcr,
)
