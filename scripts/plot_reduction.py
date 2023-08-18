#!/usr/bin/env python
"""Python 3 script plotting a stacked line chart of the reads in each sample.

This Python script was used to generate Figure 4 in the Cock *et al.* (2023)
paper in PeerJ, using matplotlib. This figure replaced the equivalent figure
in the preprint https://doi.org/10.1101/2023.03.24.534090 created with
Microsoft Excel.
"""
import argparse
import sys

import matplotlib as mpl
import matplotlib.pylab as plt
import numpy as np

plt.style.use("tableau-colorblind10")

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.2")
    sys.exit(0)

# Parse Command Line
usage = """\
This Python script uses matplotlib to draw a stacked line chart (with a line
for each sample) from a THAPBI PICT sample report TSV input file. It plots

* Raw FASTQ read count
* Merged read count (after applying flash)
* Primer matched (after applying cutadapt)
* Without singletons (reads seen only once are discard by default)
* Accepted reads (after applying abundance thresholds)
* Unique reads (number of unique reads accepted, not their total)

The final drop off shown, usually from millions to hundreds of sequences, is
to illustrate switching from counting reads to counting unique sequences as a
tally table. See example figure in Cock *et al.* (2023).

Example usage:

$ ./plot_reduction.py -i thapbi-pict.ITS1.samples.onebp.tsv \
                      -o thapbi-pict.ITS1.reduction.pdf -c 1

By default it draws one line per sample, and the column argument can be set
to another column for sample labels. However, this will also merge on that
column allowing you to plot things like biological samples vs controls, or
different sample groups. This is particularly useful when the number of
samples is too large for the legend to be meaningful. You can even group
on multiple columns (e.g. sequencing platform and amplicon name).

The default mode is a stacked line graph of raw read counts, but it can also
draw non-stacked raw counts, or in percentage mode shows each sample (or
group of samples) as percentages of its raw FASTQ read count. In this mode
you may wish to exclude any negative controls (by subsetting your input).
"""

parser = argparse.ArgumentParser(
    prog="plot_reduction.py",
    description="Stacked line plot showing read data reduction per sample.",
    epilog=usage,
)
parser.add_argument(
    "-i",
    "--input",
    default="/dev/stdin",
    metavar="TSV",
    help="THAPBI PICT sample report TSV input file. Required, default stdin.",
)
parser.add_argument(
    "-c",
    "--column",
    type=str,
    default="0",
    metavar="COLUMN",
    help="Which column(s) in the input table to pool on and use as the caption. "
    "Use 0 (default) for the sample FASTQ stem in the 'Sequencing Sample' "
    "column. Use a comma-separated list of integers for a column numbers, or a "
    "single column header name.",
)
parser.add_argument(
    "-m",
    "--mode",
    choices=["stacked", "counts", "percent"],
    default="stacked",
    help="Plot stacked raw counts (default), raw counts, or percentages.",
)
parser.add_argument(
    "-t",
    "--threshold",
    type=int,
    default=10000,
    metavar="COUNT",
    help="For percentage graphs, minimum number of raw reads to include a sample. "
    "Default 10,000 reads is suitable for MiSeq data.",
)
parser.add_argument(
    "-o",
    "--output",
    dest="output",
    metavar="IMAGE",
    default=None,
    help="Output image filename, ending '.png' or '.pdf'. Default live plot.",
)
if len(sys.argv) == 1:
    sys.exit("ERROR: Invalid command line, try -h or --help.")
options = parser.parse_args()


def load_samples(input_sample_report_tsv, caption_column=0, sample_threshold=0):
    """Load a THAPBI PICT TSV sample report."""
    # The key data is all in the headers of the THAPBI PICT tally file but
    # that lacks any user-supplied metadata which we want for sample names.
    data = {}  # key on sample/group caption via specified column
    captions = [
        "Raw FASTQ",
        "Flash",
        "Cutadapt",
        "Without singletons",
        "Accepted reads",
        "Accepted unique",
    ]
    try:
        caption_column = tuple(int(_) for _ in caption_column.split(","))
    except ValueError:
        pass
    with open(input_sample_report_tsv) as handle:
        line = handle.readline()
        if not line.startswith("#"):
            sys.exit("ERROR - Input TSV file did not start with #")
        parts = line[1:].rstrip("\n").split("\t")
        try:
            idn_col = (parts.index("Sequencing sample"),)  # default caption
            raw_col = parts.index("Raw FASTQ")
            merged_col = parts.index("Flash")
            primer_matched_col = parts.index("Cutadapt")
            singletons_col = parts.index("Singletons")
            accepted_total_col = parts.index("Accepted")
            accepted_unique_col = parts.index("Unique")
        except IndexError:
            sys.exit("ERROR - Did not find all expected columns in TSV header")
        if isinstance(caption_column, tuple):
            if caption_column == (0,):
                # Default, use inferred idn_col
                pass
            elif max(caption_column) > len(parts):
                sys.exit(
                    f"ERROR - Only {len(parts)} columns, "
                    f"can't use {max(caption_column)}"
                )
            else:
                idn_col = tuple(v - 1 for v in caption_column)
        elif caption_column in parts:
            idn_col = (parts.index(caption_column),)
        else:
            sys.exit(f"ERROR - Did not find this in header columns: {caption_column}")
        sys.stderr.write(
            f"Using column(s) {','.join(str(v+1) for v in idn_col)}, "
            f"{' - '.join(parts[v] for v in idn_col)}, for captions/grouping\n"
        )
        count = 0
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if int(parts[raw_col]) < sample_threshold:
                sys.stderr.write(
                    "WARNING: Ignoring low abundance sample "
                    f"{' - '.join(parts[v] for v in idn_col)}\n"
                )
                continue
            try:
                idn = tuple(parts[v] for v in idn_col)
                if idn in data:
                    # Append - TODO refactor, maybe use numpy array?
                    old = data[idn]
                    data[idn] = (
                        old[0] + int(parts[raw_col]),
                        old[1] + int(parts[merged_col]),
                        old[2] + int(parts[primer_matched_col]),
                        old[3]
                        + int(parts[primer_matched_col])
                        - int(parts[singletons_col]),
                        old[4] + int(parts[accepted_total_col]),
                        old[5] + int(parts[accepted_unique_col]),
                    )
                    del old
                else:
                    data[idn] = (
                        int(parts[raw_col]),
                        int(parts[merged_col]),
                        int(parts[primer_matched_col]),
                        int(parts[primer_matched_col]) - int(parts[singletons_col]),
                        int(parts[accepted_total_col]),
                        int(parts[accepted_unique_col]),
                    )
            except IndexError:
                sys.exit("ERROR - Not enough fields in line:\n" + repr(line))
            count += 1
        if len(data) < count:
            sys.stderr.write(
                f"Loaded read counts for {count} samples into {len(data)} groups\n"
            )
        else:
            sys.stderr.write(f"Loaded read counts for {len(data)} samples\n")
    return captions, [" - ".join(_) for _ in data.keys()], list(data.values())


def plot_read_reduction(
    input_sample_report_tsv,
    output_stacked_plot,
    caption_column=0,
    mode="stacked",
    percent_threshold=10000,
):
    """Load a THAPBI PICT TSV sample report, and plot read reduction."""
    # The key data is all in the headers of the THAPBI PICT tally file but
    # that lacks any user-supplied metadata which we want for sample names.
    captions, labels, data = load_samples(
        input_sample_report_tsv,
        caption_column,
        percent_threshold if mode == "percent" else 0,
    )

    color_count = len(plt.rcParams["axes.prop_cycle"].by_key()["color"])
    line_styles = ("solid", "dotted", "dashed", "dashdot")
    marker_styles = (".", "o", "s")

    fig, ax = plt.subplots(figsize=(12, 6))
    # ax.stackplot(captions, data, labels=labels)
    if mode == "stacked":
        line_values = np.zeros(len(captions), dtype=np.uint64)
    for idx, (sample, values) in enumerate(zip(labels, data)):
        if mode == "percent":
            # Convert to a percentage of this group's raw reads
            line_values = 100.0 * np.array(values, dtype=np.float64) / values[0]
        elif mode == "counts":
            line_values = np.array(values, dtype=np.uint64)
        elif mode == "stacked":
            # Add to previous values for a stacked read-count plot
            line_values += np.array(values, dtype=np.uint64)
        else:
            sys.exit(f"ERROR - Unsupported mode {mode}")
        ax.plot(
            captions,
            line_values,
            label=sample,
            linestyle=line_styles[(idx // color_count) % len(line_styles)],
            linewidth=2,
            marker=marker_styles[(idx // color_count) % len(marker_styles)],
        )
    # ax.set_ylim(0)
    ax.set_xlim(-0.1, 5.1)
    ax.xaxis.set_ticks_position("top")
    legend_cols = 1 + len(labels) // 28
    if mode == "percent":
        # Sub-plot rectangle (left, bottom, right, top),
        # generous space on the right for legens with long sample names:
        plt.tight_layout(rect=[0, 0, 0.85, 1])
        ax.set_ylim(0, 100 + 5)  # top margin for x-axis labels
        ax.yaxis.set_ticks(range(0, 110, 10))  # every 10%
        ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter())
        # Generous space on the right for legens with long sample names:

    else:
        # Leaving a little extra space on left if read counts are excessive,
        # and generous space on the right for legends with long sample names:
        plt.tight_layout(rect=[0.05, 0, 0.85, 1])
        # Force read counts to be comma-separated thousands (not scientific notation)
        # ax.yaxis.get_major_formatter().set_useOffset(False)
        # ax.yaxis.get_major_formatter().set_scientific(False)
        ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter("{x:,.0f}"))
    ax.legend(
        reverse=True,
        loc="center left",
        # Left position affect by layout and xlim?
        bbox_to_anchor=(1, 0.5),
        ncol=legend_cols,
    )
    ax.set_frame_on(False)
    ax.grid(axis="y", which="major")

    # Display
    if output_stacked_plot:
        plt.savefig(output_stacked_plot, dpi=300, orientation="landscape")
        sys.stderr.write(f"Drew {output_stacked_plot} from {input_sample_report_tsv}\n")
    else:
        # Interactive
        sys.stderr.write(f"Displaying image from {input_sample_report_tsv}\n")
        plt.show()


plot_read_reduction(
    options.input, options.output, options.column, options.mode, options.threshold
)
