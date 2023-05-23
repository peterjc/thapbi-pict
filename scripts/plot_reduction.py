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
    print("v0.0.1")
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
    help="Which column in the input table to use as the caption. Use 0 "
    "(default) for the sample FASTQ stem in the 'Sequencing Sample' column. "
    "Use integer 1 or more for a column number, or the column header name.",
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


def load_samples(input_sample_report_tsv, caption_column=0):
    """Load a THAPBI PICT TSV sample report."""
    # The key data is all in the headers of the THAPBI PICT tally file but
    # that lacks any user-supplied metadata which we want for sample names.
    data = []
    labels = []
    captions = [
        "Raw FASTQ",
        "Flash",
        "Cutadapt",
        "Without singletons",
        "Accepted reads",
        "Accepted unique",
    ]
    try:
        caption_column = int(caption_column)
    except ValueError:
        pass
    with open(input_sample_report_tsv) as handle:
        line = handle.readline()
        if not line.startswith("#"):
            sys.exit("ERROR - Input TSV file did not start with #")
        parts = line[1:].rstrip("\n").split("\t")
        try:
            idn_col = parts.index("Sequencing sample")  # default caption
            raw_col = parts.index("Raw FASTQ")
            merged_col = parts.index("Flash")
            primer_matched_col = parts.index("Cutadapt")
            singletons_col = parts.index("Singletons")
            accepted_total_col = parts.index("Accepted")
            accepted_unique_col = parts.index("Unique")
        except IndexError:
            sys.exit("ERROR - Did not find all expected columns in TSV header")
        if isinstance(caption_column, int):
            if caption_column == 0:
                # Default, use inferred idn_col
                pass
            elif caption_column > len(parts):
                sys.exit(
                    f"ERROR - Only {len(parts)} columns, can't use {caption_column}"
                )
            else:
                if caption_column > idn_col:
                    sys.stderr.write(
                        f"WARNING - Selected caption column {parts[caption_column-1]} "
                        "is not metadata\n"
                    )
                idn_col = caption_column - 1
        elif caption_column in parts:
            idn_col = parts.index(caption_column)
        else:
            sys.exit(f"ERROR - Did not find this in header columns: {caption_column}")
        sys.stderr.write(
            f"Using column {idn_col+1}, '{parts[idn_col]}', for sample captions\n"
        )
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            try:
                labels.append(parts[idn_col])
                data.append(
                    (
                        int(parts[raw_col]),
                        int(parts[merged_col]),
                        int(parts[primer_matched_col]),
                        int(parts[primer_matched_col]) - int(parts[singletons_col]),
                        int(parts[accepted_total_col]),
                        int(parts[accepted_unique_col]),
                    )
                )
            except IndexError:
                sys.exit("ERROR - Not enough fields in line:\n" + repr(line))
        sys.stderr.write(f"Loaded read counts for {len(data)} samples\n")
    return captions, labels, data


def plot_read_reduction(input_sample_report_tsv, output_stacked_plot, caption_column=0):
    """Load a THAPBI PICT TSV sample report, and plot read reduction."""
    # The key data is all in the headers of the THAPBI PICT tally file but
    # that lacks any user-supplied metadata which we want for sample names.
    captions, labels, data = load_samples(input_sample_report_tsv, caption_column)

    color_count = len(plt.rcParams["axes.prop_cycle"].by_key()["color"])
    line_styles = ("solid", "dotted", "dashed", "dashdot")
    marker_styles = (".", "o", "s")

    fig, ax = plt.subplots(figsize=(12, 6))
    # ax.stackplot(captions, data, labels=labels)
    stacked = np.zeros(len(captions), dtype=np.uint64)
    for idx, (sample, values) in enumerate(zip(labels, data)):
        stacked += np.array(values, dtype=np.uint64)
        ax.plot(
            captions,
            stacked,
            label=sample,
            linestyle=line_styles[(idx // color_count) % len(line_styles)],
            linewidth=2,
            marker=marker_styles[(idx // color_count) % len(marker_styles)],
        )
    # Leaving a little extra space on left if read counts are excessive,
    # and generous space on the right for legends with long sample names:
    plt.tight_layout(rect=[0.05, 0, 0.85, 1])
    # ax.set_ylim(0)
    ax.set_xlim(-0.1, 5.1)
    ax.xaxis.set_ticks_position("top")
    ax.legend(
        reverse=True,
        loc="center left",
        # Left position affect by layout and xlim?
        bbox_to_anchor=(1, 0.5),
    )
    ax.set_frame_on(False)
    ax.grid(axis="y", which="major")
    # Force read counts to be comma-separated thousands (not scientific notation)
    # ax.yaxis.get_major_formatter().set_useOffset(False)
    # ax.yaxis.get_major_formatter().set_scientific(False)
    ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter("{x:,.0f}"))

    # Display
    if output_stacked_plot:
        plt.savefig(output_stacked_plot, dpi=300, orientation="landscape")
        sys.stderr.write(f"Drew {output_stacked_plot} from {input_sample_report_tsv}\n")
    else:
        # Interactive
        sys.stderr.write(f"Displaying image from {input_sample_report_tsv}\n")
        plt.show()


plot_read_reduction(options.input, options.output, options.column)
