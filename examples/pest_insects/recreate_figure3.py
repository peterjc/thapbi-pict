#!/usr/bin/env python3
"""Script to reproduce Batovska *et al.* (2021) Figure 3.

See ``figure3.R`` and its output ``figure3.tsv`` which exports the ratios
underlying Batovska *et al.* (2021) Figure 3 from their R-based analysis:
https://github.com/alexpiper/HemipteraMetabarcodingMS/blob/master/hemiptera_metabarcoding.Rmd

This Python script loads ``figure3.tsv`` (original data) or an equivalent file
(our data), and recreates Figure 3 using matplot lib.

Figure 3 consists of five sub-figures:

* Expected (i.e. what was put into the control mixes)
* All genes (i.e. pooling the markers, mean ratios)
* COI
* 18S
* 12S

Each of these five sub-figures is a stacked bar chart breaking down by species
using the following colour key:

* Dark blue: Acizzia alternata/solanicola, aggregated for display
* Dark red: Bactericera cockerelli (TPP)
* Orange: Diuraphis noxia (RWA)
* Sea green: Metopolophium dirhodum
* Light green: Rhopalosiphum padi

The stacked-bars are horizontal, relative abundance labelled 0, 0.25, 0.50,
0.75 and 1. From top top bottom we have four entries (100, 250, 500 and 1000)
for pools 1 to 5:

* 100 Pool 1
* 250 Pool 1
* 500 Pool 1
* 1000 Pool 1
* ...
* 100 Pool 5
* 250 Pool 5
* 500 Pool 5
* 1000 Pool 5

Additionally some of the control samples are asterisked where all species are
correctly identified, and coloured (+) and (-) entries have been added in some
of the bars for false positives and false negatives respectively.
"""
import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

species_key = [
    # Caption, colour
    ("Rhopalosiphum padi", "#95cf77"),
    ("Metopolophium dirhodum", "#3a9e82"),
    ("Diuraphis noxia (RWA)", "#fa6e24"),
    ("Bactericera cockerelli (TPP)", "#ae0707"),
    ("Acizzia alternata/solanicola", "#0c4687"),
]

col_captions = ["Expected", "All Genes", "COI", "18S", "12S"]

# First entry appears at bottom of the plot! (note the reversal in [::-1])
row_captions = [
    f"{count} Pool {pool}"
    for pool in (1, 2, 3, 4, 5)
    for count in (100, 250, 500, 1000)
][::-1]

# Parse Command Line
usage = """\
The input file should be a tab-separated-variable plain text table, for
example 'figure3.tsv' generated from 'figure3.R' containing the data from
Batovska et al. (2021) Figure 3.
"""

parser = argparse.ArgumentParser(
    prog="recreate_figure3.py",
    description="Draw an image like Batovska et al. (2021) Figure 3.",
    epilog=usage,
)
parser.add_argument(
    "-i",
    "--input",
    default="/dev/stdin",
    metavar="TSV",
    help="Input data TSV filename, default stdin. Expects seven columns, "
    "1. Sample caption (for figure), 2. Genus (expanded to hard coded species "
    "names), 3. Expected ratio, 4. COI ratio, 5. 18S ratio, 6. 12S ratio. "
    "Counts may be used inplace of ratios, and will be summed over each sample.",
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


def load_subplot_data(tsv_filename, classifications, classes):
    """Return a dict of dicts of lists of numbers from the TSV file.

    Expects an input TSV file with rows as follows:

    * Sample caption
    * Classification, entry in the classes argument (e.g. taxonomy)
    * Columns of values for each classifier class (e.g. amplicon name)

    Returns dict[class-name][sample-caption] = value
    """
    # print(f"Loading {classes} from {tsv_filename}")
    answer = {}
    with open(tsv_filename) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        assert len(header) >= 3, "Not enough cols"
        answer = {class_name: {} for class_name in header[2:]}
        for class_name in answer:
            assert class_name in classes, f"Did not expect {class_name}"
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            caption = fields[0].replace(" ", "-")  # TODO drop this
            classification = classifications.index(fields[1])
            # print(f"{caption} {fields[2]} --> field {classification}")
            for class_name, value in zip(header[2:], fields[2:]):
                if caption not in answer[class_name]:
                    answer[class_name][caption] = np.zeros(len(classes))
                answer[class_name][caption][classification] = float(value)
    return answer


def plot(tsv_filename, img_filename=None):
    """Recreate Batovska *et al.* (2021) Figure 3."""
    # Pass list of genus names to match figure3.tsv column headers
    subplot_data = load_subplot_data(
        tsv_filename, [_[0].split(None, 1)[0] for _ in species_key], col_captions
    )

    subplot_data["All Genes"] = {
        # Can take ratio mean or sum, either works for this plotting style
        sample: sum(subplot_data[gene][sample] for gene in ("COI", "18S", "12S")) / 3
        for sample in subplot_data["Expected"]
    }

    fig, axes = plt.subplots(1, len(subplot_data))
    y_pos = np.arange(len(row_captions))
    for plot_index, (ax, caption) in enumerate(zip(axes, col_captions)):
        if plot_index == 0:
            assert caption == "Expected"
            ax.set_yticks(y_pos, labels=row_captions)
            ax.set_ylabel("Mock Community")
        else:
            # print(caption, data_dict)
            # ax.set_yticks(y_pos, labels=[""] * len(people))
            ax.set_yticks([])

        assert caption in subplot_data, f"{caption} not in {list(subplot_data.keys())}"
        marker_data = np.array(
            [subplot_data[caption][_.replace(" ", "-")] for _ in row_captions]
        )
        for sp_index, (sp_caption, sp_color) in enumerate(species_key):
            ax.barh(
                y_pos,
                marker_data[:, sp_index],
                # Stacked bar chart by offsetting left position:
                left=marker_data[:, 0:sp_index].sum(axis=1),
                # xerr=error,
                align="center",
                label=sp_caption,
                color=sp_color,
            )
            for i, sample in enumerate(row_captions):
                if (
                    subplot_data["Expected"][sample.replace(" ", "-")][sp_index]
                    and not marker_data[i, sp_index]
                ):
                    # False negative!
                    ax.plot(
                        0.95 - 0.08 * sp_index,
                        y_pos[i],
                        marker="$×$",
                        markersize=7,
                        color=sp_color,
                    )
                if (
                    not subplot_data["Expected"][sample.replace(" ", "-")][sp_index]
                    and marker_data[i, sp_index]
                ):
                    # False positive!
                    ax.plot(
                        0.5 - 0.02 * sp_index,
                        y_pos[i],
                        marker="$!$",
                        markersize=7,
                        color=sp_color,
                    )
                elif 0 < marker_data[i, sp_index] < 0.05:
                    # Highlight low but positive entries
                    ax.plot(
                        0.1 + 0.06 * sp_index,
                        y_pos[i],
                        marker=r"$\checkmark$",
                        markersize=7,
                        color=sp_color,
                    )
        # ax.invert_yaxis()  # labels read top-to-bottom

        # Allocate small visual margin:
        ax.set_ylim(-0.5, len(row_captions) - 0.5)
        ax.set_xlim(-0.025, 1.025)

        ax.set_xticks(
            [0, 0.25, 0.5, 0.75, 1],
            # Using "0.50" looks too squeezed
            labels=["0", "0.25", "0.5", "0.75", "1"],
            fontsize=6,
        )
        if plot_index == 2:
            # Once only on central of five subplots:
            ax.set_xlabel("Relative abundance")
        else:
            ax.set_xlabel(None)
        ax.set_title(caption)

    # Create the species color legend
    legend_elements = [
        Patch(facecolor=color, label=label)  # edgecolor="black"
        for label, color in species_key[::-1]  # flip to match Figure 3
    ]
    fig.legend(
        handles=legend_elements,  # The line objects
        # borderaxespad=0.1,    # Small spacing around legend box
        fontsize=5.5,
        fancybox=False,
        loc="upper center",
        edgecolor="white",  # i.e. no edge
        # mode="expand",
        ncol=len(species_key),  # i.e. one row
        columnspacing=1.1,
    )

    # Ensure can read the captions on the left:
    plt.subplots_adjust(left=0.2)

    # Display
    if img_filename:
        plt.savefig(img_filename, dpi=300, orientation="landscape")
        print(f"Drew {img_filename} from {tsv_filename}")
    else:
        # Interactive
        print(f"Displaying image from {tsv_filename}")
        plt.show()


plot(options.input, options.output)
