#!/usr/bin/env python3
"""Turn ``summary/*.samples.1s5g.tsv`` into ``figure3reproduction.tsv`` output.

The output file ``figure3reproduction.tsv`` (THAPBI PICT output) has the same
format as ``figure3.tsv`` (Batovska et al. 2021 data underlying Figure 3,
exported using our ``figure3.R`` script). Either TSV file can be plotted in
the style of Figure 3 using the ``recreate_figure3.py`` Python script.
"""

import numpy as np

species_key = [
    # Caption (from Figure 3),
    # species list (note merging two Acizzia),
    # hex colour (from Figure 3)
    (
        "Acizzia alternata/solanicola",
        ("Acizzia alternata", "Acizzia solanicola"),
        "#0c4687",  # dark blue
    ),
    (
        "Bactericera cockerelli (TPP)",
        ("Bactericera cockerelli",),
        "#ae0707",  # dark red
    ),
    (
        "Diuraphis noxia (RWA)",
        ("Diuraphis noxia",),
        "#fa6e24",  # dark orange
    ),
    (
        "Metopolophium dirhodum",
        ("Metopolophium dirhodum",),
        "#3a9e82",  # Green / Blue
    ),
    (
        "Rhopalosiphum padi",
        ("Rhopalosiphum padi",),
        "#95cf77",  # Green Yellow
    ),
]

# First entry appears at bottom of the plot!
row_captions = [
    f"{count} Pool {pool}"
    for pool in (1, 2, 3, 4, 5)
    for count in (100, 250, 500, 1000)
][::-1]

# Expected ratios extracted from here:
# https://github.com/alexpiper/HemipteraMetabarcodingMS/blob/master/sample_data/expected/exp_seqtab.csv
# https://github.com/alexpiper/HemipteraMetabarcodingMS/blob/master/sample_data/expected/exp_taxtab.csv
# Uses same species order as our variable species_key
# Community order matches the display in Figure 3
expected = {
    # Mock community: A. solanicola, B. cockerelli, D. noxia, M. dirhodum, R. padi
    "100-Pool-1": (99, 1, 0, 0, 0),
    "250-Pool-1": (249, 1, 0, 0, 0),
    "500-Pool-1": (499, 1, 0, 0, 0),
    "1000-Pool-1": (999, 1, 0, 0, 0),
    "100-Pool-2": (0, 0, 1, 5, 94),
    "250-Pool-2": (0, 0, 1, 14, 235),
    "500-Pool-2": (0, 0, 1, 25, 474),
    "1000-Pool-2": (0, 0, 1, 50, 949),
    "100-Pool-3": (45, 5, 1, 15, 34),
    "250-Pool-3": (119, 5, 1, 39, 86),
    "500-Pool-3": (236, 5, 1, 79, 179),
    "1000-Pool-3": (476, 5, 1, 159, 359),
    "100-Pool-4": (49, 1, 5, 30, 15),
    "250-Pool-4": (125, 1, 5, 79, 40),
    "500-Pool-4": (258, 1, 5, 157, 79),
    "1000-Pool-4": (518, 1, 5, 317, 159),
    "100-Pool-5": (50, 0, 0, 25, 25),
    "250-Pool-5": (125, 0, 0, 63, 62),
    "500-Pool-5": (250, 0, 0, 125, 125),
    "1000-Pool-5": (500, 0, 0, 250, 250),
}


def load_species_counts(tsv_filename):
    """Load the five classification counts from a per-marker summary file."""
    species_of_interest = [_[1] for _ in species_key]
    samples_of_interest = [_.replace(" ", "-") for _ in row_captions]
    answer = {}
    with open(tsv_filename) as handle:
        header = handle.readline()
        assert header.startswith("#")
        fields = header[1:].rstrip("\n").split("\t")
        sample_name_col = fields.index("sample_alias")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts[sample_name_col] in samples_of_interest:
                values = [0] * len(species_of_interest)
                for i, species_list in enumerate(species_of_interest):
                    for j, field in enumerate(fields):
                        # print(f"Does field {field} match {species_list}?")
                        if set(field.split(";")).intersection(species_list):
                            values[i] += int(parts[j])
                answer[parts[sample_name_col]] = values
                assert sum(values), line
    assert len(answer) == len(
        samples_of_interest
    ), f"Looking for {samples_of_interest}, got {answer}"
    return answer


def make_into_ratios(data_dict):
    """Convert dict of counts into a dict of ratios."""
    return {key: np.array(value) / sum(value) for key, value in data_dict.items()}


ratio_expt = make_into_ratios(expected)
ratio_COI = make_into_ratios(load_species_counts("summary/COI.samples.onebp.tsv"))
ratio_18S = make_into_ratios(load_species_counts("summary/18S.samples.onebp.tsv"))
ratio_12S = make_into_ratios(load_species_counts("summary/12S.samples.onebp.tsv"))

# Output one row per sample/genus,
# 1. Sample caption (from original Figure 3)
# 2. Genus (abbreviated from hard coded species names)
# 3. Expected ratio, or read count"
# 4. COI ratio, or read count
# 5. 18S ratio, or read count
# 6. 12S ratio, or read count

with open("figure3reproduction.tsv", "w") as handle:
    handle.write("Caption\tGenus\tExpected\tCOI\t18S\t12S\n")
    for sample in expected:
        for i, genus in enumerate([_[0].split()[0] for _ in species_key]):
            handle.write(
                f"{sample.replace('-', ' ')}\t{genus}"
                f"\t{ratio_expt[sample][i]}\t{ratio_COI[sample][i]}"
                f"\t{ratio_18S[sample][i]}\t{ratio_12S[sample][i]}\n"
            )
print("Wrote figure3reproduction.tsv")
