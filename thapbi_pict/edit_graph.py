# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

"""Generate edit-distance network graph from FASTA files.

This implementes the ``thapbi_pict edit-graph ...`` command.
"""

import os
import sys

from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from Levenshtein import distance as levenshtein

from sqlalchemy.orm import aliased, contains_eager

from Bio.SeqIO.FastaIO import SimpleFastaParser

from .db_orm import ITS1, SequenceSource, Taxonomy, connect_to_db

from .utils import find_requested_files
from .utils import genus_species_name
from .utils import md5seq
from .utils import species_level
from .utils import split_read_name_abundance


genus_color = {
    # From the VGA colors, in order of DB abundance,
    # dark red for other, grey for none, dark orange for conflicts
    "Phytophthora": "#FF0000",  # Red
    "Peronospora": "#00FF00",  # Lime
    "Hyaloperonospora": "#0000FF",  # Blue
    "Bremia": "#FFFF00",  # Yellow
    "Pseudoperonospora": "#00FFFF",  # Cyan
    "Plasmopara": "#FF00FF",  # Magenta
    "Nothophytophthora": "#800000",  # Maroon
    "Peronosclerospora": "#808000",  # Olive
    "Perofascia": "#008000",  # Green
    "Paraperonospora": "#800080",  # Purple
    "Protobremia": "#008080",  # Teal
    # Basidiophora
    # Calycofera
    # Plasmoverna
    "synthetic": "#FFA500",  # Orange
}


def write_pdf(G, filename):
    """Render NetworkX graph to PDF using GraphViz fdp."""
    # TODO: Try "sfdp" but need GraphViz built with triangulation library
    default = G.graph["node_default"]["color"]
    node_colors = [G.node[_].get("color", default) for _ in G]

    default = 1.0
    node_sizes = [G.node[_].get("size", default) for _ in G]

    default = ""
    node_labels = {_: G.node[_].get("label", default) for _ in G}

    default = G.graph["edge_default"]["style"]
    edge_styles = [G.edges[_].get("style", default) for _ in G.edges()]

    default = G.graph["edge_default"]["color"]
    edge_colors = [G.edges[_].get("color", default) for _ in G.edges()]

    default = G.graph["edge_default"]["width"]
    edge_widths = [G.edges[_].get("width", default) for _ in G.edges()]

    placement = nx.drawing.nx_pydot.graphviz_layout(G, "fdp")
    nx.draw_networkx_nodes(G, placement, node_color=node_colors, node_size=node_sizes)
    nx.draw_networkx_edges(
        G,
        placement,
        style=edge_styles,
        width=edge_widths,
        edge_color=edge_colors,
        alpha=0.5,
    )
    nx.draw_networkx_labels(G, placement, node_labels, font_size=4)
    plt.axis("off")
    if filename in ["-", "/dev/stdout"]:
        plt.savefig(sys.stdout.buffer, format="pdf")
    else:
        plt.savefig(filename, format="pdf")


def write_xgmml(G, filename, name="THAPBI PICT edit-graph"):
    """Save graph in XGMML format suitable for Cytoscape import."""
    # Not currently supported in NetworkX, and third party
    # package networkxgmml is not up to date (Python 3,
    # setting graphical properties on edges). So, DIY time!
    with open(filename, "w") as handle:
        handle.write('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n')
        handle.write(
            '<graph directed="0"  xmlns:dc="http://purl.org/dc/elements/1.1/" '
            'xmlns:xlink="http://www.w3.org/1999/xlink" '
            'xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" '
            'xmlns="http://www.cs.rpi.edu/XGMML">\n'
        )
        for n in G:
            node = G.nodes[n]
            try:
                label = node["label"].replace("\n", ";")  # Undo newline for Graphviz
            except KeyError:
                label = None
            if not label:
                label = "{%s}" % n[:6]  # start of MD5, prefixed for sorting
            # weight = node["weight"]
            handle.write('  <node id="%s" label="%s">\n' % (n, label))
            color = node["color"]
            # Size 1 to 100 works fine in PDF output, not so good in Cytoscape!
            # Rescale to use range 5 to 50.
            size = (node["size"] * 0.45) + 5.0
            handle.write(
                '    <graphics type="CIRCLE" fill="%s" outline="#000000" '
                'h="%0.2f" w="%0.2f"/>\n' % (color, size, size)
            )
            # Cytoscape hides the node ID (presumably assumes not usually user facing):
            handle.write('    <att type="string" name="MD5" value="%s"/>\n' % n)
            handle.write(
                '    <att type="integer" name="Total-abundance" value="%i"/>\n'
                % node.get("abundance", 0)
            )
            handle.write(
                '    <att type="integer" name="Sample-count" value="%i"/>\n'
                % node.get("sample_count", 0)
            )
            if node["genus"]:
                handle.write(
                    '    <att type="string" name="Genus" value="%s"/>\n' % node["genus"]
                )
            if node["taxonomy"]:
                handle.write(
                    '    <att type="string" name="Taxonomy" value="%s"/>\n'
                    % node["taxonomy"]
                )
            handle.write("  </node>\n")
        for n1, n2 in G.edges():
            edge = G.edges[n1, n2]
            handle.write(
                '  <edge source="%s" target="%s" weight="%0.2f">\n'
                % (n1, n2, edge["weight"])
            )
            # edge["style"]  # Not in XGMML?
            handle.write(
                '    <graphics fill="%s" width="%0.1f"/>\n'
                % (edge["color"], edge["width"])
            )
            if "edit_dist" in edge:
                handle.write(
                    '    <att type="integer" name="Edit-distance" value="%i"/>\n'
                    % edge["edit_dist"]
                )
            # Seems Cytoscape does not expose the weight attribute above
            # in the user-facing edge table, so doing it again here:
            handle.write(
                '    <att type="integer" name="Edit-distance-weight" value="%i"/>\n'
                % edge["weight"]
            )
            handle.write("  </edge>\n")
        handle.write("</graph>\n")


def main(
    graph_output,
    graph_format,
    db_url,
    inputs,
    min_abundance=100,
    always_show_db=False,
    total_min_abundance=0,
    max_edit_dist=3,
    debug=False,
):
    """Run the edit-graph command with arguments from the command line.

    Plan is to show sequences from a database (possibly with species/genus
    limits) and/or selected FASTA files (possibly with predictions or other
    metadata, and minimum abundance limits).
    """
    if inputs is None:
        inputs = []
    assert isinstance(inputs, list)

    if 3 < max_edit_dist:
        sys.exit("ERROR: Maximum supported edit distance is 3bp.")

    samples = set()
    md5_abundance = Counter()
    md5_sample_count = Counter()
    abundance_by_samples = {}
    md5_to_seq = {}
    md5_species = {}
    md5_in_db = set()
    md5_in_fasta = set()

    if not (inputs or db_url):
        sys.exit("Require -d / --database and/or -i / --input argument.")

    if not inputs and not always_show_db:
        sys.exit("If not using -i / --input argument, require -s / --showdb.")

    if inputs:
        if debug:
            sys.stderr.write("DEBUG: Loading FASTA sequences and abundances\n")
        for fasta_file in find_requested_files(inputs, ".fasta", debug):
            # TODO: Refactor this shared code with sample-summary?
            sample = os.path.basename(fasta_file).rsplit(".", 1)[0]
            if sample.startswith("Undetermined"):
                sys.stderr.write("WARNING: Ignoring %s\n" % fasta_file)
                continue
            if sample in samples:
                sys.exit("Duplicate sample name %s" % sample)
            samples.add(sample)
            with open(fasta_file) as handle:
                md5_warn = False
                for title, seq in SimpleFastaParser(handle):
                    idn, abundance = split_read_name_abundance(title.split(None, 1)[0])
                    md5 = md5seq(seq)
                    if idn != md5:
                        md5_warn = True
                    if min_abundance > 1 and abundance < min_abundance:
                        continue
                    md5_in_fasta.add(md5)
                    abundance_by_samples[md5, sample] = abundance
                    md5_abundance[md5] += abundance
                    md5_sample_count[md5] += 1
                    if md5 in md5_to_seq:
                        assert md5_to_seq[md5] == seq
                    else:
                        md5_to_seq[md5] = seq
                if md5_warn:
                    sys.stderr.write(
                        "WARNING: Sequence(s) in %s not using MD5_abundance naming\n"
                        % fasta_file
                    )
        sys.stderr.write(
            "Loaded %i unique sequences from %i FASTA files.\n"
            % (len(md5_in_fasta), len(samples))
        )
        # Drop low total abundance FASTA sequences now (before compute distances)
        if total_min_abundance:
            for md5, total in md5_abundance.items():
                if total < total_min_abundance:
                    # Remove it!
                    md5_in_fasta.remove(md5)
                    del md5_to_seq[md5]
            sys.stderr.write(
                "Minimum total abundance threshold %i "
                "left %i sequences from FASTA files.\n"
                % (total_min_abundance, len(md5_in_fasta))
            )

    if db_url:
        if debug:
            sys.stderr.write("DEBUG: Connecting to database %s\n" % db_url)
        # Connect to the DB,
        Session = connect_to_db(db_url, echo=False)  # echo=debug
        session = Session()

        # Doing a join to pull in the ITS1 and Taxonomy tables too:
        cur_tax = aliased(Taxonomy)
        its1_seq = aliased(ITS1)
        view = (
            session.query(SequenceSource)
            .join(its1_seq, SequenceSource.its1)
            .join(cur_tax, SequenceSource.current_taxonomy)
            .options(contains_eager(SequenceSource.its1, alias=its1_seq))
            .options(contains_eager(SequenceSource.current_taxonomy, alias=cur_tax))
        )
        # Sorting for reproducibility
        view = view.order_by(SequenceSource.id)
        # TODO - Copy genus/species filtering behvaiour from dump command?

        for seq_source in view:
            md5 = seq_source.its1.md5
            if not always_show_db and md5 not in md5_in_fasta:
                # Low abundance or absenst from FASTA files, ignore it
                continue
            md5_to_seq[md5] = seq_source.its1.sequence
            genus_species = genus_species_name(
                seq_source.current_taxonomy.genus, seq_source.current_taxonomy.species
            )
            try:
                md5_species[md5].add(genus_species)
            except KeyError:
                md5_species[md5] = {genus_species}

        for md5 in md5_species:
            for genus_species in [_ for _ in md5_species[md5] if species_level(_)]:
                genus = genus_species.split(None, 1)[0]
                if genus in md5_species[md5]:
                    # If have species level, discard genus level only
                    md5_species[md5].remove(genus)
        md5_in_db = set(md5_species)
        if always_show_db:
            sys.stderr.write(
                "Loaded %i unique sequences from database\n" % len(md5_in_db)
            )
        else:
            sys.stderr.write(
                "Matched %i unique sequences in database\n" % len(md5_in_db)
            )

    if db_url and inputs and always_show_db:
        sys.stderr.write(
            "DB had %i sequences (%i not in FASTA), "
            "FASTA had %i sequences (%i not in DB).\n"
            % (
                len(md5_in_db),
                len(md5_in_db.difference(md5_in_fasta)),
                len(md5_in_fasta),
                len(md5_in_fasta.difference(md5_in_db)),
            )
        )
        sys.stderr.write(
            "DB and FASTA had %i sequences in common; %i combined.\n"
            % (
                len(md5_in_db.intersection(md5_in_fasta)),
                len(md5_in_db.union(md5_in_fasta)),
            )
        )

    if not md5_to_seq:
        sys.exit("ERROR: No sequences to plot.")

    # For drawing performance reasons, calculate the distances, and then may
    # drop nodes with no edges (unless for example DB entry at species level,
    # or for environmental sequences at high abundance)
    md5_list = list(md5_to_seq)
    wanted = set()
    n = len(md5_list)
    distances = np.zeros((n, n), np.uint)
    for i, check1 in enumerate(md5_list):
        seq1 = md5_to_seq[check1]
        for j, check2 in enumerate(md5_list):
            if i < j:
                seq2 = md5_to_seq[check2]
                distances[i, j] = distances[j, i] = d = levenshtein(seq1, seq2)
                if d and d <= max_edit_dist:
                    wanted.add(check1)
                    wanted.add(check2)
    sys.stderr.write(
        "Computed %i Levenshtein edit distances between %i sequences.\n"
        % (n * (n - 1), n)
    )
    sys.stderr.write(
        "Will draw %i nodes with at least one edge (%i are isolated sequences).\n"
        % (len(wanted), n - len(wanted))
    )

    # Matrix computation of multi-step paths vs edit distances, e.g.
    # will use fact A-B is 1bp and B-C is 2bp to skip drawing A-C of 3bp.
    one_bp = distances == 1  # boolean
    two_step = np.dot(one_bp, one_bp)  # matrix multiply
    two_bp = (distances == 2) | two_step
    three_step = (
        np.dot(one_bp, two_bp) | np.dot(two_bp, one_bp) | np.dot(one_bp, one_bp, one_bp)
    )
    del one_bp, two_bp

    for md5 in md5_list:
        if md5 not in wanted:
            # Will include high abundance singletons too
            if total_min_abundance <= md5_abundance.get(md5, 0):
                wanted.add(md5)
    if inputs:
        sys.stderr.write(
            "Including high abundance isolated sequences, will draw %i nodes.\n"
            % len(wanted)
        )

    if md5_sample_count:
        # scaling factor
        SIZE = 100.0 / max(md5_sample_count.values())
    else:
        # Happens with DB only graph,
        SIZE = 1.0
    G = nx.Graph()
    G.graph["node_default"] = {"color": "#8B0000", "size": 1.0}
    G.graph["edge_default"] = {
        "color": "#FF0000",
        "weight": 1.0,
        "width": 1.0,
        "style": "solid",
    }
    for md5 in md5_list:
        if md5 not in wanted:
            continue
        sp = md5_species.get(md5, [])
        genus = sorted({_.split(None, 1)[0] for _ in sp})
        if not genus:
            node_color = "#808080"  # grey
        elif len(genus) > 1:
            node_color = "#FF8C00"  # dark orange
        elif genus[0] in genus_color:
            node_color = genus_color[genus[0]]
        else:
            node_color = "#8B0000"  # dark red
        if any(species_level(_) for _ in sp):
            node_label = "\n".join(sorted(sp)).replace("Phytophthora", "P.")
        else:
            # Genus only
            node_label = ""
        genus = ";".join(sorted({_.split(None, 1)[0] for _ in sp}))
        # DB only entries get size one, FASTA entries can be up to 100.
        abundance = md5_abundance.get(md5, 0)
        node_size = max(1, SIZE * md5_sample_count.get(md5, 0))
        G.add_node(
            md5,
            color=node_color,
            size=node_size,
            label=node_label,
            abundance=abundance,
            sample_count=md5_sample_count.get(md5, 0),
            genus=genus,
            taxonomy=";".join(sorted(sp)),
        )

    edge_count = 0
    edge_count1 = edge_count2 = edge_count3 = 0
    edge_style = []
    edge_width = []
    edge_color = []
    redundant = 0
    for i, check1 in enumerate(md5_list):
        if check1 not in wanted:
            continue
        # seq1 = md5_to_seq[check1]
        for j, check2 in enumerate(md5_list):
            if i < j and check2 in wanted:
                # seq2 = md5_to_seq[check2]
                # dist = levenshtein(seq1, seq2)
                dist = distances[i, j]
                # Some graph layout algorithms can use weight attr; some want int
                # Larger weight makes it closer to the requested length.
                # fdp default length is 0.3, neato is 1.0

                # i.e. edit distance 1, 2, 3 becomes distance 0.1, 0.2 and 0.3
                edge_length = 0.3 * dist / max_edit_dist
                # i.e. edit distance 1, 2, 3 get weights 3, 2, 1
                edge_weight = max_edit_dist - dist + 1

                if dist > max_edit_dist:
                    continue
                if (dist == 2 and two_step[i, j]) or (dist == 3 and three_step[i, j]):
                    # Redundant edge, if dist=2, two 1bp edges exist
                    # Or, if dist=3, three 1bp edges exist, or 1bp+2bp
                    redundant += 1
                    # edge_style = "invis"
                    # ValueError: Unrecognized linestyle: invis
                    edge_style = "dotted"
                    edge_width = 0.1
                    edge_color = "#0000FF9F"  # transparent blue for debug
                    # edge_color = "#000000FF"  # fully transparent
                    continue
                else:
                    # Some graph layout algorithms can use weight attr; some want int
                    # Larger weight makes it closer to the requested length.
                    # fdp default length is 0.3, neato is 1.0
                    edge_count += 1
                    if dist <= 1:
                        edge_count1 += 1
                        edge_style = "solid"
                        edge_width = 1.0
                        edge_color = "#404040"
                    elif dist <= 2:
                        edge_count2 += 1
                        edge_style = "dashed"
                        edge_width = 0.33
                        edge_color = "#707070"
                    else:
                        edge_count3 += 1
                        edge_style = "dotted"
                        edge_width = 0.25
                        edge_color = "#808080"
                G.add_edge(
                    check1,
                    check2,
                    len=edge_length,
                    edit_dist=dist,
                    K=edge_length,
                    weight=edge_weight,
                    style=edge_style,
                    width=edge_width,
                    color=edge_color,
                )

    if debug:
        sys.stderr.write(
            "DEBUG: %i edges up to maximum edit distance %i\n"
            % (edge_count, max_edit_dist)
        )
        sys.stderr.write(
            "DEBUG: %i one-bp edges; %i two-bp edges; %i three-bp edges.\n"
            % (edge_count1, edge_count2, edge_count3)
        )
        assert edge_count == edge_count1 + edge_count2 + edge_count3
        sys.stderr.write(
            "DEBUG: Dropped %i redundant 2-bp or 3-bp edges.\n" % redundant
        )

    render = {
        "gexf": nx.write_gexf,
        "gml": nx.write_gml,
        "graphml": nx.write_graphml,
        "xgmml": write_xgmml,
        "pdf": write_pdf,
    }
    try:
        write_fn = render[graph_format]
    except KeyError:
        # Typically this would be caught in __main__.py
        sys.exit("ERROR: Unexpected graph output format: %s" % graph_format)

    write_fn(G, "/dev/stdout" if graph_output == "-" else graph_output)

    return 0
