"""Generate edit-distance network graph from FASTA files.

This implementes the ``thapbi_pict edit-graph ...`` command.
"""

import sys

import matplotlib.pyplot as plt
import networkx as nx

from Levenshtein import distance as levenshtein

from sqlalchemy.orm import aliased, contains_eager

from .db_orm import ITS1, SequenceSource, Taxonomy, connect_to_db

from .utils import genus_species_name, species_level


def main(
    graph_output,
    db_url,
    # inputs,
    # graph_output,
    # method="identity",
    # min_abundance=100,
    # total_min_abundance=1000,
    max_edit_dist=3,
    debug=False,
):
    """Run the edit-graph command with arguments from the command line.

    Plan is to show sequences from a database (possibly with species/genus
    limits) and/or selected FASTA files (possibly with predictions or other
    metadata, and minimum abundance limits).
    """
    md5_to_seq = {}
    md5_species = {}

    # Connect to the DB,
    Session = connect_to_db(db_url, echo=debug)
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
        md5_to_seq[md5] = seq_source.its1.sequence
        genus_species = genus_species_name(
            seq_source.current_taxonomy.genus, seq_source.current_taxonomy.species
        )
        try:
            md5_species[md5].add(genus_species)
        except KeyError:
            md5_species[md5] = {genus_species}

    # SIZE = 100 / (total_max_abundance - total_min_abundance)  # scaling factor
    graph = nx.Graph()
    node_colors = []
    node_labels = {}
    node_sizes = []
    for check1 in md5_to_seq:
        graph.add_node(check1)
        sp = md5_species[check1]
        if not sp:
            node_colors.append("#808080")
        elif species_level(sp):
            node_colors.append("#ff0000")
            node_labels[check1] = "\n".join(sorted(sp)).replace("Phytophthora", "P.")
        else:
            # Genus only, dark red
            node_colors.append("#600000")
        # node_sizes.append(SIZE * (md5_abundance[check1] - total_min_abundance))
        node_sizes.append(1)
    if debug:
        sys.stderr.write(
            "Node sizes %0.2f to %0.2f\n" % (min(node_sizes), max(node_sizes))
        )

    # for check1 in md5_abundance:
    #    if check1 in node_labels:
    #        if debug:
    #            node_labels[check1] += "\n" + check1[:6]
    #    else:
    #        node_labels[check1] = check1[:6]

    if debug:
        sys.stderr.write(
            "DEBUG: About to compute edit distances between %i unique sequences\n"
            % len(md5_to_seq)
        )
    edge_count = 0
    edge_style = []
    edge_width = []
    edge_color = []
    for i, check1 in enumerate(md5_to_seq):
        seq1 = md5_to_seq[check1]
        for j, check2 in enumerate(md5_to_seq):
            if i < j:
                seq2 = md5_to_seq[check2]
                dist = levenshtein(seq1, seq2)
                if dist <= max_edit_dist:
                    # Some graph layout algorithms can use weight attr
                    graph.add_edge(check1, check2, weight=1.0 / dist)
                    edge_count += 1
                    if dist <= 1:
                        edge_style.append("solid")
                        edge_width.append(1.0)
                        edge_color.append("#404040")
                    elif dist <= 2:
                        edge_style.append("dashed")
                        edge_width.append(0.33)
                        edge_color.append("#707070")
                    else:
                        edge_style.append("dotted")
                        edge_width.append(0.25)
                        edge_color.append("#808080")
                    # if debug:
                    #    sys.stderr.write("%s\t%s\t%i\n" % (check1, check2, dist))

    if debug:
        sys.stderr.write(
            "DEBUG: %i edges up to maximum edit distance %i\n"
            % (edge_count, max_edit_dist)
        )

    placement = nx.drawing.nx_pydot.graphviz_layout(graph, "fdp")
    nx.draw_networkx_nodes(
        graph,
        placement,
        node_color=node_colors,
        # node_labels=node_labels,
        node_size=node_sizes,
    )
    nx.draw_networkx_edges(
        graph,
        placement,
        style=edge_style,
        width=edge_width,
        edge_color=edge_color,
        alpha=0.5,
    )
    nx.draw_networkx_labels(graph, placement, node_labels, font_size=4)
    plt.axis("off")
    plt.savefig(graph_output)

    return 0
