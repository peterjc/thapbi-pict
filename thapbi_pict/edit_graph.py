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
from .utils import species_level
from .utils import split_read_name_abundance


genus_color = {
    # From the VGA colors, in order of DB abundance,
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
}


def main(
    graph_output,
    db_url,
    inputs,
    min_abundance=100,
    total_min_abundance=1000,
    max_edit_dist=3,
    debug=False,
):
    """Run the edit-graph command with arguments from the command line.

    Plan is to show sequences from a database (possibly with species/genus
    limits) and/or selected FASTA files (possibly with predictions or other
    metadata, and minimum abundance limits).
    """
    assert isinstance(inputs, list)

    samples = set()
    md5_abundance = Counter()
    abundance_by_samples = {}
    md5_to_seq = {}
    md5_species = {}
    md5_in_db = set()
    md5_in_fasta = set()

    if not (inputs or db_url):
        sys.exit("Require -d / --database and/or -i / --input argument.")

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
            # if debug and len(md5_to_seq) > 1000:
            #    break
            md5 = seq_source.its1.md5
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
        sys.stderr.write("Loaded %i unique sequences from database\n" % len(md5_in_db))

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
                for title, seq in SimpleFastaParser(handle):
                    md5, abundance = split_read_name_abundance(title.split(None, 1)[0])
                    if min_abundance > 1 and abundance < min_abundance:
                        continue
                    md5_in_fasta.add(md5)
                    abundance_by_samples[md5, sample] = abundance
                    md5_abundance[md5] += abundance
                    if md5 in md5_to_seq:
                        assert md5_to_seq[md5] == seq
                    else:
                        md5_to_seq[md5] = seq
        sys.stderr.write(
            "Loaded %i unique sequences from %i FASTA files.\n"
            % (len(md5_in_fasta), len(samples))
        )
        # Drop low total abundance FASTA sequences now (before compute distances)
        for md5, total in md5_abundance.items():
            if total < total_min_abundance and md5 not in md5_in_db:
                # Remove it!
                md5_in_fasta.remove(md5)
                del md5_to_seq[md5]
        sys.stderr.write(
            "Minimum total abundance threshold %i left %i sequences from FASTA files.\n"
            % (total_min_abundance, len(md5_in_fasta))
        )

    if db_url and inputs:
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

    # For drawing performance reasons, calculate the distances, and then may
    # drop nodes with no edges (unless for example DB entry at species level,
    # or for environmental sequences at high abundance)
    n = len(md5_to_seq)
    distances = np.zeros((n, n), np.uint)
    for i, check1 in enumerate(md5_to_seq):
        seq1 = md5_to_seq[check1]
        for j, check2 in enumerate(md5_to_seq):
            if i < j:
                seq2 = md5_to_seq[check2]
                distances[i, j] = distances[j, i] = levenshtein(seq1, seq2)
    sys.stderr.write(
        "Computed %i Levenshtein edit distances between %i sequences.\n"
        % (n * (n - 1), n)
    )
    dropped = set()
    for i, md5 in enumerate(md5_to_seq):
        if total_min_abundance <= md5_abundance.get(md5, 0):
            # High abundance, include it (not applied to DB sequences)
            continue
        # Ignore self-vs-self which will be zero distance
        if i > 1 and min(distances[i, 0 : i - 1]) <= max_edit_dist:
            # Good, will draw this
            continue
        elif i + 1 < n and min(distances[i, i + 1 :]) <= max_edit_dist:
            # Good, will draw this
            continue
        # No reason to draw this:
        dropped.add(md5)
    sys.stderr.write(
        "Dropped %i sequences with no siblings within maximum edit distance %i.\n"
        % (len(dropped), max_edit_dist)
    )

    # SIZE = 100 / (total_max_abundance - total_min_abundance)  # scaling factor
    graph = nx.Graph()
    node_colors = []
    node_labels = {}
    node_sizes = []
    for check1 in md5_to_seq:
        if check1 in dropped:
            continue
        graph.add_node(check1)
        sp = md5_species.get(check1, [])
        genus = sorted({_.split(None, 1)[0] for _ in sp})
        if not genus:
            node_colors.append("#808080")  # grey
        elif len(genus) > 1:
            node_colors.append("#FF8C00")  # dark orange
        elif genus[0] in genus_color:
            node_colors.append(genus_color[genus[0]])
        else:
            node_colors.append("#8B0000")  # dark red
        if any(species_level(_) for _ in sp):
            node_labels[check1] = "\n".join(sorted(sp)).replace("Phytophthora", "P.")
        else:
            # Genus only
            pass
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

    edge_count = 0
    edge_style = []
    edge_width = []
    edge_color = []
    for i, check1 in enumerate(md5_to_seq):
        if check1 in dropped:
            continue
        seq1 = md5_to_seq[check1]
        for j, check2 in enumerate(md5_to_seq):
            if i < j and check2 not in dropped:
                seq2 = md5_to_seq[check2]
                # dist = levenshtein(seq1, seq2)
                dist = distances[i, j]
                if dist <= max_edit_dist:
                    # Some graph layout algorithms can use weight attr; some want int
                    # Larger weight makes it closer to the requested length.
                    # fdp default length is 0.3, neato is 1.0
                    graph.add_edge(check1, check2, len=0.3 * dist / max_edit_dist)
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
