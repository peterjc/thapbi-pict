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


def connected_components(boolean_edge_matrix):
    """Given an N x N boolean symetric edge matrix, return a partition of the N nodes.

    For example, this 3 x 3 matrix breaks down into two connected components,
    node 0 and 2, and node 1 on its own:

    >>> import numpy as np
    >>> connected_components(np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], np.bool))
    [(0, 2), (1, )]

    The partition is sorted by size (largest component first)
    """
    # Could have built a NetworkX object, and used their code...
    n = boolean_edge_matrix.shape[0]
    if boolean_edge_matrix.shape != (n, n):
        raise ValueError("Expected a square adjacency matrix")
    partition = [{_} for _ in range(n)]
    for i in range(n):
        # So, what is node i connected to?
        for j in range(n):
            if boolean_edge_matrix[i, j] != boolean_edge_matrix[j, i]:
                raise ValueError("Expected a symmetric adjacency matrix")
            if i < j:
                continue
            if boolean_edge_matrix[i, j]:
                # Need to merge the partitions containing i and j
                i_old = [_ for _ in partition if i in _]
                assert len(i_old) == 1
                i_old = i_old[0]
                j_old = [_ for _ in partition if j in _]
                assert len(j_old) == 1, (
                    "%i found in %i members of the partitioning %r"
                    % (j, len(j_old), partition)
                )
                j_old = j_old[0]
                if i_old == j_old:
                    # Already in same partition
                    continue
                # Need to merge the partitions containing i and j
                partition = [_ for _ in partition if i not in _ and j not in _]
                partition.append(i_old.union(j_old))
    return [tuple(sorted(_)) for _ in sorted(partition, key=len, reverse=True)]


assert connected_components(np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], np.bool)) == [
    (0, 2),
    (1,),
]
assert connected_components(
    np.array([[0, 0, 1, 1], [0, 0, 0, 1], [1, 0, 0, 1], [1, 1, 1, 0]], np.bool)
) == [(0, 1, 2, 3)]
assert connected_components(
    np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 1], [0, 1, 1, 0]], np.bool)
) == [(0, 1, 2, 3)]
assert connected_components(
    np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]], np.bool)
) == [(0, 2), (1, 3)]


def write_pdf(graph, filename):
    """Render NetworkX graph to PDF using GraphViz fdp."""
    # TODO: Try "sfdp" but need GraphViz built with triangulation library
    default = graph["node_default"]["color"]
    node_colors = [graph.node[_].get("color", default) for _ in graph]

    default = 1.0
    node_sizes = [graph.node[_].get("size", default) for _ in graph]

    default = ""
    node_labels = [graph.node[_].get("label", default) for _ in graph]

    default = graph["edge_default"]["color"]
    edge_colors = [graph.edge[_].get("color", default) for _ in graph.edges]

    placement = nx.drawing.nx_pydot.graphviz_layout(graph, "fdp")
    nx.draw_networkx_nodes(
        graph, placement, node_color=node_colors, node_size=node_sizes
    )
    nx.draw_networkx_edges(
        graph,
        placement,
        # style=edge_styles,
        # width=edge_widths,
        edge_color=edge_colors,
        alpha=0.5,
    )
    nx.draw_networkx_labels(graph, placement, node_labels, font_size=4)
    plt.axis("off")
    plt.savefig(filename)


def main(
    graph_output,
    graph_format,
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
    if inputs is None:
        inputs = []
    assert isinstance(inputs, list)

    MIN_CLUMP = 150  # command line option?

    if 3 < max_edit_dist:
        sys.exit("ERROR: Maximum supported edit distance is 3bp.")

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

    if not md5_to_seq:
        sys.exit("ERROR: No sequences to plot.")

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

    clumps = connected_components(distances <= max_edit_dist)
    sys.stderr.write(
        "Max edit distance %i gave %i connected components (size %i to %i).\n"
        % (max_edit_dist, len(clumps), len(clumps[0]), len(clumps[-1]))
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

    md5_list = list(md5_to_seq)
    dropped = set()
    for clump in clumps:
        if MIN_CLUMP <= len(clump):
            # Keep them all
            continue
        for index in clump:
            md5 = md5_list[index]
            if total_min_abundance <= md5_abundance.get(md5, 0):
                # High abundance, include it (not applied to DB sequences)
                continue
            # # Ignore self-vs-self which will be zero distance
            # if i > 1 and min(distances[i, 0 : i - 1]) <= max_edit_dist:
            #    # Good, will draw this
            #    continue
            # elif i + 1 < n and min(distances[i, i + 1 :]) <= max_edit_dist:
            #    # Good, will draw this
            #    continue
            # No reason to draw this:
            dropped.add(md5)
    sys.stderr.write(
        "Dropped %i sequences with no siblings within maximum edit distance %i.\n"
        % (len(dropped), max_edit_dist)
    )

    if md5_abundance:
        SIZE = 100 / (
            max(md5_abundance.values()) - total_min_abundance
        )  # scaling factor
    else:
        # Happens with DB only graph,
        SIZE = 100
    graph = nx.Graph()
    graph.graph["node_default"] = {"color": "#8B0000", "size": 1.0}
    graph.graph["edge_default"] = {"color": "#808080", "weight": 1.0}
    for clump in clumps:
        if len(clump) < MIN_CLUMP:
            continue
        for index in clump:
            check1 = md5_list[index]
            if check1 in dropped:
                continue
            sp = md5_species.get(check1, [])
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
            # DB only entries get size one, FASTA entries can be up to 100.
            node_size = max(
                1, SIZE * (md5_abundance.get(check1, 0) - total_min_abundance)
            )
            graph.add_node(check1, color=node_color, size=node_size, label=node_label)

    edge_count = 0
    edge_count1 = edge_count2 = edge_count3 = 0
    edge_style = []
    edge_width = []
    edge_color = []
    redundant = 0
    for i, check1 in enumerate(md5_to_seq):
        if check1 in dropped:
            continue
        seq1 = md5_to_seq[check1]
        for j, check2 in enumerate(md5_to_seq):
            if i < j and check2 not in dropped:
                seq2 = md5_to_seq[check2]
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
                    graph.add_edge(
                        check1,
                        check2,
                        len=edge_length,
                        K=edge_length,
                        weight=edge_weight,
                    )
                    # edge_style.append("invis")
                    # ValueError: Unrecognized linestyle: invis
                    edge_style.append("dotted")
                    edge_width.append(0.1)
                    edge_color.append("#0000FF9F")  # blue for debug
                    # edge_color.append("#000000FF")  # fully transparent
                else:
                    # Some graph layout algorithms can use weight attr; some want int
                    # Larger weight makes it closer to the requested length.
                    # fdp default length is 0.3, neato is 1.0
                    graph.add_edge(
                        check1,
                        check2,
                        len=edge_length,
                        K=edge_length,
                        weight=edge_weight,
                    )
                    edge_count += 1
                    if dist <= 1:
                        edge_count1 += 1
                        edge_style.append("solid")
                        edge_width.append(1.0)
                        edge_color.append("#404040")
                    elif dist <= 2:
                        edge_count2 += 1
                        edge_style.append("dashed")
                        edge_width.append(0.33)
                        edge_color.append("#707070")
                    else:
                        edge_count3 += 1
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
        "pdf": write_pdf,
    }
    try:
        write_fn = render[graph_format]
    except KeyError:
        # Typically this would be caught in __main__.py
        sys.exit("ERROR: Unexpected graph output format: %s" % graph_format)

    write_fn(graph, "/dev/stdout" if graph_output == "-" else graph_output)

    return 0
