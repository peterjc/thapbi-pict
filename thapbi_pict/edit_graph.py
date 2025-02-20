# Copyright 2019-2024 by Peter Cock, The James Hutton Institute.
# Revisions copyright 2024 by Peter Cock, University of Strathclyde.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Generate edit-distance network graph from FASTA files.

This implements the ``thapbi_pict edit-graph ...`` command.
"""

from __future__ import annotations

import sys
from collections import Counter

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from rapidfuzz.distance import Levenshtein
from rapidfuzz.process import cdist
from sqlalchemy.orm import aliased
from sqlalchemy.orm import contains_eager

from .db_orm import connect_to_db
from .db_orm import MarkerDef
from .db_orm import MarkerSeq
from .db_orm import SeqSource
from .db_orm import Taxonomy
from .utils import genus_species_name
from .utils import md5seq
from .utils import parse_sample_tsv
from .utils import species_level
from .versions import check_rapidfuzz
from .versions import check_tools

genus_color = {
    # From the VGA colors, in order of DB abundance,
    # dark red for other, grey for none, dark orange for conflicts
    #
    # Oomycetes; Peronosporales:
    # ==========================
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
    #
    # Nematoda:
    # =========
    "Globodera": "#FF0000",  # Red
    "Heterodera": "#00FF00",  # Lime
    #
    # Special cases:
    # ==============
    "synthetic": "#FFA500",  # Orange
}


def write_pdf(G, handle) -> None:
    """Render NetworkX graph to PDF using GraphViz fdp."""
    # TODO: Try "sfdp" but need GraphViz built with triangulation library
    check_tools(["fdp"], debug=False)

    default = G.graph["node_default"]["color"]
    node_colors = [G.nodes[_].get("color", default) for _ in G]

    default = 1.0
    node_sizes = [G.nodes[_].get("size", default) for _ in G]

    default = ""
    node_labels = {_: G.nodes[_].get("label", default) for _ in G}

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
    plt.savefig(handle, format="pdf")


def write_xgmml(G, handle, name: str = "THAPBI PICT edit-graph") -> None:
    """Save graph in XGMML format suitable for Cytoscape import."""
    # Not currently supported in NetworkX, and third party
    # package networkxgmml is not up to date (Python 3,
    # setting graphical properties on edges). So, DIY time!
    handle.write(b'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n')
    handle.write(
        b'<graph directed="0"  xmlns:dc="http://purl.org/dc/elements/1.1/" '
        b'xmlns:xlink="http://www.w3.org/1999/xlink" '
        b'xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" '
        b'xmlns="http://www.cs.rpi.edu/XGMML">\n'
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
        handle.write(
            b'  <node id="%b" label="%b">\n'
            % (n.encode("ascii"), label.encode("ascii"))
        )
        color = node["color"]
        # Size 1 to 100 works fine in PDF output, not so good in Cytoscape!
        # Rescale to use range 5 to 50.
        size = (node["size"] * 0.45) + 5.0
        handle.write(
            b'    <graphics type="CIRCLE" fill="%b" outline="#000000" '
            b'h="%0.2f" w="%0.2f"/>\n' % (color.encode("ascii"), size, size)
        )
        # Cytoscape hides the node ID (presumably assumes not usually user facing):
        handle.write(
            b'    <att type="string" name="MD5" value="%b"/>\n' % n.encode("ascii")
        )
        handle.write(
            b'    <att type="integer" name="Total-abundance" value="%i"/>\n'
            % node.get("total_abundance", 0)
        )
        handle.write(
            b'    <att type="integer" name="Max-sample-abundance" value="%i"/>\n'
            % node.get("max_sample_abundance", 0)
        )
        handle.write(
            b'    <att type="integer" name="Sample-count" value="%i"/>\n'
            % node.get("sample_count", 0)
        )
        handle.write(
            b'    <att type="integer" name="in-db" value="%i"/>\n'
            % (1 if node.get("in_db", False) else 0)
        )
        if node["genus"]:
            # Are there any non-ASCII genus names?
            handle.write(
                b'    <att type="string" name="Genus" value="%b"/>\n'
                % node["genus"].encode("ascii")
            )
        if node["taxonomy"]:
            # TODO - how to encode any non-ASCII species names?
            # e.g. For hybrids the NCBI avoids the multiplication sign
            # "×" in favour of the letter "x".
            # https://doi.org/10.1093/database/baaa062
            handle.write(
                b'    <att type="string" name="Taxonomy" value="%b"/>\n'
                % node["taxonomy"].encode("ascii")
            )
        handle.write(b"  </node>\n")
    for n1, n2 in G.edges():
        edge = G.edges[n1, n2]
        handle.write(
            b'  <edge source="%b" target="%b" weight="%0.2f">\n'
            % (n1.encode("ascii"), n2.encode("ascii"), edge["weight"])
        )
        # edge["style"]  # Not in XGMML?
        handle.write(
            b'    <graphics fill="%b" width="%0.1f"/>\n'
            % (edge["color"].encode("ascii"), edge["width"])
        )
        if "edit_dist" in edge:
            handle.write(
                b'    <att type="integer" name="Edit-distance" value="%i"/>\n'
                % edge["edit_dist"]
            )
        # Seems Cytoscape does not expose the weight attribute above
        # in the user-facing edge table, so doing it again here:
        handle.write(
            b'    <att type="integer" name="Edit-distance-weight" value="%i"/>\n'
            % edge["weight"]
        )
        handle.write(b"  </edge>\n")
    handle.write(b"</graph>\n")


def main(
    graph_output: str,
    graph_format: str,
    db_url: str,
    input_file: str,
    min_abundance: int = 100,
    show_db_marker: str | None = None,
    total_min_abundance: int = 0,
    min_samples: int = 0,
    max_edit_dist: int = 3,
    ignore_prefixes: tuple[str, ...] | None = None,
    debug: bool = False,
) -> int:
    """Run the edit-graph command with arguments from the command line.

    This shows sequences from a database (possibly filtered with species/genus
    limits) and/or selected sample-tally TSV file (optionally with classifier
    output, and possibly with a minimum abundance limit set here).

    Computes a Levenshtein edit-distance matrix from the selected sequences,
    which can be exported as a matrix, but is usually converted into a graph
    of unique sequences as nodes, with short edit distances as edges.

    Graph node size is scaled by sample count (number of FASTA files that it
    appears in), and colored by assigned species (from a classifier TSV file).
    """
    assert isinstance(input_file, str) or input_file is None

    if 3 < max_edit_dist:
        sys.exit("ERROR: Maximum supported edit distance is 3bp.")

    samples = set()
    md5_abundance: dict[str, int] = Counter()
    md5_sample_count: dict[str, int] = Counter()
    abundance_by_samples = {}
    max_sample_abundance: dict[str, int] = {}
    md5_to_seq: dict[str, str] = {}
    md5_species = {}
    md5_in_db: set[str] = set()
    md5_in_fasta: set[str] = set()

    if not (input_file or db_url):
        sys.exit("Require -d / --database and/or -i / --input argument.")

    if not input_file and not show_db_marker:
        sys.exit(
            "If not using -i / --input argument, require -k / --marker to use DB only."
        )

    if input_file:
        if not input_file.endswith(".tsv"):
            sys.exit(f"ERROR: Expected a .tsv file, not {input_file}")
        else:
            filename = input_file
            if debug:
                sys.stderr.write(
                    f"DEBUG: Loading sequences sample tallies from {filename}\n"
                )
            # sample = file_to_sample_name(filename)
            # samples.add(sample)
            seqs, seq_meta, sample_headers, counts = parse_sample_tsv(
                filename, min_abundance=min_abundance, debug=debug
            )
            md5_warn = False
            for (marker, md5), seq in seqs.items():
                if show_db_marker and marker != show_db_marker:
                    continue
                if md5 != md5seq(seq):
                    md5_warn = True
                for sample in sample_headers:
                    samples.add(sample)
                    abundance = counts.get((marker, md5, sample), 0)
                    if not abundance:
                        continue
                    assert min_abundance <= abundance, (marker, md5, sample)
                    md5_to_seq[md5] = seq
                    md5_in_fasta.add(md5)
                    abundance_by_samples[md5, sample] = abundance
                    max_sample_abundance[md5] = max(
                        abundance, max_sample_abundance.get(md5, 0)
                    )
                    md5_abundance[md5] += abundance
                    md5_sample_count[md5] += 1
                    if md5 in md5_to_seq:
                        assert md5_to_seq[md5] == seq, f"{md5} vs {seq}"
                    else:
                        md5_to_seq[md5] = seq
            del seqs
            if md5_warn:
                sys.stderr.write(
                    f"WARNING: Sequence(s) in {filename}"
                    " not using MD5_abundance naming\n"
                )
            for (marker, idn), meta in seq_meta.items():
                if show_db_marker and marker != show_db_marker:
                    continue
                sp = meta.get("genus-species", "")
                if sp:
                    md5_species[idn] = set(sp.split(";"))
                del sp
            if debug:
                sys.stderr.write(
                    f"Have species assignments for {len(md5_species)} unique seqs\n"
                )
        sys.stderr.write(
            f"Loaded {len(md5_in_fasta)} unique sequences"
            f" from {len(samples)} samples.\n"
        )

        # Drop low total abundance FASTA sequences now (before compute distances)
        if total_min_abundance:
            for md5, total in md5_abundance.items():
                if total < total_min_abundance:
                    # Remove it!
                    md5_in_fasta.remove(md5)
                    del md5_to_seq[md5]
            sys.stderr.write(
                f"Minimum total abundance threshold {total_min_abundance}"
                f" left {len(md5_in_fasta)} sequences from input files.\n"
            )
        if min_samples:
            for md5, total in md5_sample_count.items():
                if total < min_samples:
                    # Remove it!
                    md5_in_fasta.remove(md5)
                    del md5_to_seq[md5]
            sys.stderr.write(
                f"Minimum sample threshold {min_samples}"
                f" left {len(md5_in_fasta)} sequences from input files.\n"
            )
        if len(md5_to_seq) > 6000:
            sys.stderr.write(
                "WARNING: Over 6000 sequences to plot; aborting edit-graph\n"
            )
            # Special return value for use within pipeline
            return 2

    if db_url:
        if debug:
            sys.stderr.write(f"DEBUG: Connecting to database {db_url}\n")
        # Connect to the DB,
        session = connect_to_db(db_url, echo=False)  # echo=debug

        # Doing a join to pull in the marker and taxonomy tables too:
        cur_tax = aliased(Taxonomy)
        marker_seq = aliased(MarkerSeq)
        view = (
            session.query(SeqSource)
            .join(marker_seq, SeqSource.marker_seq)
            .join(cur_tax, SeqSource.taxonomy)
            .options(contains_eager(SeqSource.marker_seq, alias=marker_seq))
            .options(contains_eager(SeqSource.taxonomy, alias=cur_tax))
        )
        if show_db_marker:
            # TODO - Check this marker is actually in the DB?
            # Note if marker not specified, will use all the DB entries to
            # label nodes regardless of which marker they are for - probably
            # fine, but could have same sequence against multiple markers.
            view = view.join(MarkerDef, SeqSource.marker_definition).filter(
                MarkerDef.name == show_db_marker
            )
        # Sorting for reproducibility
        view = view.order_by(SeqSource.id)
        # TODO - Copy genus/species filtering behaviour from dump command?

        for seq_source in view:
            md5 = seq_source.marker_seq.md5
            if not show_db_marker and md5 not in md5_in_fasta:
                # Low abundance or absent from FASTA files, ignore it
                continue
            md5_in_db.add(md5)
            md5_to_seq[md5] = seq_source.marker_seq.sequence
            genus_species = genus_species_name(
                seq_source.taxonomy.genus, seq_source.taxonomy.species
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
        if show_db_marker:
            sys.stderr.write(
                f"Loaded {len(md5_in_db)} unique {show_db_marker} sequences from DB.\n"
            )
        else:
            sys.stderr.write(
                f"Matched {len(md5_in_db)} unique sequences in database.\n"
            )

    if db_url and input_file and show_db_marker:
        sys.stderr.write(
            f"DB had {len(md5_in_db)} sequences"
            f" ({len(md5_in_db.difference(md5_in_fasta))} not in FASTA),"
            f" FASTA had {len(md5_in_fasta)} sequences"
            f" ({len(md5_in_fasta.difference(md5_in_db))} not in DB).\n"
        )
        sys.stderr.write(
            f"DB and FASTA had {len(md5_in_db.intersection(md5_in_fasta))} sequences"
            f" in common; {len(md5_in_db.union(md5_in_fasta))} combined.\n"
        )

    if not md5_to_seq:
        sys.exit("ERROR: No sequences to plot.")

    # For drawing performance reasons, calculate the distances, and then may
    # drop nodes with no edges (unless for example DB entry at species level,
    # or for environmental sequences at high abundance)
    check_rapidfuzz()
    n = len(md5_to_seq)
    md5_list = sorted(md5_to_seq)
    seq_list: list[str] = [md5_to_seq[_] for _ in md5_list]
    # Will get values 0, 1, ..., max_edit_dist, or
    # max_edit_dist+1 if distance is higher (was -1 prior to rapidfuzz v2.0.0)
    distances = cdist(
        seq_list,
        seq_list,
        scorer=Levenshtein.distance,
        dtype=np.int16 if graph_format == "matrix" else np.int8,
        score_cutoff=None if graph_format == "matrix" else max_edit_dist,
    )
    del seq_list
    sys.stderr.write("Computed Levenshtein edit distances.\n")
    assert min(min(_) for _ in distances) == 0, (
        f"Possible overflow, min distance {min(min(_) for _ in distances)} not zero."
    )

    if graph_format == "matrix":
        # Report all nodes, even if isolated and low abundance
        # i.e. ignores the wanted list used for plotting
        if graph_output in ("-", "/dev/stdout"):
            handle = sys.stdout
        else:
            handle = open(graph_output, "w")
        cols = "\t".join(md5_list)
        handle.write(f"MD5\tSpecies\t{cols}\n")
        del cols
        for i, md5 in enumerate(md5_list):
            sp = ",".join(sorted(md5_species.get(md5, [])))
            dists = "\t".join(str(_) for _ in distances[i])
            handle.write(f"{md5}\t{sp}\t{dists}\n")
            del sp, dists
        if graph_output != "-":
            handle.close()
        return 0

    # Isolated node's distances will be a single 0 (self) and (n-1) of (max_edit_dist+1)
    wanted = {
        md5
        for i, md5 in enumerate(md5_list)
        if set(distances[i]) != {0, max_edit_dist + 1}
    }
    sys.stderr.write(
        f"Will draw {len(wanted)} nodes with at least one edge"
        f" ({n - len(wanted)} are isolated sequences).\n"
    )
    for md5 in md5_list:
        if md5 not in wanted:
            # Will include high abundance singletons too
            if total_min_abundance <= md5_abundance.get(md5, 0):
                wanted.add(md5)
    if input_file:
        sys.stderr.write(
            "Including high abundance isolated sequences,"
            f" will draw {len(wanted)} nodes.\n"
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
        sp_list: list[str] = sorted(set(md5_species.get(md5, [])))
        genera = {_.split(None, 1)[0] for _ in sp_list}
        if md5 not in md5_in_db or not genus:
            node_color = "#808080"  # grey
        elif len(genera) > 1:
            node_color = "#FF8C00"  # dark orange
        elif genus[0] in genus_color:
            node_color = genus_color[genus[0]]
        else:
            node_color = "#8B0000"  # dark red
        if sp_list:
            # TODO - Remove this Phytophthora specific hack, or automate it?
            node_label = "\n".join(sp_list).replace("Phytophthora", "P.")
        else:
            # No species, not even genus only - fall back on MD5 as node ID
            node_label = ""
        # DB only entries get size zero, FASTA entries can be up to 100.
        node_size = max(1, SIZE * md5_sample_count.get(md5, 0))
        G.add_node(
            md5,
            color=node_color,
            size=node_size,
            label=node_label,
            total_abundance=md5_abundance.get(md5, 0),
            max_sample_abundance=max_sample_abundance.get(md5, 0),
            sample_count=md5_sample_count.get(md5, 0),
            genus=";".join(sorted(genera)),
            taxonomy=";".join(sp_list),
            in_db=md5 in md5_in_db,
        )
        del sp_list

    edge_count = 0
    edge_count1 = edge_count2 = edge_count3 = 0
    edge_style = ""
    edge_width = 0.0
    edge_color = ""
    redundant = 0
    for i, check1 in enumerate(md5_list):
        if check1 not in wanted:
            continue
        # seq1 = md5_to_seq[check1]
        for j, check2 in enumerate(md5_list):
            if i < j and check2 in wanted:
                # seq2 = md5_to_seq[check2]
                # dist = levenshtein(seq1, seq2)
                dist = int(distances[i, j])  # casting to drop numpy dtype
                if dist > max_edit_dist:
                    continue

                # Some graph layout algorithms can use weight attr; some want int
                # Larger weight makes it closer to the requested length.
                # fdp default length is 0.3, neato is 1.0

                # i.e. edit distance 1, 2, 3 becomes distance 0.1, 0.2 and 0.3
                edge_length = 0.3 * dist / max_edit_dist
                # i.e. edit distance 1, 2, 3 get weights 3, 2, 1
                edge_weight = max_edit_dist - dist + 1

                if (
                    dist == 2
                    and np.logical_and(distances[i, :] == 1, distances[:, j] == 1).any()
                ):
                    # Redundant edge dist=2 with a path with two 1bp edges
                    redundant += 1
                    continue
                elif dist == 3 and (
                    np.logical_and(distances[i, :] == 1, distances[:, j] == 2).any()
                    or np.logical_and(distances[i, :] == 2, distances[:, j] == 1).any()
                ):
                    # Redundant edge dist=3 where there is a path of 1bp and 2bp edges
                    # (if there is a path of 1bp, 1bp, 1bp, then there are also two
                    # routes of 1bp, 2bp and 2bp, 1bp as well which we'll find)
                    redundant += 1
                    continue
                else:
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
                    edit_dist=int(dist),  # cast numpy.uint64 to int
                    K=edge_length,
                    weight=edge_weight,
                    style=edge_style,
                    width=edge_width,
                    color=edge_color,
                )

    if debug:
        sys.stderr.write(
            f"DEBUG: {edge_count} edges up to maximum edit distance {max_edit_dist}\n"
        )
        sys.stderr.write(
            f"DEBUG: {edge_count1} one-bp edges; {edge_count2} two-bp edges;"
            f" {edge_count3} three-bp edges.\n"
        )
        assert edge_count == edge_count1 + edge_count2 + edge_count3
        sys.stderr.write(f"DEBUG: Dropped {redundant} redundant 2-bp or 3-bp edges.\n")

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
        sys.exit(f"ERROR: Unexpected graph output format: {graph_format}")

    # Output is in bytes (following NetworkX functions):
    if graph_output in ("-", "/dev/stdout"):
        write_fn(G, sys.stdout.buffer)
    else:
        with open(graph_output, "wb") as output_handle:
            write_fn(G, output_handle)

    return 0
