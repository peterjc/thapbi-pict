# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Code for THAPBI PICT to deal with NCBI taxonomy dumps.

The code is needed initially for loading an NCBI taxdump folder (files
``names.dmp``, ``nodes.dmp`` etc) into a marker database.
"""
import os
import sys

from .db_orm import connect_to_db
from .db_orm import Synonym
from .db_orm import Taxonomy


def load_nodes(nodes_dmp, wanted_ranks=None):
    """Load the NCBI taxdump nodes.dmp file.

    Returns two dicts, the parent/child relationships, and the ranks (values
    are lists of taxids).

    Default is all ranks, can provide a possibly empty list/set of ranks of
    interest.
    """
    tree = {}
    ranks = {}  # keys are rank names, values are lists of ids
    with open(nodes_dmp) as handle:
        for line in handle:
            parts = line.split("\t|\t", 3)
            taxid = int(parts[0].strip())
            parent = int(parts[1].strip())
            rank = parts[2].strip()
            tree[taxid] = parent
            if wanted_ranks is None or rank in wanted_ranks:
                try:
                    ranks[rank].append(taxid)
                except KeyError:
                    ranks[rank] = [taxid]
    return tree, ranks


def load_names(names_dmp, wanted=None):
    """Load scientific names of species from NCBI taxdump names.dmp file."""
    names = {}
    synonym = {}
    with open(names_dmp) as handle:
        for line in handle:
            parts = line.split("\t|\t", 3)
            taxid = int(parts[0].strip())
            if wanted and taxid not in wanted:
                continue
            name = parts[1].strip()
            if line.endswith("\t|\tscientific name\t|\n"):
                names[taxid] = name
            elif line.endswith("\t|\tsynonym\t|\n") or line.endswith(
                "\t|\tincludes\t|\n"
            ):
                # e.g. Phytophthora aquimorbida 'includes' Phytophthora sp. CCH-2009b
                # which we want to treat like a synonym
                try:
                    synonym[taxid].append(name)
                except KeyError:
                    synonym[taxid] = [name]
    return names, synonym


def filter_tree(tree, ranks, ancestors):
    """Return a filtered version of the tree & ranks dict.

    NOTE: Does NOT preserve the original dict order.
    """
    count = len(tree)
    wanted = {a: a for a in ancestors}  # new tree
    reject = set()
    while tree:
        taxid, parent = tree.popitem()
        if taxid in wanted:
            # sys.stderr.write(f"Already took {taxid}\n")
            continue
        if taxid in reject:
            # sys.stderr.write(f"Already rejected {taxid}\n")
            continue
        if parent in wanted:
            # sys.stderr.write(f"Taking {taxid} via {parent}\n")
            wanted[taxid] = parent
            continue
        if taxid == parent or parent in reject:
            # sys.stderr.write(f"Reject {taxid} via {parent}\n")
            reject.add(taxid)
            continue
        t = parent
        batch = [parent]
        while True:
            t = tree[t]
            if t in wanted:
                # sys.stderr.write(f"Taking {taxid} via {batch}\n")
                wanted[taxid] = parent  # Not in batch, already popped
                wanted.update((_, tree[_]) for _ in batch)
                break
            batch.append(t)
            if t in reject or t == tree[t]:
                # sys.stderr.write(f"Rejecting {taxid} via {batch}\n")
                reject.add(taxid)
                reject.update(batch)
                break
    assert count == len(wanted) + len(reject)

    return wanted, {
        rank: sorted(set(values).intersection(wanted)) for rank, values in ranks.items()
    }


def get_ancestor(taxid, tree, stop_nodes):
    """Walk up tree until reach a stop node, or root."""
    t = taxid
    while True:
        if t in stop_nodes:
            return t
        elif not t or t == tree[t]:
            return t  # root
        else:
            t = tree[t]


def top_level_species(tree, ranks, names):
    """Find taxids for species under the genus_list.

    Our "genus" list matches the NCBI rank "genus", and includes child nodes
    as aliases (unless they fall on our "species" list or reject list of
    "environmental samples" or "unclassified <genus>").

    However, our "species" list are either NCBI rank "species" or "species
    group", where the parent is "genus" or "subgenus".

    This is to exclude NCBI rank "species" elements under "no rank" nodes like
    "environmental samples" or "unclassified Phytophthora" (taxid 211524),
    which we want to treat as genus aliases.

    Yields (species taxid, genus taxid) tuples.
    """
    genus_list = set(ranks["genus"])
    species_or_spgroup = set(ranks["species"] + ranks["species group"])
    genus_or_subgenus = set(ranks["genus"] + ranks["subgenus"])
    for taxid in tree:
        if taxid in species_or_spgroup and tree[taxid] in genus_or_subgenus:
            # Want this as a "species"
            assert taxid not in genus_or_subgenus
            genus_taxid = get_ancestor(taxid, tree, genus_list)
            yield taxid, genus_taxid


def not_top_species(tree, ranks, names, synonyms, top_species):
    """Find all 'minor' species, takes set of species taxid to ignore.

    Intended usage is to map minor-species as genus aliases, for
    example all the species under unclassified Phytophthora will
    be treated as synonyms of the genus Phytophthora.

    Yields (genus taxid, alias name) tuples.
    """
    stop_nodes = set(top_species).union(ranks["genus"])
    for taxid in tree:
        if taxid in stop_nodes:
            continue
        parent = get_ancestor(taxid, tree, stop_nodes)
        if taxid != parent:
            aliases = {names[taxid]}
            if taxid in synonyms:
                aliases.update(synonyms[taxid])
            for name in sorted(aliases):
                yield parent, name


def reject_name(species):
    """Species name be rejected, and not recorded in the DB."""
    return (
        species.split(None, 1)[0] in ("unclassified", "uncultured", "unidentified")
        or species == "environmental samples"
    )


def main(tax, db_url, ancestors, debug=True):
    """Load an NCBI taxdump into a database."""
    if not os.path.isdir(tax):
        sys.exit(f"ERROR: Could not find taxdump directory: {tax!r}\n")
    for filename in ("names.dmp", "nodes.dmp"):
        if not os.path.isfile(os.path.join(tax, filename)):
            sys.exit(f"ERROR: Missing {filename} in the taxdump directory: {tax!r}\n")

    try:
        ancestors = [int(_) for _ in ancestors.split(",")]
    except ValueError:
        sys.exit(f"ERROR: Invalid ancestors argument: {ancestors!r}\n")

    tree, ranks = load_nodes(
        os.path.join(tax, "nodes.dmp"),
        ("genus", "subgenus", "species group", "species"),
    )
    if debug:
        sys.stderr.write(f"Loaded {len(tree)} nodes from nodes.dmp\n")
    tree, ranks = filter_tree(tree, ranks, ancestors)
    if debug:
        sys.stderr.write(f"Reduced to {len(tree)} nodes under ancestors\n")
    genus_list = sorted(ranks["genus"])
    if not genus_list:
        sys.exit("ERROR: Could not identify any genus names under the given nodes\n")
    tree, ranks = filter_tree(tree, ranks, genus_list)
    if debug:
        sys.stderr.write(
            f"Reduced to {len(tree)} nodes under {len(genus_list)} genera\n"
        )
    assert genus_list == sorted(ranks["genus"])

    names, synonyms = load_names(os.path.join(tax, "names.dmp"), tree)
    if debug:
        sys.stderr.write(f"Loaded {len(names)} scientific names from names.dmp\n")

    if debug:
        sys.stderr.write(
            f"Identified {len(genus_list)} genera under specified ancestor node:"
            f" {', '.join(sorted(names[_] for _ in genus_list))}\n"
        )

    genus_species = sorted(top_level_species(tree, ranks, names))
    if debug:
        sys.stderr.write(f"Filtered down to {len(genus_species)} species names\n")

    minor_species = sorted(
        not_top_species(tree, ranks, names, synonyms, {_[0] for _ in genus_species})
    )
    if debug:
        sys.stderr.write(
            f"Treating {len(minor_species)} minor species names as genus aliases\n"
        )

    # Connect to the DB,
    Session = connect_to_db(db_url, echo=False)  # echo=debug
    session = Session()

    g_old = 0
    g_new = 0
    for taxid in genus_list:
        genus = names[taxid]
        # Is genus already there? e.g. prior import
        taxonomy = (
            session.query(Taxonomy)
            .filter_by(genus=genus, species="", ncbi_taxid=taxid)
            .one_or_none()
        )
        if taxonomy is None:
            # Is genus there without an NCBI taxid?
            taxonomy = (
                session.query(Taxonomy).filter_by(genus=genus, species="").one_or_none()
            )
            if taxonomy and taxonomy.ncbi_taxid != taxid:
                # Prior entry had missing/different taxid, must update it
                if debug or taxonomy.ncbi_taxid != 0:
                    sys.stderr.write(
                        f"WARNING: {genus} {''},"
                        f" updating NCBI taxid {taxonomy.ncbi_taxid} -> {taxid}\n"
                    )
                taxonomy.ncbi_taxid = taxid
                session.add(taxonomy)
        if taxonomy is None:
            g_new += 1
            taxonomy = Taxonomy(genus=genus, species="", ncbi_taxid=taxid)
            session.add(taxonomy)
        else:
            g_old += 1

    old = 0
    new = 0
    s_old = 0
    s_new = 0
    for taxid, genus_taxid in genus_species:
        genus = names[genus_taxid]
        species = names[taxid]
        aliases = set(synonyms.get(taxid, []))

        first_word = species.split(" ", 1)[0]
        assert not reject_name(species), species
        if first_word == genus:
            # We are storing "Phytophthora infestans" as
            # genus="Phytophthora", species="infestans"
            species = species.split(" ", 1)[1]

        if f" x {genus} " in species:
            # Want to turn "Phytophthora medicaginis x Phytophthora cryptogea"
            # into "Phytophthora medicaginis x cryptogea", and then take as the
            # species just "medicaginis x cryptogea" (as in older NCBI taxonomy)
            aliases.add(genus + " " + species)
            species = species.replace(f" x {genus} ", " x ")

        # Is species already there? e.g. prior import
        taxonomy = (
            session.query(Taxonomy)
            .filter_by(genus=genus, species=species)
            .one_or_none()
        )
        if taxonomy and taxonomy.ncbi_taxid != taxid:
            # Prior entry had missing/different taxid, must update it
            if debug or taxonomy.ncbi_taxid != 0:
                sys.stderr.write(
                    f"WARNING: {genus} {species},"
                    f" updating NCBI taxid {taxonomy.ncbi_taxid} -> {taxid}\n"
                )
            taxonomy.ncbi_taxid = taxid
            session.add(taxonomy)
        if taxonomy is None:
            new += 1
            taxonomy = Taxonomy(genus=genus, species=species, ncbi_taxid=taxid)
            session.add(taxonomy)
        else:
            old += 1

        for name in sorted(aliases):
            minor_species.append((taxid, name))

    del genus_species, ranks

    # Record aliases/synonyms
    # Treat species under 'unclassified GenusX' as aliases for 'GenusX'
    for taxid, name in sorted(minor_species):
        # Would this actually be useful? e.g. genus matches first word
        if reject_name(name) or name.split(None, 1)[0] == names[taxid]:
            continue
        # Is target already there?
        taxonomy = session.query(Taxonomy).filter_by(ncbi_taxid=taxid).one_or_none()
        assert taxonomy is not None
        assert taxonomy.id is not None, taxonomy
        synonym = session.query(Synonym).filter_by(name=name).one_or_none()
        if synonym is None:
            synonym = Synonym(taxonomy_id=taxonomy.id, name=name)
            session.add(synonym)
            s_new += 1
        else:
            if debug:
                # See this with Hyaloperonospora parasitica species group,
                # treated as a species with entries under it already set as aliases
                sys.stderr.write(
                    f"Minor species {name!r} -> NCBI {genus_taxid},"
                    f" table {taxonomy.id} pre-existing -> {synonym.taxonomy_id}\n"
                )
            s_old += 1

    del names, synonyms

    session.commit()

    sys.stderr.write(
        f"Loaded {g_new} new genera, and {new} new species entries"
        f" with {s_new} synonyms into DB"
        f" ({g_old}, {old} and {s_old} already there)\n"
    )

    return 0
