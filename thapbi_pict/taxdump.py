# Copyright 2018-2020 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Code for THAPBI PICT to deal with NCBI taxonomy dumps.

The code is needed initially for loading an NCBI taxdump folder (files
``names.dmp``, ``nodes.dmp`` etc) into an ITS1 database.
"""
import os
import sys

from .db_orm import connect_to_db
from .db_orm import Synonym
from .db_orm import Taxonomy


def load_nodes(nodes_dmp):
    """Load the NCBI taxdump nodes.dmp file.

    Returns three dicts, the parent/child relationships (two directions),
    and the rank of each node.
    """
    tree = {}
    ranks = {}  # child points at parent
    descendants = {}  # parent points at children and grandchildren etc
    with open(nodes_dmp) as handle:
        for line in handle:
            parts = line.split("\t|\t", 3)
            taxid = int(parts[0].strip())
            parent = int(parts[1].strip())
            rank = parts[2].strip()
            tree[taxid] = parent
            ranks[taxid] = rank
            try:
                descendants[parent].add(taxid)
            except KeyError:
                descendants[parent] = {taxid}
    return tree, descendants, ranks


def load_names(names_dmp):
    """Load scientific names of species from NCBI taxdump names.dmp file."""
    names = {}
    synonym = {}
    with open(names_dmp) as handle:
        for line in handle:
            parts = line.split("\t|\t", 3)
            taxid = int(parts[0].strip())
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


def check_ancestor(tree, ancestors, taxid):
    """Return True taxid is descended from any of the ancestors.

    Argument ancestors is a list of integer taxids.
    """
    t = taxid
    while t not in ancestors and t != tree[t]:
        t = tree[t]
    return t in ancestors


def genera_under_ancestors(tree, ranks, ancestors):
    """Return genus taxid found under ancestors."""
    for taxid, rank in ranks.items():
        if rank == "genus" and check_ancestor(tree, ancestors, taxid):
            yield taxid


def get_children(children, taxid):
    """Return all taxid descended from given entry."""
    answer = set()
    for child in children.get(taxid, []):
        answer.add(child)
        answer.update(get_children(children, child))  # recurse!
    return answer


def synonyms_and_variants(children, ranks, names, synonyms, species_taxid):
    """Return all scientific names and synonyms of any variants etc under species."""
    assert ranks[species_taxid] in ("species", "species group", "no rank"), ranks[
        species_taxid
    ]
    variants = set(synonyms.get(species_taxid, []))  # include own synonyms
    for taxid in get_children(children, species_taxid):
        variants.add(names[taxid])
        variants.update(synonyms.get(taxid, []))  # values are already names
    return sorted(variants)


def top_level_species(children, ranks, names, genus_list):
    """Find taxids for species under the genus_list.

    Rather than just taking all taxids with rank species, we
    are taking all the immediate children of the genus - this
    is specifically to treat "unclassified Phytophthora"
    (taxid 211524), which has no rank, as a species level ID
    AND to ignore all the species rank entries under it.
    """
    for genus_taxid in genus_list:
        assert ranks[genus_taxid] in ("genus", "subgenus")
        if genus_taxid not in children:
            sys.stderr.write(
                f"WARNING: Genus {names[genus_taxid]} ({genus_taxid})"
                " has no children\n"
            )
            continue
        for taxid in children[genus_taxid]:
            name = names[taxid]
            if ranks[taxid] == "species":
                yield taxid, names[genus_taxid], name
            elif ranks[taxid] == "subgenus":
                if taxid in children:
                    sys.stderr.write(
                        f"WARNING: Collapsing sub-genus {name} into parent\n"
                    )
                    for _ in top_level_species(children, ranks, names, [taxid]):
                        yield _
                else:
                    sys.stderr.write(
                        f"WARNING: Ignoring sub-genus {name} with no children\n"
                    )
            else:
                if name.split(None, 1)[0].lower() in ("unclassified", "unidentified"):
                    # Not worth including, nor giving a warning about
                    continue
                if name.lower() in ("environmental samples"):
                    # Again silently ignore
                    continue
                # e.g. Hyaloperonospora parasitica species group
                sys.stderr.write(
                    f"WARNING: Treating {ranks[taxid]} '{name}' (txid{taxid})"
                    " as a species.\n"
                )
                yield taxid, names[genus_taxid], name


def not_top_species(top_species, children, ranks, names, genus_list):
    """Find all 'minor' species, takes set of species taxid to ignore.

    Intended usage is to map minor-species as genus aliases, for
    example all the species under unclassified Phytophthora will
    be treated as synonyms of the genus Phytophthora.
    """
    for genus_taxid in genus_list:
        if genus_taxid not in children:
            continue
        for taxid in children[genus_taxid]:
            if ranks[taxid] == "species":
                if taxid not in top_species:
                    yield genus_taxid, names[taxid]
            if taxid in children:
                # recurse...
                for _ in not_top_species(top_species, children, ranks, names, [taxid]):
                    yield genus_taxid, _[1]


def main(tax, db_url, ancestors, debug=True):
    """Load an NCBI taxdump into an ITS1 database."""
    if not os.path.isdir(tax):
        sys.exit(f"ERROR: Could not find taxdump directory: {tax!r}\n")
    for filename in ("names.dmp", "nodes.dmp"):
        if not os.path.isfile(os.path.join(tax, filename)):
            sys.exit(f"ERROR: Missing {filename} in the taxdump directory: {tax!r}\n")

    try:
        ancestors = [int(_) for _ in ancestors.split(",")]
    except ValueError:
        sys.exit(f"ERROR: Invalid ancestors argument: {ancestors!r}\n")

    tree, children, ranks = load_nodes(os.path.join(tax, "nodes.dmp"))
    if debug:
        sys.stderr.write(f"Loaded {len(tree)} nodes from nodes.dmp\n")

    names, synonyms = load_names(os.path.join(tax, "names.dmp"))
    if debug:
        sys.stderr.write(f"Loaded {len(names)} scientific names from names.dmp\n")

    genus_list = list(genera_under_ancestors(tree, ranks, ancestors))
    if not genus_list:
        sys.exit("ERROR: Could not identify any genus names under the given nodes\n")
    if debug:
        sys.stderr.write(
            f"Identified {len(genus_list)} genera under specified ancestor node:"
            f" {', '.join(sorted(names[_] for _ in genus_list))}\n"
        )

    genus_species = list(top_level_species(children, ranks, names, genus_list))
    if debug:
        sys.stderr.write(f"Filtered down to {len(genus_species)} species names\n")

    minor_species = list(
        not_top_species(
            {_[0] for _ in genus_species}, children, ranks, names, genus_list
        )
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
    for taxid, genus, species in genus_species:
        aliases = []

        if species.split(" ", 1)[0] == genus:
            # We are storing "Phytophthora infestans" as
            # genus="Phytophthora", species="infestans"
            species = species.split(" ", 1)[1]

        if f" x {genus} " in species:
            # Want to turn "Phytophthora medicaginis x Phytophthora cryptogea"
            # into "Phytophthora medicaginis x cryptogea", and then take as the
            # species just "medicaginis x cryptogea" (as in older NCBI taxonomy)
            aliases.append(genus + " " + species)
            species = species.replace(f" x {genus} ", " x ")

        if species == "unclassified " + genus:
            # Another special case
            species = "unclassified"

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

        for name in (
            synonyms_and_variants(children, ranks, names, synonyms, taxid) + aliases
        ):
            # Is it already there?
            synonym = session.query(Synonym).filter_by(name=name).one_or_none()
            if synonym is None:
                synonym = Synonym(taxonomy_id=taxonomy.id, name=name)
                session.add(synonym)
                s_new += 1
            elif name in aliases:
                # Don't double count it
                pass
            else:
                s_old += 1

    # Treat species under 'unclassified GenusX' as aliases for 'GenusX'
    for genus_taxid, name in minor_species:
        # Is it already there?
        taxonomy = (
            session.query(Taxonomy).filter_by(ncbi_taxid=genus_taxid).one_or_none()
        )
        assert taxonomy is not None
        assert taxonomy.id is not None, taxonomy
        synonym = session.query(Synonym).filter_by(name=name).one_or_none()
        if synonym is None:
            synonym = Synonym(taxonomy_id=taxonomy.id, name=name)
            session.add(synonym)
            s_new += 1
        elif name in aliases:
            # Don't double count it
            pass
        else:
            if debug:
                # See this with Hyaloperonospora parasitica species group,
                # treated as a species with entries under it already set as aliases
                sys.stderr.write(
                    f"Minor species {name!r} -> NCBI {genus_taxid},"
                    f" table {taxonomy.id} pre-existing -> {synonym.taxonomy_id}\n"
                )
            s_old += 1

    session.commit()

    sys.stderr.write(
        f"Loaded {g_new} new genera, and {new} new species entries"
        f" with {s_new} synonyms into DB"
        f" ({g_old}, {old} and {s_old} already there)\n"
    )

    return 0
