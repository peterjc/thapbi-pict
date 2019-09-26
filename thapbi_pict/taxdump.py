# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
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

from .db_orm import Synonym
from .db_orm import Taxonomy
from .db_orm import connect_to_db


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
            elif line.endswith("\t|\tsynonym\t|\n"):
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


def top_level_species(tree, ranks, names, genus_list):
    """Find taxids for species under the genus_list.

    Rather than just taking all taxids with rank species, we
    are taking all the immediate children of the genus - this
    is specifically to treat "unclassified Phytophthora"
    (taxid 211524), which has no rank, as a species level ID
    AND to ignore all the species rank entries under it.
    """
    for taxid, name in names.items():
        if tree[taxid] in genus_list:
            if ranks[taxid] != "species":
                if name.split(None, 1)[0].lower() == "unclassified":
                    # Not worth including, nor giving a warning about
                    continue
                sys.stderr.write(
                    "WARNING: Treating %s '%s' (txid%i) as a species.\n"
                    % (ranks[taxid], name, taxid)
                )
            yield taxid, names[tree[taxid]], name


def main(tax, db_url, ancestors, debug=True):
    """Load an NCBI taxdump into an ITS1 database."""
    if not os.path.isdir(tax):
        sys.exit("ERROR: Could not find taxdump directory: %r\n" % tax)
    for filename in ("names.dmp", "nodes.dmp"):
        if not os.path.isfile(os.path.join(tax, filename)):
            sys.exit(
                "ERROR: Missing %s in the taxdump directory: %r\n" % (filename, tax)
            )

    try:
        ancestors = [int(_) for _ in ancestors.split(",")]
    except ValueError:
        sys.exit("ERROR: Invalid ancestors argument: %r\n" % ancestors)

    tree, children, ranks = load_nodes(os.path.join(tax, "nodes.dmp"))
    if debug:
        sys.stderr.write("Loaded %i nodes from nodes.dmp\n" % len(tree))

    names, synonyms = load_names(os.path.join(tax, "names.dmp"))
    if debug:
        sys.stderr.write("Loaded %i scientific names from names.dmp\n" % len(names))

    genus_list = list(genera_under_ancestors(tree, ranks, ancestors))
    if not genus_list:
        sys.exit("ERROR: Could not identify any genus names under the given nodes\n")
    if debug:
        sys.stderr.write(
            "Identified %i genera under specified ancestor node: %s\n"
            % (len(genus_list), ", ".join(sorted(names[_] for _ in genus_list)))
        )

    genus_species = list(top_level_species(tree, ranks, names, genus_list))
    if debug:
        sys.stderr.write("Filtered down to %i species names\n" % len(genus_species))

    # Connect to the DB,
    Session = connect_to_db(db_url, echo=False)  # echo=debug
    session = Session()

    g_old = 0
    g_new = 0
    for taxid in genus_list:
        genus = names[taxid]
        # Is is already there? e.g. prior import
        taxonomy = (
            session.query(Taxonomy)
            .filter_by(genus=genus, species="", ncbi_taxid=taxid)
            .one_or_none()
        )
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

        if " x %s " % genus in species:
            # Want to turn "Phytophthora medicaginis x Phytophthora cryptogea"
            # into "Phytophthora medicaginis x cryptogea", and then take as the
            # species just "medicaginis x cryptogea" (as in older NCBI taxonomy)
            aliases.append(genus + " " + species)
            species = species.replace(" x %s " % genus, " x ")

        if species == "unclassified " + genus:
            # Another special case
            species = "unclassified"

        # Is is already there? e.g. prior import
        taxonomy = (
            session.query(Taxonomy)
            .filter_by(genus=genus, species=species)
            .one_or_none()
        )
        if taxonomy and taxonomy.ncbi_taxid != taxid:
            # Prior entry had missing/different taxid, must update it
            if debug or taxonomy.ncbi_taxid != 0:
                sys.stderr.write(
                    "WARNING: %s %s, updating NCBI taxid %i -> %i\n"
                    % (genus, species, taxonomy.ncbi_taxid, taxid)
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

    session.commit()

    sys.stderr.write(
        "Loaded %i new genera, and %i new species entries with %i synonyms into DB "
        "(%i, %i and %i already there)\n" % (g_new, new, s_new, g_old, old, s_old)
    )

    return 0
