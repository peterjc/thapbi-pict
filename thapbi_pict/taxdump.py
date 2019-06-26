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

from .db_orm import Taxonomy
from .db_orm import connect_to_db


def load_nodes(nodes_dmp):
    """Load the NCBI taxdump nodes.dmp file.

    Returns two dicts, the parent/child relationships, and the rank of
    each node.
    """
    tree = {}
    ranks = {}
    with open(nodes_dmp) as handle:
        for line in handle:
            parts = line.split("\t|\t", 3)
            taxid = int(parts[0].strip())
            parent = int(parts[1].strip())
            rank = parts[2].strip()
            tree[taxid] = parent
            ranks[taxid] = rank
    return tree, ranks


def load_names(names_dmp):
    """Load scientific names of species from NCBI taxdump names.dmp file."""
    names = {}
    with open(names_dmp) as handle:
        for line in handle:
            if not line.endswith("\t|\tscientific name\t|\n"):
                continue
            parts = line.split("\t|\t", 3)
            taxid = int(parts[0].strip())
            sci_name = parts[1].strip()
            names[taxid] = sci_name
    return names


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

    tree, ranks = load_nodes(os.path.join(tax, "nodes.dmp"))
    if debug:
        sys.stderr.write("Loaded %i nodes from nodes.dmp\n" % len(tree))

    names = load_names(os.path.join(tax, "names.dmp"))
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
    Session = connect_to_db(db_url, echo=debug)
    session = Session()

    g_old = 0
    g_new = 0
    for taxid in genus_list:
        genus = names[taxid]
        # Is is already there? e.g. prior import
        taxonomy = (
            session.query(Taxonomy)
            .filter_by(clade="", genus=genus, species="", ncbi_taxid=taxid)
            .one_or_none()
        )
        if taxonomy is None:
            g_new += 1
            taxonomy = Taxonomy(clade="", genus=genus, species="", ncbi_taxid=taxid)
            session.add(taxonomy)
        else:
            g_old += 1

    old = 0
    new = 0
    for taxid, genus, species in genus_species:
        clade = ""

        if species.split(" ", 1)[0] == genus:
            # We are storing "Phytophthora infestans" as
            # genus="Phytophthora", species="infestans"
            species = species.split(" ", 1)[1]

        if species == "unclassified " + genus:
            # Another special case
            species = "unclassified"

        # Is is already there? e.g. prior import
        taxonomy = (
            session.query(Taxonomy)
            .filter_by(clade=clade, genus=genus, species=species, ncbi_taxid=taxid)
            .one_or_none()
        )
        if taxonomy is None:
            new += 1
            taxonomy = Taxonomy(
                clade=clade, genus=genus, species=species, ncbi_taxid=taxid
            )
            session.add(taxonomy)
        else:
            old += 1

    session.commit()

    sys.stderr.write(
        "Loaded %i new genera, and %i new species entries into DB "
        "(%i and %i already there)\n" % (g_new, new, g_old, old)
    )

    return 0
