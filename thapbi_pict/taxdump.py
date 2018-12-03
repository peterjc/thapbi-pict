"""Code for THAPBI PICT to deal with NCBI taxonomy dumps.

The code is needed initially for loading an NCBI taxdump folder (files
``names.dmp``, ``nodes.dmp`` etc) into an ITS1 database.
"""

import os
import sys


def load_nodes(nodes_dmp):
    """Load parent/child relationships from NCBI taxdump nodes.dmp file."""
    tree = dict()
    with open(nodes_dmp) as handle:
        for line in handle:
            parts = line.split("\t|\t", 2)
            taxid = int(parts[0].strip())
            parent = int(parts[1].strip())
            tree[taxid] = parent
    return tree


def load_names(names_dmp):
    """Load scientific names of species from NCBI taxdump names.dmp file."""
    names = dict()
    with open(names_dmp) as handle:
        for line in handle:
            if not line.endswith("\t|\tscientific name\t|\n"):
                continue
            parts = line.split("\t|\t", 3)
            taxid = int(parts[0].strip())
            sci_name = parts[1].strip()
            names[taxid] = sci_name
    return names


def main(tax, db_url, debug=True):
    """Load an NCBI taxdump into an ITS1 database."""
    if not os.path.isdir(tax):
        sys.exit("Could not find taxdump directory: %r\n" % tax)
    for filename in ("names.dmp", "nodes.dmp"):
        if not os.path.isfile(os.path.join(tax, filename)):
            sys.exit("Missing %s in the taxdump directory: %r\n"
                     % (filename, tax))

    tree = load_nodes(os.path.join(tax, "nodes.dmp"))
    if debug:
        sys.stderr.write("Loaded %i nodes from nodes.dmp\n" % len(tree))

    names = load_names(os.path.join(tax, "names.dmp"))
    if debug:
        sys.stderr.write("Loaded %i relevant scientific names from names.dmp\n"
                         % len(names))

    return 0
