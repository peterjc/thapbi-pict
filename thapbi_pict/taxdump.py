"""Code for THAPBI PICT to deal with NCBI taxonomy dumps.

The code is needed initially for loading an NCBI taxdump folder (files
``names.dmp``, ``nodes.dmp`` etc) into an ITS1 database.
"""

import os
import sys


def main(tax, db_url, debug=True):
    """Load an NCBI taxdump into an ITS1 database."""
    if not os.path.isdir(tax):
        sys.exit("Could not find taxdump directory: %r\n" % tax)
    for filename in ("names.dmp", "nodes.dmp"):
        if not os.path.isfile(os.path.join(tax, filename)):
            sys.exit("Missing %s in the taxdump directory: %r\n"
                     % (filename, tax))
    return 0
