"""Classifying prepared ITS1 reads using an ITS1 database.

This implementes the ``thapbi_pict classify-reads ...`` command.
"""

import os
import sys

from .db_orm import connect_to_db


def method_identity(filename, session, debug=False):
    """Classify using perfect identity.

    This is a deliberately simple approach, in part for testing
    purposes.
    """
    return 0


methods = {
    "identity": method_identity
}


def find_fasta_files(filenames_or_folders, ext=".fasta", debug=False):
    """Interpret a list of filenames and/or foldernames."""
    answer = []
    for x in filenames_or_folders:
        if os.path.isdir(x):
            if debug:
                sys.stderr.write("Walking directory %r\n" % x)
            for f in os.listdir(x):
                if f.endswith(ext):
                    # Check not a directory?
                    answer.append(os.path.join(x, f))
        elif os.path.isfile(x):
            answer.append(x)
        else:
            sys.exit("ERROR: %r is not a file or a directory\n" % x)
    # Warn if there were duplicates?
    return sorted(set(answer))


def main(fasta, db_url, method, debug=False):
    """Implement the thapbi_pict classify-reads command."""
    assert isinstance(fasta, list)

    if method not in methods:
        sys.exit(
            "ERROR: Invalid method name %r, should be one of: %s\n"
            % (method, ", ".join(sorted(methods))))
    method_fn = methods[method]

    fasta_files = find_fasta_files(fasta, debug=debug)
    if debug:
        sys.stderr.write(
            "Classifying %i input FASTA files\n" % len(fasta_files))

    # Connect to the DB,
    Session = connect_to_db(db_url, echo=debug)
    session = Session()

    read_count = 0
    for filename in fasta_files:
        if debug:
            sys.stderr.write("%s classifer on %s\n" % (method, filename))
        read_count += method_fn(filename, session, debug)

    sys.stderr.write(
        "%s classifier ran on %i reads from %i files\n"
        % (method, read_count, len(fasta_files)))

    return 0
