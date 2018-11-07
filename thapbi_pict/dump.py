"""Dumping out ITS1 database to text files.

This implementes the ``thapbi_pict dump ...`` command.
"""

import sys

from sqlalchemy.orm import joinedload

from .db_orm import SequenceSource, connect_to_db


def main(db_url, output_txt, debug=True):
    """Run the database dump with arguments from the command line."""
    # Connect to the DB,
    Session = connect_to_db(db_url, echo=debug)
    session = Session()

    entry_count = 0

    if output_txt == "-":
        out_handle = sys.stdout
    else:
        out_handle = open(output_txt, "w")

    for seq_source in session.query(SequenceSource).options(
            joinedload(SequenceSource.its1_seq)).order_by(
            SequenceSource.date_added):
        entry_count += 1
        out_handle.write("%s\t%s\t%s\t%s\n"
                         % (seq_source.accession,
                            seq_source.current_species,
                            seq_source.its1_seq.sequence,
                            seq_source.sequence))

    if output_txt == "-":
        sys.stderr.write("Wrote %i entries\n" % entry_count)
    else:
        out_handle.close()
        sys.stderr.write("Wrote %i entries to %r\n"
                         % (entry_count, output_txt))
