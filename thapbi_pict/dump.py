"""Dumping out ITS1 database to text files.

This implementes the ``thapbi_pict dump ...`` command.
"""

import sys

from sqlalchemy.orm import joinedload, aliased

from .db_orm import SequenceSource, Taxonomy, connect_to_db


def main(db_url, output_txt, clade="", debug=True):
    """Run the database dump with arguments from the command line."""
    # Connect to the DB,
    Session = connect_to_db(db_url, echo=debug)
    session = Session()

    entry_count = 0

    if output_txt == "-":
        out_handle = sys.stdout
    else:
        out_handle = open(output_txt, "w")

    # Doing a join to pull in the ITS1 and Taxonomy tables too:
    cur_tax = aliased(Taxonomy)
    view = session.query(SequenceSource).join(
        cur_tax, SequenceSource.current_taxonomy).options(
        joinedload(SequenceSource.current_taxonomy)).options(
        joinedload(SequenceSource.its1))
    # Sorting for reproducibility
    view = view.order_by(SequenceSource.id)

    clade_list = None
    if clade:
        # Split on commas, convert "-" into "" meaning no entry
        clade_list = ["" if _ == "-" else _ for _ in clade.split(",")]
        view = view.filter(cur_tax.clade.in_(clade_list))

    for seq_source in view:
        entry_count += 1
        out_handle.write("%s\t%s\t%s\t%s\t%s\t%s\n"
                         % (seq_source.source_accession,
                            seq_source.current_taxonomy.clade,
                            seq_source.current_taxonomy.genus,
                            seq_source.current_taxonomy.species,
                            seq_source.its1.sequence,
                            seq_source.sequence))
        if clade_list:
            assert seq_source.current_taxonomy.clade in clade_list, \
                seq_source.current_taxonomy

    if output_txt == "-":
        sys.stderr.write("Wrote %i entries\n" % entry_count)
    else:
        out_handle.close()
        sys.stderr.write("Wrote %i entries to %r\n"
                         % (entry_count, output_txt))
