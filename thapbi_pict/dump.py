"""Dumping out ITS1 database to text files.

This implementes the ``thapbi_pict dump ...`` command.
"""

import sys

from sqlalchemy.orm import aliased, contains_eager

from .db_orm import ITS1, SequenceSource, Taxonomy, connect_to_db


def main(db_url, output_filename, output_format,
         clade="", genus="", species="", debug=True):
    """Run the database dump with arguments from the command line."""
    # Connect to the DB,
    Session = connect_to_db(db_url, echo=debug)
    session = Session()

    entry_count = 0

    if output_filename == "-":
        out_handle = sys.stdout
    else:
        out_handle = open(output_filename, "w")

    # Doing a join to pull in the ITS1 and Taxonomy tables too:
    cur_tax = aliased(Taxonomy)
    its1_seq = aliased(ITS1)
    view = session.query(SequenceSource).join(
        its1_seq, SequenceSource.its1).join(
        cur_tax, SequenceSource.current_taxonomy).options(
        contains_eager(SequenceSource.its1, alias=its1_seq)).options(
        contains_eager(SequenceSource.current_taxonomy, alias=cur_tax))
    # Sorting for reproducibility
    view = view.order_by(SequenceSource.id)

    clade_list = []
    if clade:
        # Split on commas, convert "-" into "" meaning no entry
        clade_list = ["" if _ == "-" else _ for _ in clade.split(",")]
        for x in clade_list:
            if not session.query(Taxonomy).filter_by(clade=x).count():
                sys.stderr.write("WARNING: Clade %r not in database\n" % x)
        view = view.filter(cur_tax.clade.in_(clade_list))
    del clade

    genus_list = []
    if genus.strip():
        # Split on commas, strip white spaces
        genus_list = [_.strip() for _ in genus.strip().split(",")]
        for x in genus_list:
            if not session.query(Taxonomy).filter_by(genus=x).count():
                sys.stderr.write("WARNING: Genus %r not in database\n" % x)
        view = view.filter(cur_tax.genus.in_(genus_list))
    del genus

    sp_list = []
    if species.strip():
        # Split on commas, strip white spaces
        sp_list = [_.strip() for _ in species.strip().split(",")]
        view = view.filter(cur_tax.species.in_(sp_list))
    del species

    if sp_list and len(genus_list) != 1:
        # This is to avoid ambiguity as some species names are used
        # in more than one genus.
        sys.exit("Using -s/--species requires a single genus via -g/--genus\n")

    for x in sp_list:
        if not session.query(Taxonomy).filter_by(species=x,
                                                 genus=genus_list[0]).count():
            sys.stderr.write(
                "WARNING: '%s %s' not in database\n" % (genus_list[0], x))

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
        if genus_list:
            assert seq_source.current_taxonomy.genus in genus_list, \
                seq_source.current_taxonomy
        if sp_list:
            assert seq_source.current_taxonomy.species in sp_list, \
                seq_source.current_taxonomy

    if output_filename == "-":
        sys.stderr.write("Wrote %i %s format entries\n"
                         % (entry_count, output_format))
    else:
        out_handle.close()
        sys.stderr.write("Wrote %i %s format entries to %r\n"
                         % (entry_count, output_format, output_filename))
