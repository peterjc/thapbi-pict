"""Explore conflicts at species and genus level."""

import sys

from sqlalchemy.orm import aliased
from sqlalchemy.orm import contains_eager

from thapbi_pict.db_orm import connect_to_db
from thapbi_pict.db_orm import ITS1
from thapbi_pict.db_orm import SequenceSource
from thapbi_pict.db_orm import Taxonomy
from thapbi_pict.utils import genus_species_name


def main(db_url, output_filename, debug=False):
    """Implement the conflicts subcommand.

    Look for taxonomy conflicts at genus or species level, with the
    number of genus level conflicts used as the return code. i.e.
    Unix failure (non-zero) when there are genus conflicts.
    """
    if output_filename == "-":
        out_handle = sys.stdout
    else:
        out_handle = open(output_filename, "w")

    # Connect to the DB,
    Session = connect_to_db(db_url, echo=False)
    session = Session()

    # Doing a join to pull in the ITS1 and Taxonomy tables too:
    cur_tax = aliased(Taxonomy)
    its1_seq = aliased(ITS1)
    view = (
        session.query(SequenceSource)
        .join(its1_seq, SequenceSource.its1)
        .join(cur_tax, SequenceSource.current_taxonomy)
        .options(contains_eager(SequenceSource.its1, alias=its1_seq))
        .options(contains_eager(SequenceSource.current_taxonomy, alias=cur_tax))
    )
    md5_to_seq = {}
    md5_to_genus = {}
    md5_to_species = {}
    for seq_source in view:
        md5 = seq_source.its1.md5
        seq = seq_source.its1.sequence
        genus = seq_source.current_taxonomy.genus
        md5_to_seq[md5] = seq
        if genus:
            try:
                md5_to_genus[md5].add(genus)
            except KeyError:
                md5_to_genus[md5] = {genus}
            if seq_source.current_taxonomy.species:
                genus_species = genus_species_name(
                    genus, seq_source.current_taxonomy.species
                )
                try:
                    md5_to_species[md5].add(genus_species)
                except KeyError:
                    md5_to_species[md5] = {genus_species}

    sys.stderr.write("Loaded taxonomy for %i sequences from DB\n" % len(md5_to_seq))

    genus_conflicts = 0
    out_handle.write("#MD5\tLevel\tConflicts\n")
    for md5, genus in md5_to_genus.items():
        if len(genus) > 1:
            out_handle.write("%s\t%s\t%s\n" % (md5, "genus", ";".join(sorted(genus))))
            genus_conflicts += 1
    for md5, species in md5_to_species.items():
        if len(species) > 1:
            out_handle.write(
                "%s\t%s\t%s\n" % (md5, "species", ";".join(sorted(species)))
            )

    if output_filename != "-":
        out_handle.close()

    if debug:
        sys.stderr.write("%i genus level conflicts\n" % genus_conflicts)

    return genus_conflicts
