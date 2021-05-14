# Copyright 2019-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Explore conflicts at species and genus level."""
import sys

from sqlalchemy.orm import aliased
from sqlalchemy.orm import contains_eager

from .db_orm import connect_to_db
from .db_orm import RefMarker
from .db_orm import SequenceSource
from .db_orm import Taxonomy
from .utils import genus_species_name


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

    # Doing a join to pull in the marker and taxonomy tables too:
    cur_tax = aliased(Taxonomy)
    marker_seq = aliased(RefMarker)
    view = (
        session.query(SequenceSource)
        .join(marker_seq, SequenceSource.marker)
        .join(cur_tax, SequenceSource.taxonomy)
        .options(contains_eager(SequenceSource.marker, alias=marker_seq))
        .options(contains_eager(SequenceSource.taxonomy, alias=cur_tax))
    )
    md5_to_seq = {}
    md5_to_genus = {}
    md5_to_species = {}
    for seq_source in view:
        md5 = seq_source.marker.md5
        seq = seq_source.marker.sequence
        genus = seq_source.taxonomy.genus
        md5_to_seq[md5] = seq
        if genus:
            try:
                md5_to_genus[md5].add(genus)
            except KeyError:
                md5_to_genus[md5] = {genus}
            if seq_source.taxonomy.species:
                genus_species = genus_species_name(genus, seq_source.taxonomy.species)
                try:
                    md5_to_species[md5].add(genus_species)
                except KeyError:
                    md5_to_species[md5] = {genus_species}

    sys.stderr.write(f"Loaded taxonomy for {len(md5_to_seq)} sequences from DB\n")

    genus_conflicts = 0
    out_handle.write("#MD5\tLevel\tConflicts\n")
    for md5, genus in sorted(md5_to_genus.items()):
        if len(genus) > 1:
            out_handle.write(f"{md5}\tgenus\t{';'.join(sorted(genus))}\n")
            genus_conflicts += 1
    for md5, species in sorted(md5_to_species.items()):
        if len(species) > 1:
            out_handle.write(f"{md5}\tspecies\t{';'.join(sorted(species))}\n")

    if output_filename != "-":
        out_handle.close()

    if debug:
        sys.stderr.write(f"{genus_conflicts} genus level conflicts\n")

    return genus_conflicts
