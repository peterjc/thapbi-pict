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
from .db_orm import MarkerDef
from .db_orm import MarkerSeq
from .db_orm import SeqSource
from .db_orm import Taxonomy
from .utils import genus_species_name


def main(db_url, output_filename, debug=False):
    """Implement the ``thapbi_pict conflicts`` subcommand.

    Looks for taxonomy conflicts at marker, genus or species level, with the
    number of marker or genus level conflicts used as the return code. i.e.
    Unix failure (non-zero) when there are marker or genus level conflicts.

    A marker level conflict is when a unique sequence appears in the DB under
    more than one marker name (e.g. both COI and ITS1), which is most likely
    an error in the DB construction.

    Genus level conflicts are where a unique sequence in the DB is reported
    from more than one genus, which is considered undesirable. Similarly for
    species level conflicts, but for some markers this is sadly common and not
    considered to be an error.
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
    marker_seq = aliased(MarkerSeq)
    marker_def = aliased(MarkerDef)
    view = (
        session.query(SeqSource)
        .join(marker_seq, SeqSource.marker_seq)
        .join(marker_def, SeqSource.marker_definition)
        .join(cur_tax, SeqSource.taxonomy)
        .options(contains_eager(SeqSource.marker_seq, alias=marker_seq))
        .options(contains_eager(SeqSource.marker_definition, alias=marker_def))
        .options(contains_eager(SeqSource.taxonomy, alias=cur_tax))
    )
    md5_to_seq = {}
    md5_to_marker = {}
    md5_to_genus = {}
    md5_to_species = {}
    for seq_source in view:
        md5 = seq_source.marker_seq.md5
        seq = seq_source.marker_seq.sequence
        genus = seq_source.taxonomy.genus
        md5_to_seq[md5] = seq
        try:
            md5_to_marker[md5].add(seq_source.marker_definition.name)
        except KeyError:
            md5_to_marker[md5] = {seq_source.marker_definition.name}
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

    if debug:
        sys.stderr.write(f"Loaded taxonomy for {len(md5_to_seq)} sequences from DB\n")

    marker_conflicts = 0
    genus_conflicts = 0
    out_handle.write("#MD5\tLevel\tConflicts\n")
    for md5, markers in sorted(md5_to_marker.items()):
        if len(markers) > 1:
            out_handle.write(f"{md5}\tmarker\t{';'.join(sorted(markers))}\n")
            marker_conflicts += 1
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
        sys.stderr.write(f"{marker_conflicts} marker level conflicts\n")
        sys.stderr.write(f"{genus_conflicts} genus level conflicts\n")

    return marker_conflicts + genus_conflicts  # non-zero error
