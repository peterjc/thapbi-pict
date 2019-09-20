# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

"""Dumping out ITS1 database to text files.

This implements the ``thapbi_pict dump ...`` command.
"""

import sys

from sqlalchemy.orm import aliased, contains_eager

from .db_orm import ITS1, SequenceSource, Taxonomy, connect_to_db

from .utils import genus_species_name

from .classify import taxid_and_sp_lists


def none_str(value, none_value=""):
    """Turn value into a string, special case None to empty string."""
    if value is None:
        return none_value
    else:
        return str(value)


def main(
    db_url,
    output_filename,
    output_format,
    minimal=False,
    genus="",
    species="",
    debug=True,
):
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
    view = (
        session.query(SequenceSource)
        .join(its1_seq, SequenceSource.its1)
        .join(cur_tax, SequenceSource.current_taxonomy)
        .options(contains_eager(SequenceSource.its1, alias=its1_seq))
        .options(contains_eager(SequenceSource.current_taxonomy, alias=cur_tax))
    )
    # Sorting for reproducibility
    view = view.order_by(its1_seq.sequence, SequenceSource.id)

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
        sys.exit("ERROR: Using -s/--species requires a single genus via -g/--genus\n")

    for x in sp_list:
        if (
            not session.query(Taxonomy)
            .filter_by(species=x, genus=genus_list[0])
            .count()
        ):
            sys.stderr.write(
                "WARNING: '%s' not in database\n" % genus_species_name(genus_list[0], x)
            )

    if output_format == "fasta":
        # no header
        pass
    elif minimal:
        out_handle.write("#MD5\tSpecies\tSequence\n")
    else:
        out_handle.write(
            "#Identifier\tGenus\tSpecies\tTaxID\tITS1-MD5\tITS1-seq\tSequence\n"
        )

    if minimal:
        md5_seq = {}
        md5_sp = {}
        for seq_source in view:
            md5 = seq_source.its1.md5
            # genus_species = genus_species_name(
            #                            seq_source.current_taxonomy.genus,
            #                            seq_source.current_taxonomy.species,
            #                            )
            if md5 in md5_seq:
                assert md5_seq[md5] == seq_source.its1.sequence
                md5_sp[md5].add(seq_source.current_taxonomy)
            else:
                md5_seq[md5] = seq_source.its1.sequence
                md5_sp[md5] = set([seq_source.current_taxonomy])  # noqa: C405
        for md5, seq in md5_seq.items():
            _, genus_species, _ = taxid_and_sp_lists(list(md5_sp[md5]))
            if output_format == "fasta":
                template = ">%s %s\n%s\n"
            else:
                template = "%s\t%s\t%s\n"
            try:
                out_handle.write(template % (md5, genus_species, seq))
            except BrokenPipeError:
                # Likely writing to stdout | head, or similar
                # If so, stdout has been closed
                sys.stderr.write("Aborting with broken pipe\n")
                sys.stderr.close()
                sys.exit(1)
        entry_count = len(md5_seq)
    else:
        for seq_source in view:
            entry_count += 1
            taxid = (
                str(seq_source.current_taxonomy.ncbi_taxid)
                if seq_source.current_taxonomy.ncbi_taxid
                else ""
            )
            try:
                if output_format == "fasta":
                    genus_species = genus_species_name(
                        seq_source.current_taxonomy.genus,
                        seq_source.current_taxonomy.species,
                    )
                    out_handle.write(
                        ">%s [species=%s] [taxid=%s]\n%s\n"
                        % (
                            seq_source.source_accession,
                            genus_species,
                            taxid,
                            seq_source.its1.sequence,
                        )
                    )
                else:
                    out_handle.write(
                        "%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                        % (
                            seq_source.source_accession,
                            none_str(seq_source.current_taxonomy.genus),
                            none_str(seq_source.current_taxonomy.species),
                            taxid,
                            seq_source.its1.md5,
                            seq_source.its1.sequence,
                            seq_source.sequence,
                        )
                    )
            except BrokenPipeError:
                # Likely writing to stdout | head, or similar
                # If so, stdout has been closed
                sys.stderr.write("Aborting with broken pipe\n")
                sys.stderr.close()
                sys.exit(1)

            if genus_list:
                assert (
                    seq_source.current_taxonomy.genus in genus_list
                ), seq_source.current_taxonomy
            if sp_list:
                assert (
                    seq_source.current_taxonomy.species in sp_list
                ), seq_source.current_taxonomy

    if output_filename == "-":
        sys.stderr.write("Wrote %i %s format entries\n" % (entry_count, output_format))
    else:
        out_handle.close()
        sys.stderr.write(
            "Wrote %i %s format entries to %r\n"
            % (entry_count, output_format, output_filename)
        )
