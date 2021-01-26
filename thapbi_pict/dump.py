# Copyright 2018-2020 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Dumping out ITS1 database to text files.

This implements the ``thapbi_pict dump ...`` command.
"""
import sys

from sqlalchemy.orm import aliased
from sqlalchemy.orm import contains_eager

from .classify import taxid_and_sp_lists
from .db_orm import connect_to_db
from .db_orm import ITS1
from .db_orm import SequenceSource
from .db_orm import Taxonomy
from .utils import genus_species_name


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
    sep=None,
    debug=True,
):
    """Run the database dump with arguments from the command line."""
    if not sep:
        # TODO - Use this argument for tab/comma/etc in txt output?
        sep = ";"

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
                sys.stderr.write(f"WARNING: Genus {x!r} not in database\n")
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
                f"WARNING: '{genus_species_name(genus_list[0], x)}' not in database\n"
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
    elif output_format == "fasta":
        seq_entry = {}
        for seq_source in view:
            seq = seq_source.its1.sequence
            entry = genus_species_name(
                seq_source.current_taxonomy.genus, seq_source.current_taxonomy.species
            )
            entry = f"{seq_source.source_accession} {entry}"
            if seq in seq_entry:
                seq_entry[seq].add(entry)
            else:
                seq_entry[seq] = {entry}
        for seq, entries in sorted(seq_entry.items()):
            # entry_count += len(entries)
            try:
                out_handle.write(f">{sep.join(sorted(entries))}\n{seq}\n")
            except BrokenPipeError:
                # Likely writing to stdout | head, or similar
                # If so, stdout has been closed
                sys.stderr.write("Aborting with broken pipe\n")
                sys.stderr.close()
                sys.exit(1)
        entry_count = len(seq_entry)  # number of FASTA entries
    else:
        for seq_source in view:
            entry_count += 1
            taxid = (
                str(seq_source.current_taxonomy.ncbi_taxid)
                if seq_source.current_taxonomy.ncbi_taxid
                else ""
            )
            try:
                out_handle.write(
                    f"{seq_source.source_accession}"
                    f"\t{none_str(seq_source.current_taxonomy.genus)}"
                    f"\t{none_str(seq_source.current_taxonomy.species)}"
                    f"\t{taxid}"
                    f"\t{seq_source.its1.md5}"
                    f"\t{seq_source.its1.sequence}"
                    f"\t{seq_source.sequence}\n"
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
        sys.stderr.write(f"Wrote {entry_count} {output_format} format entries\n")
    else:
        out_handle.close()
        sys.stderr.write(
            f"Wrote {entry_count} {output_format} format entries to"
            f" {output_filename!r}\n"
        )
