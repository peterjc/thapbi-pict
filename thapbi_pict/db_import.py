# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

"""Shared code for THAPBI PICT to import FASTA into our database.

This code is used both for importing NCBI formatted FASTA files, and also
importing our legacy ITS1 sequence FASTA file databases - see ``ncbi.py``
and ``legacy.py`` which contain specific meta-data handling code for the
different naming conventions.
"""

import os
import sys

from . import __version__
from .db_orm import DataSource
from .db_orm import ITS1
from .db_orm import SequenceSource
from .db_orm import Synonym
from .db_orm import Taxonomy
from .db_orm import connect_to_db
from .hmm import filter_for_ITS1
from .utils import genus_species_name
from .utils import genus_species_split
from .utils import md5seq
from .utils import md5_hexdigest
from .versions import check_tools


def load_taxonomy(session):
    """Pre-load all the species and synonym names as a set."""
    names = set()
    view = session.query(Taxonomy).distinct(Taxonomy.genus, Taxonomy.species)
    for taxonomy in view:
        names.add(genus_species_name(taxonomy.genus, taxonomy.species))
    for synonym in session.query(Synonym):
        if synonym.name in names:
            sys.stderr.write("WARNING: Synonym %s duplicated?\n" % synonym.name)
        else:
            names.add(synonym.name)
    return names


def lookup_species(session, name):
    """Find this species entry in the taxonomy/synonym table (if present)."""
    assert isinstance(name, str), name
    genus, species = genus_species_split(name)
    # Try main table
    taxonomy = (
        session.query(Taxonomy).filter_by(genus=genus, species=species).one_or_none()
    )
    if taxonomy:
        return taxonomy
    # Try synonyms
    return (
        session.query(Taxonomy).join(Synonym).filter(Synonym.name == name).one_or_none()
    )


def lookup_genus(session, name):
    """Find genus entry via taxonomy/synonym table (if present)."""
    # Apply synonym (which might change the genus)
    taxonomy = (
        session.query(Taxonomy).join(Synonym).filter(Synonym.name == name).one_or_none()
    )
    if taxonomy:
        genus = taxonomy.genus
    else:
        genus = genus_species_split(name)[0]
    return session.query(Taxonomy).filter_by(genus=genus, species="").one_or_none()


def import_fasta_file(
    fasta_file,
    db_url,
    name=None,
    debug=True,
    fasta_entry_fn=None,
    entry_taxonomy_fn=None,
    validate_species=False,
    genus_only=False,
):
    """Import a FASTA file into the database.

    For ``thapbi_pict legacy-import`` some FASTA sequences are
    treated as multiple entries sharing that same sequence. For
    ``thapbi_pict ncbi-import``, each FASTA sequence is treated
    as a single entry. For ``thapbi_pict seq-import`` again each
    FASTA sequence is treated as a single entry, but depending on
    the meta-data it may or may not be imported into the DB.

    That behaviour is controlled by the optional argument
    fasta_entry_fn which is a function which will be given the
    full FASTA title line, and should return a list of sub-entries
    (by default, a single entry list) to be imported. This can
    return an empty list if the FASTA entry is to be ignored.

    Required argument entry_taxonomy_fn is a function which will
    be given the entries from function fasta_split_fn, and should
    return the associated taxonomy information. This can be done
    by parsing the string, or looking it up in another source.

    In ``thapbi_pict legacy-import`` and ``thapbi_pict ncbi-import``
    the species metadata is recorded directly in the FASTA title
    lines. However, the metadata for ``thapbi_pict seq-import``
    comes from a sister TSV file, and is cross-referenced by the
    FASTA sequence identifier.
    """
    check_tools(["hmmscan"], debug)

    # Argument validation,
    if fasta_entry_fn is None:
        if debug:
            sys.stderr.write("DEBUG: Treating each FASTA entry as a singleton.\n")

        def fasta_entry_fn(text):
            """Treat all FASTA entries as singletons.

            Default is not to support merged FASTA entries, and accept
            all the sequences.
            """
            return [text]

    if entry_taxonomy_fn is None:
        raise ValueError("Need function to get species from FASTA title.")

    if os.stat(fasta_file).st_size == 0:
        if debug:
            sys.stderr.write("Ignoring empty FASTA file %s\n" % fasta_file)
        return

    # Connect to the DB,
    Session = connect_to_db(db_url, echo=False)  # echo=debug
    session = Session()

    additional_taxonomy = {}  # any entries added this session
    preloaded_taxonomy = load_taxonomy(session)
    if validate_species and not preloaded_taxonomy:
        sys.exit("ERROR: Taxonomy table empty, cannot validate species.\n")
    if debug:
        sys.stderr.write(
            "Taxonomy/synonym tables contains %i distinct species names\n"
            % len(preloaded_taxonomy)
        )

    if not name:
        name = "Import of %s" % os.path.basename(fasta_file)

    md5 = md5_hexdigest(fasta_file)

    # TODO - explicit check for reusing name, and/or unique in schema
    # TODO - explicit check for reusing MD5 (not just DB schema check)
    db_source = DataSource(
        name=name,
        uri=fasta_file,
        md5=md5,
        notes="Imported with thapbi_pict legacy_import v%s" % __version__,
    )
    session.add(db_source)

    seq_count = 0
    its1_seq_count = 0
    good_seq_count = 0

    entry_count = 0
    bad_entries = 0
    bad_sp_entries = 0
    good_entries = 0
    idn_set = set()
    multiple_its1 = False

    for title, seq, hmm_name, its1_seqs in filter_for_ITS1(fasta_file, cache_dir=None):
        seq_count += 1
        if not its1_seqs:
            assert not hmm_name, hmm_name
            if debug:
                sys.stderr.write(
                    "DEBUG: Ignoring entry with no HMM matches: %s\n" % title
                )
            continue
        its1_seq_count += 1

        # One sequence can have multiple entries
        idn = title.split(None, 1)[0]
        if idn in idn_set:
            sys.stderr.write("WARNING: Duplicated identifier %r\n" % idn)
        idn_set.add(idn)

        if len(its1_seqs) > 1:
            sys.stderr.write("WARNING: %i HMM matches in %s\n" % (len(its1_seqs), idn))
            multiple_its1 = True

        entries = fasta_entry_fn(title)
        if not entries:
            sys.stderr.write(
                "WARNING: Ignoring %r\n"
                % (title if len(title) < 70 else title[:66] + "...")
            )
            continue

        accepted = False
        for entry in entries:
            entry_count += 1
            try:
                taxid, name = entry_taxonomy_fn(entry, preloaded_taxonomy)
            except ValueError:
                bad_entries += 1
                sys.stderr.write(
                    "WARNING: Ignoring %r\n"
                    % (entry if len(entry) < 60 else entry[:67] + "...")
                )
                continue

            assert isinstance(name, str), name

            if name.lower().startswith("uncultured "):
                bad_entries += 1
                sys.stderr.write(
                    "WARNING: Ignoring %r\n"
                    % (entry if len(entry) < 60 else entry[:67] + "...")
                )
                continue

            if not taxid and not name:
                bad_entries += 1
                sys.stderr.write("WARNING: No species information: %r\n" % idn)
                continue

            # Load into the DB
            #
            # Store "Phytophthora aff infestans" as
            # genus "Phytophthora", species "aff infestans"
            #
            # Note even for genus only, must check synonyms,
            # e.g. "Pythium undulatum" -> "Phytophthora undulatum"
            if debug and not name:
                sys.stderr.write("WARNING: No species information from %r\n" % entry)

            assert not name.startswith("P."), title

            if genus_only:
                taxonomy = lookup_genus(session, name)
            else:
                taxonomy = lookup_species(session, name)
            if not taxonomy:
                if debug:
                    sys.stderr.write(
                        "WARNING: Could not validate species %r from %r\n"
                        % (name, entry)
                    )
                if validate_species:
                    bad_sp_entries += 1
                    continue
                if name in additional_taxonomy:
                    # Appeared earlier in this import
                    taxonomy = additional_taxonomy[name]
                else:
                    # Must add this now
                    genus, species = genus_species_split(name)
                    taxonomy = Taxonomy(
                        genus=genus, species="" if genus_only else species, ncbi_taxid=0
                    )
                    session.add(taxonomy)
                    additional_taxonomy[name] = taxonomy

            assert taxonomy is not None
            for its1_seq in its1_seqs:
                # its1_seq_count += 1
                its1_md5 = md5seq(its1_seq)

                # Is sequence already there? e.g. duplicate sequences in FASTA file
                its1 = (
                    session.query(ITS1)
                    .filter_by(md5=its1_md5, sequence=its1_seq)
                    .one_or_none()
                )
                if its1 is None:
                    its1 = ITS1(md5=its1_md5, sequence=its1_seq)
                    session.add(its1)
                record_entry = SequenceSource(
                    source_accession=entry.split(None, 1)[0],
                    source=db_source,
                    its1=its1,
                    sequence=seq,
                    original_taxonomy=taxonomy,
                    current_taxonomy=taxonomy,
                    seq_strategy=0,
                    seq_platform=0,
                    curated_trust=0,
                )
                session.add(record_entry)
            good_entries += 1  # count once?
            accepted = True
        if accepted:
            good_seq_count += 1

    session.commit()
    sys.stderr.write(
        "File %s had %i sequences. Found %i with ITS1, of which %i accepted.\n"
        % (fasta_file, seq_count, its1_seq_count, good_seq_count)
    )
    assert multiple_its1 or its1_seq_count <= seq_count, (its1_seq_count, seq_count)
    assert bad_entries <= entry_count, (bad_entries, entry_count)
    assert multiple_its1 or good_entries <= entry_count, (good_entries, entry_count)
    if validate_species:
        sys.stderr.write(
            "Of %i potential entries, %i unparsable, %i failed sp. validation, %i OK.\n"
            % (entry_count, bad_entries, bad_sp_entries, good_entries)
        )
        assert entry_count == good_entries + bad_entries + bad_sp_entries
    else:
        sys.stderr.write(
            "Of %i potential entries, loaded %i entries, %i failed parsing.\n"
            % (entry_count, good_entries, bad_entries)
        )
        assert bad_sp_entries == 0, bad_sp_entries
        assert entry_count == good_entries + bad_entries
