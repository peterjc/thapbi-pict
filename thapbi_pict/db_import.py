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

import hashlib
import os
import sys

from . import __version__
from .db_orm import DataSource, ITS1, SequenceSource
from .db_orm import Taxonomy
from .db_orm import connect_to_db
from .hmm import filter_for_ITS1
from .utils import genus_species_name
from .utils import md5seq
from .versions import check_tools


def md5_hexdigest(filename, chunk_size=1024):
    """Return the MD5 hex-digest of the given file."""
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                # EOF
                break
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def lookup_genus_taxonomy(session, genus):
    """Find this genus in the taxonomy table (if present)."""
    assert isinstance(genus, str), genus
    # Can we find a match without knowing the taxid?
    taxonomy = session.query(Taxonomy).filter_by(genus=genus, species="").one_or_none()
    if taxonomy is not None:
        # There was a unique entry already, use it.
        # It may even have an NCBI taxid?
        return taxonomy
    # Can we find a match with taxid=0?
    taxonomy = session.query(Taxonomy).filter_by(genus=genus, species="").one_or_none()
    if taxonomy is not None:
        # There was a unique entry already, use it.
        return taxonomy


def lookup_species_taxonomy(session, clade, genus, species):
    """Find this species entry in the taxonomy table (if present)."""
    assert isinstance(clade, str), clade
    assert isinstance(genus, str), genus
    assert isinstance(species, str), species
    # Can we find a match without knowing the taxid?
    taxonomy = (
        session.query(Taxonomy)
        .filter_by(clade=clade, genus=genus, species=species)
        .one_or_none()
    )
    if taxonomy is not None:
        # There was a unique entry already, use it.
        # It may even have an NCBI taxid?
        return taxonomy

    # Can we find a match with taxid=0?
    taxonomy = (
        session.query(Taxonomy)
        .filter_by(clade=clade, genus=genus, species=species, ncbi_taxid=0)
        .one_or_none()
    )
    if taxonomy is not None:
        # There was a unique entry already, use it.
        return taxonomy

    # Can we find a match without the clade?
    if clade:
        # If there is a unique match without the clade,
        # use it to get the NCBI taxid:
        taxonomy = (
            session.query(Taxonomy)
            .filter_by(clade="", genus=genus, species=species)
            .one_or_none()
        )
        if taxonomy is not None and taxonomy.ncbi_taxid:
            # There was a unique entry already, use as template!
            taxonomy = Taxonomy(
                clade=clade,
                genus=genus,
                species=species,
                ncbi_taxid=taxonomy.ncbi_taxid,
            )
            session.add(taxonomy)  # Can we refactor this?
            return taxonomy


def find_taxonomy(session, taxid, clade, sp_name, sp_name_etc, validate_species):
    """Fuzzy search for this entry in the taxonomy table (if present)."""
    assert isinstance(clade, str), clade
    assert isinstance(sp_name, str), sp_name
    assert isinstance(sp_name_etc, str), sp_name_etc

    if taxid:
        # Perfect match?
        genus, species = sp_name.split(" ", 1) if sp_name else ("", "")
        taxonomy = (
            session.query(Taxonomy)
            .filter_by(ncbi_taxid=taxid, clade=clade, genus=genus, species=species)
            .one_or_none()
        )
        if taxonomy is not None:
            return taxonomy
        # Ignoring clade?
        taxonomy = (
            session.query(Taxonomy)
            .filter_by(ncbi_taxid=taxid, genus=genus, species=species)
            .one_or_none()
        )
        if taxonomy is not None:
            return taxonomy
        # Ignoring species name
        taxonomy = session.query(Taxonomy).filter_by(ncbi_taxid=taxid).one_or_none()
        if taxonomy is not None and taxonomy.species:
            sys.stderr.write(
                "WARNING: Using taxid %i, mapped %r to %s"
                % (taxid, sp_name, genus_species_name(taxonomy.genus, taxonomy.species))
            )
            return taxonomy
        sys.exit("WARNING: Could not uniquely match taxid %i\n" % taxid)

    # First loop removes words (until just one word for species)
    # in the hope of finding a match. This is for our legacy files.
    # Second loop adds words (from the species_etc string) in the
    # hope of finding a match. This is for the NCBI FASTA files.
    genus, species = sp_name.split(" ", 1) if sp_name else ("", "")
    while True:
        # sys.stderr.write("DEBUG: Trying genus=%r, species=%r\n"
        #                  % (genus, species))
        taxonomy = lookup_species_taxonomy(session, clade, genus, species)
        if taxonomy:
            return taxonomy

        if not validate_species:
            # If the DB has not been preloaded, we don't
            # want to try trimming the species - it would
            # make the import results order dependent etc.
            return None

        # No match! The DB has been preloaded, so can we trim
        # unwanted text off species name to get a match?
        if not species:
            return None
        words = species.split()
        if len(words) == 1:
            # Nope, give up the trimming approach
            break
        # Remove last word, loop and try again
        species = " ".join(words[:-1])

    # Removing words failed, can we try adding words?
    if not sp_name or not sp_name_etc:
        return None
    extra = sp_name_etc.split(None, 6)
    for i in range(1, max(5, len(extra))):
        genus, species = sp_name.split(" ", 1)
        species = " ".join([species] + extra[:i])
        # sys.stderr.write(
        #     "DEBUG: Extending species name, trying genus=%r, species=%r\n"
        #     % (genus, species))
        taxonomy = lookup_species_taxonomy(session, clade, genus, species)
        if taxonomy:
            return taxonomy

    # Failed to find a match
    return None


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

        def fasta_entry_fn(text):
            """Treat all FASTA entries as singletons.

            Default is not to support merged FASTA entries, and accept
            all the sequences.
            """
            return [text]

    if entry_taxonomy_fn is None:
        raise ValueError("Need function to get meta-data from FASTA title.")

    if os.stat(fasta_file).st_size == 0:
        if debug:
            sys.stderr.write("Ignoring empty FASTA file %s\n" % fasta_file)
        return

    # Connect to the DB,
    Session = connect_to_db(db_url, echo=debug)
    session = Session()

    if validate_species:
        count = (
            session.query(Taxonomy).distinct(Taxonomy.genus, Taxonomy.species).count()
        )
        if debug:
            sys.stderr.write("Taxonomy table contains %i distinct species\n" % count)
        if not count:
            sys.exit("ERROR: Taxonomy table empty, cannot validate species.\n")

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

        # One sequence can have multiple entries
        idn = title.split(None, 1)[0]
        if idn in idn_set:
            sys.stderr.write("WARNING: Duplicated identifier %r\n" % idn)
        idn_set.add(idn)

        entries = fasta_entry_fn(title)
        accepted_entries = []  # Some or all may fail species validation
        for entry in entries:
            entry_count += 1
            try:
                taxid, clade, name, name_etc = entry_taxonomy_fn(entry)
            except ValueError as e:
                bad_entries += 1
                sys.stderr.write("WARNING: %s - Can't parse: %r\n" % (e, idn))
                continue
            assert isinstance(name, str), name
            assert isinstance(name_etc, str), name_etc
            if not taxid and not clade and not name:
                bad_entries += 1
                sys.stderr.write("WARNING: No species information: %r\n" % idn)
                continue
            # Load into the DB
            # Store "Phytophthora aff infestans" as
            # genus "Phytophthora", species "aff infestans"
            if debug and not name:
                sys.stderr.write("WARNING: No species information from %r\n" % entry)
            genus, species = name.split(None, 1) if name else ("", "")
            assert genus != "P.", title
            if genus_only:
                # Genus is assumes to be one word, no fuzzy matching needed
                # TODO: Take any taxid and walk up tree to genus?
                # (useful in corner case of a genus being renamed)
                taxonomy = lookup_genus_taxonomy(session, genus)
            else:
                # Time for some fuzzy matching...
                taxonomy = find_taxonomy(
                    session, taxid, clade, name, name_etc, validate_species
                )
            if taxonomy is None and validate_species:
                bad_sp_entries += 1
                if name and debug:
                    sys.stderr.write(
                        "WARNING: Could not validate species %r from %r\n"
                        % (name, entry)
                    )
                if not name:
                    sys.stderr.write(
                        "WARNING: Could not determine species from %r\n" % entry
                    )
                # Do NOT write it to the DB!
            else:
                accepted_entries.append((entry, clade, genus, species, taxid, taxonomy))

        if not accepted_entries:
            continue

        if len(its1_seqs) > 1:
            sys.stderr.write("WARNING: %i HMM matches in %s\n" % (len(its1_seqs), idn))
            multiple_its1 = True

        for its1_seq in its1_seqs:
            its1_seq_count += 1
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
            good_seq_count += 1

            for entry, clade, genus, species, taxid, taxonomy in accepted_entries:
                if taxonomy is None:
                    assert not validate_species
                    # A previous entry in this match could have had same sp,
                    assert not taxid
                    if genus_only:
                        taxonomy = lookup_genus_taxonomy(session, genus)
                        if taxonomy is None:
                            clade = ""
                            species = ""
                            taxid = None
                            taxonomy = Taxonomy(
                                clade=clade,
                                genus=genus,
                                species=species,
                                ncbi_taxid=taxid,
                            )
                            session.add(taxonomy)
                    else:
                        taxonomy = lookup_species_taxonomy(
                            session, clade, genus, species
                        )
                        if taxonomy is None:
                            taxonomy = Taxonomy(
                                clade=clade,
                                genus=genus,
                                species=species,
                                ncbi_taxid=taxid,
                            )
                            session.add(taxonomy)

                assert taxonomy is not None
                # Note we use the original FASTA identifier for traceablity
                # but means the multi-entries get the same source accession
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
                # print(clade, species, acc)
        good_entries += len(accepted_entries)  # single counting, misleading?

    session.commit()
    sys.stderr.write(
        "File %s had %i sequences. Found %i ITS1, of which %i accepted.\n"
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
        assert bad_sp_entries == 0
        assert entry_count == good_entries + bad_entries
