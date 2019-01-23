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


def lookup_taxonomy(session, clade, genus, species):
    """Find this entry in the taxonomy table (if present)."""
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


def find_taxonomy(session, clade, sp_name, sp_name_etc, validate_species):
    """Fuzzy search for this entry in the taxonomy table (if present)."""
    assert isinstance(clade, str), clade
    assert isinstance(sp_name, str), sp_name
    assert isinstance(sp_name_etc, str), sp_name_etc
    # First loop removes words (until just one word for species)
    # in the hope of finding a match. This is for our legacy files.
    # Second loop adds words (from the species_etc string) in the
    # hope of finding a match. This is for the NCBI FASTA files.
    genus, species = sp_name.split(" ", 1) if sp_name else ("", "")
    while True:
        # sys.stderr.write("DEBUG: Trying genus=%r, species=%r\n"
        #                  % (genus, species))
        taxonomy = lookup_taxonomy(session, clade, genus, species)
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
        taxonomy = lookup_taxonomy(session, clade, genus, species)
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
    entry_count = 0
    bad_entries = 0
    bad_sp_entries = 0
    good_entries = 0
    idn_set = set()

    for title, seq, its1_seq in filter_for_ITS1(fasta_file):
        seq_count += 1
        if title.startswith("Control_"):
            if debug:
                sys.stderr.write("DEBUG: Ignoring control entry: %s\n" % title)
            continue
        if not its1_seq:
            if debug:
                sys.stderr.write("DEBUG: Ignoring non-ITS entry: %s\n" % title)
            continue

        its1_seq_count += 1

        its1_md5 = hashlib.md5(its1_seq.upper().encode("ascii")).hexdigest()

        # Is is already there? e.g. duplicate sequences in FASTA file
        its1 = (
            session.query(ITS1).filter_by(md5=its1_md5, sequence=its1_seq).one_or_none()
        )
        if its1 is None:
            its1 = ITS1(md5=its1_md5, sequence=its1_seq)
            session.add(its1)

        # One sequence can have multiple entries
        idn = title.split(None, 1)[0]
        if idn in idn_set:
            sys.stderr.write("WARNING: Duplicated identifier %r\n" % idn)
        idn_set.add(idn)

        entries = fasta_entry_fn(title)
        for entry in entries:
            entry_count += 1
            try:
                clade, name, name_etc = entry_taxonomy_fn(entry)
            except ValueError as e:
                bad_entries += 1
                sys.stderr.write("WARNING: %s - Can't parse: %r\n" % (e, idn))
                continue
            assert isinstance(name, str), name
            assert isinstance(name_etc, str), name_etc
            # Load into the DB
            # Store "Phytophthora aff infestans" as
            # genus "Phytophthora", species "aff infestans"
            if debug and not name:
                sys.stderr.write("WARNING: No species information from %r\n" % entry)
            genus, species = name.split(None, 1) if name else ("", "")
            assert genus != "P.", title
            # Is is already there? e.g. duplicate sequences in FASTA file
            # Note even if have no species text, still do the DB lookup!
            taxonomy = find_taxonomy(session, clade, name, name_etc, validate_species)
            if taxonomy is None:
                if validate_species:
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
                    # Do NOT write it to the DB
                    continue
                taxonomy = Taxonomy(
                    clade=clade, genus=genus, species=species, ncbi_taxid=0
                )
                session.add(taxonomy)

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
            good_entries += 1
            # print(clade, species, acc)
    session.commit()
    sys.stderr.write(
        "%i sequences, %i of which have ITS1, giving %i potential entries.\n"
        % (seq_count, its1_seq_count, entry_count)
    )
    assert its1_seq_count <= seq_count, (its1_seq_count, seq_count)
    assert its1_seq_count <= entry_count, (its1_seq_count, entry_count)
    assert bad_entries <= entry_count, (bad_entries, entry_count)
    assert good_entries <= entry_count, (good_entries, entry_count)
    if validate_species:
        sys.stderr.write(
            "Loaded %i entries, %i failed species validation, %i rejected.\n"
            % (good_entries, bad_sp_entries, bad_entries)
        )
        assert entry_count == good_entries + bad_entries + bad_sp_entries
    else:
        sys.stderr.write(
            "Loaded %i entries, %i rejected.\n" % (good_entries, bad_entries)
        )
        assert bad_sp_entries == 0
        assert entry_count == good_entries + bad_entries
