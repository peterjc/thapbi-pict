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


def find_taxonomy(session, clade, genus, species, validate_species):
    """Find this entry in the taxonomy table (if present)."""
    while True:
        # Can we find a match without knowing the taxid?
        taxonomy = session.query(Taxonomy).filter_by(
            clade=clade, genus=genus, species=species).one_or_none()
        if taxonomy is not None:
            # There was a unique entry already, use it.
            # It may even have an NCBI taxid?
            return taxonomy

        # Can we find a match with taxid=0?
        taxonomy = session.query(Taxonomy).filter_by(
            clade=clade, genus=genus, species=species,
            ncbi_taxid=0).one_or_none()
        if taxonomy is not None:
            # There was a unique entry already, use it.
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
            # Nope, give up
            return None
        # Remove last word, loop and try again
        species = " ".join(words[:-1])


def import_fasta_file(fasta_file, db_url, name=None, debug=True,
                      fasta_split_fn=None, fasta_parse_fn=None,
                      validate_species=False):
    """Import a FASTA file into the database.

    Optional argument fasta_split_fn is given the full FASTA
    title line, and should return a list of sub-entries (by
    default, a single entry list).
    """
    # Argument validation,
    if fasta_split_fn is None:
        def fasta_split_fn(text):
            """Treat all FASTA entries as singletons.

            Default is not to support merged FASTA entries.
            """
            return [text]

    if fasta_parse_fn is None:
        raise ValueError("Need function to split FASTA title into fields.")

    # Connect to the DB,
    Session = connect_to_db(db_url, echo=debug)
    session = Session()

    if not name:
        name = "Import of %s" % os.path.basename(fasta_file)

    md5 = md5_hexdigest(fasta_file)

    # TODO - explicit check for reusing name, and/or unique in schema
    # TODO - explicit check for reusing MD5 (not just DB schema check)
    db_source = DataSource(
        name=name,
        uri=fasta_file,
        md5=md5,
        notes="Imported with thapbi_pict legacy_import v%s" % __version__)
    session.add(db_source)

    seq_count = 0
    entry_count = 0
    bad_entry_count = 0
    idn_set = set()

    for title, seq, its1_seq in filter_for_ITS1(fasta_file):
        if title.startswith("Control_"):
            if debug:
                sys.stderr.write("Ignoring control entry: %s\n"
                                 % title)
            continue
        seq_count += 1

        if not its1_seq:
            if debug:
                sys.stderr.write("Ignoring non-ITS entry: %s\n"
                                 % title)
            continue

        its1_md5 = hashlib.md5(its1_seq.upper().encode("ascii")).hexdigest()

        # Is is already there? e.g. duplicate sequences in FASTA file
        its1 = session.query(ITS1).filter_by(
            md5=its1_md5, sequence=its1_seq).one_or_none()
        if its1 is None:
            its1 = ITS1(md5=its1_md5, sequence=its1_seq)
            session.add(its1)

        # One sequence can have multiple entries
        idn = title.split(None, 1)[0]
        if idn in idn_set:
            sys.stderr.write("WARNING: Duplicated identifier %r\n"
                             % idn)
        idn_set.add(idn)

        entries = fasta_split_fn(title)
        for entry in entries:
            entry_count += 1
            try:
                clade, name, acc = fasta_parse_fn(entry)
            except ValueError as e:
                bad_entry_count += 1
                sys.stderr.write("WARNING: %s - Can't parse: %r\n"
                                 % (e, idn))
                continue
            # Load into the DB
            # Store "Phytophthora aff infestans" as
            # genus "Phytophthora", species "aff infestans"
            if debug and not name:
                sys.stderr.write(
                    "WARNING: No species information from %r\n" % entry)
            genus, species = name.split(None, 1) if name else ("", "")
            assert genus != "P.", title
            # Is is already there? e.g. duplicate sequences in FASTA file
            taxonomy = find_taxonomy(session,
                                     clade, genus, species,
                                     validate_species)
            if taxonomy is None:
                if validate_species and name:
                    sys.stderr.write("WARNING: Could not validate %r from %r\n"
                                     % (name, entry))
                    # TODO - Find any unclassified genus entry, and use that
                taxonomy = Taxonomy(
                    clade=clade, genus=genus, species=species,
                    ncbi_taxid=0)
                session.add(taxonomy)

            # Note we use the original FASTA identifier for traceablity
            # but means the multi-entries get the same source accession
            record_entry = SequenceSource(source_accession=idn,
                                          source=db_source,
                                          its1=its1,
                                          sequence=seq,
                                          original_taxonomy=taxonomy,
                                          current_taxonomy=taxonomy,
                                          seq_strategy=0,
                                          seq_platform=0,
                                          curated_trust=0)
            session.add(record_entry)
            # print(clade, species, acc)
    session.commit()
    sys.stderr.write("%i sequences, %i entries including %i bad\n"
                     % (seq_count, entry_count, bad_entry_count))
