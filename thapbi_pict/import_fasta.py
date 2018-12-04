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


def import_fasta_file(fasta_file, db_url, name=None, debug=True,
                      fasta_split_fn=None, fasta_parse_fn=None):
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
            taxid = 0
            # Is is already there? e.g. duplicate sequences in FASTA file
            taxonomy = session.query(Taxonomy).filter_by(
                clade=clade, genus=genus, species=species,
                ncbi_taxid=taxid).one_or_none()
            if taxonomy is None:
                taxonomy = Taxonomy(
                    clade=clade, genus=genus, species=species,
                    ncbi_taxid=taxid)
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
