# Copyright 2018-2020 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Shared code for THAPBI PICT to import FASTA into our database.

This code is used both for importing NCBI formatted FASTA files, and also
importing our curated ITS1 sequence FASTA file databases - see ``ncbi.py``
and ``curated.py`` which contain specific meta-data handling code for the
different naming conventions.
"""
import os
import sys
import tempfile

from Bio.Seq import reverse_complement
from Bio.SeqIO.FastaIO import SimpleFastaParser

from . import __version__
from .db_orm import connect_to_db
from .db_orm import DataSource
from .db_orm import ITS1
from .db_orm import SequenceSource
from .db_orm import Synonym
from .db_orm import Taxonomy
from .utils import genus_species_name
from .utils import genus_species_split
from .utils import md5_hexdigest
from .utils import md5seq
from .utils import run
from .versions import check_tools


def run_cutadapt_keep(
    long_in, trimmed_out, left_primer, right_primer, debug=False, cpu=0
):
    """Run cutadapt on a single file, cropping at either primer if found,.

    The input and/or output files may be compressed as long as they
    have an appropriate suffix (e.g. gzipped with ``.gz`` suffix).
    """
    cmd = ["cutadapt"]
    if cpu:
        # Not compatible with --untrimmed-output
        cmd += ["-j", str(cpu)]
    # Can't do together as while RIGHT may match, LEFT might not
    if left_primer:
        cmd += ["-g", left_primer]
    if right_primer:
        cmd += ["-a", reverse_complement(right_primer)]
    cmd += ["-o", trimmed_out, long_in]
    return run(cmd, debug=debug)


def load_taxonomy(session):
    """Pre-load all the species and synonym names as a set."""
    names = set()
    view = session.query(Taxonomy).distinct(Taxonomy.genus, Taxonomy.species)
    for taxonomy in view:
        names.add(genus_species_name(taxonomy.genus, taxonomy.species))
    for synonym in session.query(Synonym):
        if synonym.name in names:
            sys.stderr.write(f"WARNING: Synonym {synonym.name} duplicated?\n")
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
    min_length=0,
    max_length=sys.maxsize,
    name=None,
    debug=True,
    fasta_entry_fn=None,
    entry_taxonomy_fn=None,
    validate_species=False,
    genus_only=False,
    left_primer=None,
    right_primer=None,
    tmp_dir=None,
):
    """Import a FASTA file into the database.

    For ``thapbi_pict curated-import`` some FASTA sequences are
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

    In ``thapbi_pict curated-import`` and ``thapbi_pict ncbi-import``
    the species metadata is recorded directly in the FASTA title
    lines. However, the metadata for ``thapbi_pict seq-import``
    comes from a sister TSV file, and is cross-referenced by the
    FASTA sequence identifier.

    For ``thapbi_pict curated-import`` we expect pre-trimmed curated
    marker sequences. For ``thapbi_pict seq-import`` the reads have
    been primer-trimed by ``thapbi_pict prepare-reads``. However,
    for ``thapbi_pict ncbi-import`` many of the sequences will be
    longer than the marker region of interest - and thus should be
    in-silico primer trimmed.
    """
    if left_primer or right_primer:
        check_tools(["cutadapt"], debug)

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
            sys.stderr.write(f"Ignoring empty FASTA file {fasta_file}\n")
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
            "Taxonomy/synonym tables contains"
            f" {len(preloaded_taxonomy)} distinct species names\n"
        )

    md5 = md5_hexdigest(fasta_file)

    if not name:
        name = "Import of " + os.path.basename(fasta_file)

    if left_primer or right_primer:
        if tmp_dir:
            # Up to the user to remove the files
            tmp_obj = None
            shared_tmp = tmp_dir
        else:
            tmp_obj = tempfile.TemporaryDirectory()
            shared_tmp = tmp_obj.name

        if left_primer:
            trim_left = os.path.join(shared_tmp, "cutadapt_left.fasta")
            run_cutadapt_keep(fasta_file, trim_left, left_primer, None, debug)  # cpu
        else:
            trim_left = fasta_file
        if right_primer:
            trimmed_fasta = os.path.join(shared_tmp, "cutadapt.fasta")
            run_cutadapt_keep(
                trim_left, trimmed_fasta, None, right_primer, debug
            )  # cpu
        else:
            trimmed_fasta = trim_left
    else:
        trimmed_fasta = fasta_file

    # TODO - explicit check for reusing name, and/or unique in schema
    # TODO - explicit check for reusing MD5 (not just DB schema check)
    db_source = DataSource(
        name=name,
        uri=fasta_file,
        md5=md5,
        notes=f"Imported with thapbi_pict v{__version__}",
    )
    session.add(db_source)

    seq_count = 0
    good_seq_count = 0

    bad_species = set()

    entry_count = 0
    bad_entries = 0
    bad_sp_entries = 0
    good_entries = 0
    idn_set = set()

    valid_letters = set("GATCRYWSMKHBVDN")

    with open(trimmed_fasta) as handle:
        for title, seq in SimpleFastaParser(handle):
            if "-" in seq:
                sys.exit(f"ERROR: Gap in sequence for {title}")
            if set(seq.upper()).difference(valid_letters):
                bad = ", ".join(sorted(set(seq.upper()).difference(valid_letters)))
                sys.exit(
                    f"ERROR: Non-IUPAC DNA character(s) {bad} in sequence for {title}"
                )

            seq_count += 1
            idn = title.split(None, 1)[0]

            if not (min_length <= len(seq) <= max_length):
                if debug:
                    sys.stderr.write(f"DEBUG: Rejected {idn} as length {len(seq)}\n")
                continue

            # One sequence can have multiple entries
            if idn in idn_set:
                sys.stderr.write(f"WARNING: Duplicated identifier {idn}\n")
            idn_set.add(idn)

            entries = fasta_entry_fn(title)
            if not entries:
                sys.stderr.write(
                    "WARNING: Based on name, ignoring %r\n"
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
                        "WARNING: Could not parse entry %r\n"
                        % (entry if len(entry) < 60 else entry[:67] + "...")
                    )
                    continue

                assert isinstance(name, str), name

                if name.lower().startswith("uncultured "):
                    bad_sp_entries += 1
                    sys.stderr.write(
                        "WARNING: Uncultured, so ignoring %r\n"
                        % (entry if len(entry) < 60 else entry[:67] + "...")
                    )
                    continue

                if not taxid and not name:
                    bad_sp_entries += 1
                    sys.stderr.write(f"WARNING: No species information: {idn!r}\n")
                    continue

                # Load into the DB
                #
                # Store "Phytophthora aff infestans" as
                # genus "Phytophthora", species "aff infestans"
                #
                # Note even for genus only, must check synonyms,
                # e.g. "Pythium undulatum" -> "Phytophthora undulatum"
                if debug and not name:
                    sys.stderr.write(
                        f"WARNING: No species information from {entry!r}\n"
                    )

                assert not name.startswith("P."), title

                if genus_only:
                    taxonomy = lookup_genus(session, name)
                else:
                    taxonomy = lookup_species(session, name)
                if not taxonomy:
                    if preloaded_taxonomy and debug and name not in bad_species:
                        sys.stderr.write(
                            "WARNING: Could not validate species"
                            f" {name!r} from {entry!r}\n"
                        )
                    bad_species.add(name)  # To avoid repeat warnings
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
                            genus=genus,
                            species="" if genus_only else species,
                            ncbi_taxid=0,
                        )
                        session.add(taxonomy)
                        additional_taxonomy[name] = taxonomy

                assert taxonomy is not None
                # at this point we don't have the full sequence
                # prior to any primer removal...
                its1_seq = seq.upper()
                if True:
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
        f"File {fasta_file} had {seq_count} sequences, "
        f"of which {good_seq_count} accepted.\n"
    )
    assert bad_entries + bad_sp_entries <= entry_count, (
        bad_entries,
        bad_sp_entries,
        entry_count,
    )
    assert good_entries <= entry_count, (good_entries, entry_count)
    if validate_species:
        sys.stderr.write(
            f"Of {entry_count} potential entries, {bad_entries} unparsable,"
            f" {bad_sp_entries} failed sp. validation, {good_entries} OK.\n"
        )
        assert entry_count == good_entries + bad_entries + bad_sp_entries
    else:
        sys.stderr.write(
            f"Of {entry_count} potential entries, loaded {good_entries} entries,"
            f" {bad_entries} failed parsing.\n"
        )
        assert (
            entry_count == good_entries + bad_entries + bad_sp_entries
        ), f"{entry_count} != {good_entries} + {bad_entries} + {bad_sp_entries}"

    if bad_species and (validate_species or debug):
        sys.stderr.write(
            f"Could not validate {len(bad_species)} different species names\n"
        )
