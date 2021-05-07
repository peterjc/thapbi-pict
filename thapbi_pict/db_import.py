# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
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

from Bio.SeqIO.FastaIO import SimpleFastaParser

from . import __version__
from .db_orm import connect_to_db
from .db_orm import DataSource
from .db_orm import RefMarker
from .db_orm import SequenceSource
from .db_orm import Synonym
from .db_orm import Taxonomy
from .utils import find_requested_files
from .utils import genus_species_name
from .utils import genus_species_split
from .utils import md5_hexdigest
from .utils import md5seq


def parse_ncbi_fasta_entry(text, known_species=None):
    """Split an entry of Accession Genus Species-name Description.

    Returns a two-tuple: taxid (always zero), presumed genus-species
    (taken as two words by default if cannot be matched to a provided
    known species) which may be the empty string.

    >>> parse_ncbi_fasta_entry('LC159493.1 Phytophthora drechsleri genes ...')
    (0, 'Phytophthora drechsleri')
    >>> parse_ncbi_fasta_entry('A57915.1 Sequence 20 from Patent EP0751227')
    (0, '')
    >>> parse_ncbi_fasta_entry('Y08654.1 P.cambivora ribosomal internal ...')
    (0, '')

    Dividing the species name into genus, species, strain etc
    is not handled here.
    """  # noqa: E501
    parts = text.rstrip().split()
    taxid = 0
    name = parts[1:]  # ignore accession

    if known_species:
        while name and " ".join(name) not in known_species:
            name.pop()  # discard last word
        if len(name) > 1:
            # Found a perfect match
            return 0, " ".join(name)

    # Heuristics
    name = parts[1:3]  # assumes "Genus species" only (2 words)
    rest = parts[3:]
    assert name, text
    if len(name[0]) > 2 and name[0][0].isupper() and name[0][1] == ".":
        # Would need more information to infer the genus here
        # e.g. Y08654.1 P.cambivora ribosomal internal transcribed spacer, ITS1
        return 0, ""
    while rest and name[-1] in ("taxon", "aff.", "cf.", "x"):
        # Looks like species name needs at least one more word...
        # Note that sp. or sp doesn't always have another word.
        name.append(rest.pop(0))
    if name[0] == "Sequence":
        # Another special case
        # e.g. A57915.1 Sequence 20 from Patent EP0751227
        return 0, ""
    if len(rest) >= 2 and rest[0] == "x":
        # Hybrid
        if name[0] == rest[1] and len(rest) >= 3:
            # Genus repeated
            name.append("x")
            name.append(rest[2])
        else:
            name.extend(rest[:2])
    return taxid, " ".join(name)


assert parse_ncbi_fasta_entry("LC159493.1 Phytophthora drechsleri genes ...") == (
    0,
    "Phytophthora drechsleri",
)

assert parse_ncbi_fasta_entry("A57915.1 Sequence 20 from Patent EP0751227") == (0, "")

assert parse_ncbi_fasta_entry("Y08654.1 P.cambivora ribosomal internal ...") == (0, "")

assert parse_ncbi_fasta_entry(
    "MG707849.1 Phytophthora humicola x Phytophthora inundata isolate SCVWD597 internal transcribed spacer 1, ..."  # noqa: E501
) == (0, "Phytophthora humicola x inundata")

assert parse_ncbi_fasta_entry(
    "MG707849.1 Phytophthora humicola x inundata isolate SCVWD597 internal transcribed spacer 1, ..."  # noqa: E501
) == (0, "Phytophthora humicola x inundata")


def parse_curated_fasta_entry(text, known_species=None):
    """Split an entry of "Acession genus species etc" into fields.

    Returns a two-tuple of taxid (always zero), genus-species.

    >>> parse_curated_fasta_entry('HQ013219 Phytophthora arenaria')
    (0, 'Phytophthora arenaria')

    >>> parse_curated_fasta_entry('P13660 Phytophthora aff infestans')
    (0, 'Phytophthora aff infestans')
    """
    acc, sp = text.split(None, 1)
    taxid = 0
    # if sp not in known_species:
    #     sys.stderr.write(f"WARNING: Unexpected species name {sp}\n")
    while "  " in sp:
        sp = sp.replace("  ", " ")
    return taxid, sp.strip()


assert parse_curated_fasta_entry("HQ013219 Phytophthora arenaria") == (
    0,
    "Phytophthora arenaria",
)
assert parse_curated_fasta_entry("P13660 Phytophthora aff infestans") == (
    0,
    "Phytophthora aff infestans",
)


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
    fasta_entry_fn,
    entry_taxonomy_fn,
    min_length=0,
    max_length=sys.maxsize,
    name=None,
    debug=True,
    validate_species=False,
    genus_only=False,
    tmp_dir=None,
):
    """Import a FASTA file into the database."""
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

    with open(fasta_file) as handle:
        for title, seq in SimpleFastaParser(handle):
            if "-" in seq:
                sys.exit(f"ERROR: Gap in sequence for {title}")
            # NOTE: at this point don't have full sequence prior to primer removal
            seq = seq.upper()
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

                if True:
                    # marker_seq_count += 1
                    marker_md5 = md5seq(seq)

                    # Is sequence already there? e.g. duplicate sequences in FASTA file
                    marker = (
                        session.query(RefMarker)
                        .filter_by(md5=marker_md5, sequence=seq)
                        .one_or_none()
                    )
                    if marker is None:
                        marker = RefMarker(md5=marker_md5, sequence=seq)
                        session.add(marker)
                    record_entry = SequenceSource(
                        source_accession=entry.split(None, 1)[0],
                        source=db_source,
                        marker=marker,
                        taxonomy=taxonomy,
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


def main(
    fasta,
    db_url,
    min_length=0,
    max_length=sys.maxsize,
    name=None,
    ncbi_heuristics=False,
    sep=None,
    validate_species=False,
    genus_only=False,
    ignore_prefixes=None,
    tmp_dir=None,
    debug=False,
):
    """Import FASTA file(s) into the database.

    For NCBI files, recommend using ``--ncbi`` and leave ``--sep ''``
    as it is since single entries are expected.

    For curated FASTA files, omit ``--ncbi`` and use ``--sep ';'``
    or whatever multi-entry separator you are using.

    Primer trimming is done if and only if primers are given.
    """
    if sep:

        def fasta_entry_fn(text):
            """Split FASTA entries on the separator character."""
            return [_.strip() for _ in text.split(sep)]

    else:
        if debug:
            sys.stderr.write("DEBUG: Treating each FASTA entry as a singleton.\n")

        def fasta_entry_fn(text):
            """Treat all FASTA entries as singletons."""
            return [text]

    fasta_files = find_requested_files(fasta, ".fasta", ignore_prefixes, debug=debug)
    if debug:
        sys.stderr.write(f"Classifying {len(fasta_files)} input FASTA files\n")

    for fasta_file in fasta_files:
        import_fasta_file(
            fasta_file,
            db_url,
            fasta_entry_fn,
            parse_ncbi_fasta_entry if ncbi_heuristics else parse_curated_fasta_entry,
            min_length=min_length,
            max_length=max_length,
            name=name,
            validate_species=validate_species,
            genus_only=genus_only,
            debug=debug,
        )
