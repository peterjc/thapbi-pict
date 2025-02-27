# Copyright 2018-2024 by Peter Cock, The James Hutton Institute.
# Revisions copyright 2024 by Peter Cock, University of Strathclyde.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Shared code for THAPBI PICT to import FASTA into our database.

This code is used for importing NCBI formatted FASTA files, our curated ITS1
sequence FASTA file databases, and other other FASTA naming conventions.
"""

from __future__ import annotations

import os
import re
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser

from . import __version__
from .db_orm import connect_to_db
from .db_orm import DataSource
from .db_orm import MarkerDef
from .db_orm import MarkerSeq
from .db_orm import SeqSource
from .db_orm import Synonym
from .db_orm import Taxonomy
from .utils import find_requested_files
from .utils import genus_species_name
from .utils import genus_species_split
from .utils import md5_hexdigest
from .utils import md5seq
from .utils import reject_species_name
from .utils import valid_marker_name

DEF_MIN_LENGTH = 100
DEF_MAX_LENGTH = 1000

taxid_regex = re.compile(r"(ncbi|[ _:;({\[\-\t])taxid=\d+")


def parse_ncbi_fasta_entry(
    text: str, known_species: list[str] | None = None
) -> tuple[int, str]:
    """Split an entry of Accession Genus Species-name Description.

    Returns a two-tuple: taxid (always zero), presumed genus-species (may be
    the empty string).

    >>> parse_ncbi_fasta_entry("LC159493.1 Phytophthora drechsleri genes ...")
    (0, 'Phytophthora drechsleri')
    >>> parse_ncbi_fasta_entry("A57915.1 Sequence 20 from Patent EP0751227")
    (0, '')
    >>> parse_ncbi_fasta_entry("Y08654.1 P.cambivora ribosomal internal ...")
    (0, '')

    If a list of known species are used, then right most word is dropped until
    the text matches a known name. This discards any description (and strain
    level information if the list is only to species level).

    If there is no match to the provided names, heuristics are used but this
    defaults to the first two words.

    Dividing the species name into genus, species, strain etc is not handled
    here.
    """  # noqa: E501
    parts = text.rstrip().split()
    taxid = 0
    name = parts[1:]  # ignore accession

    if known_species:
        while name and " ".join(name) not in known_species:
            name.pop()  # discard last word
        if len(name) > 1:
            # Found a perfect match (even if it could be genus only)
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
# TODO: Should we stop converting "Phytophthora humicola x Phytophthora inundata"
# into "Phytophthora humicola x inundata"?


def parse_ncbi_taxid_entry(
    text: str, know_species: list[str] | None = None
) -> tuple[int, str]:
    """Find any NCBI taxid as a pattern in the text.

    Returns a two-tuple of taxid (zero if not found), and an
    empty string (use the taxonomy table in the DB to get the
    genus-species).

    Uses a regular expression based on taxid=<digits>, and
    only considers the first match:

    >>> parse_ncbi_taxid_entry("HQ013219 Phytophthora arenaria [taxid=]")
    (0, '')
    >>> parse_ncbi_taxid_entry("HQ013219 Phytophthora arenaria [taxid=123] [taxid=456]")
    (123, '')
    """
    match = taxid_regex.search(text)
    if match:
        return int(match.group().split("=", 1)[1]), ""
    return 0, ""


assert parse_ncbi_taxid_entry("HQ013219 Phytophthora arenaria") == (0, "")
assert parse_ncbi_taxid_entry(
    "HQ013219 Phytophthora arenaria [taxid=123] [taxid=456]"
) == (123, "")
assert parse_ncbi_taxid_entry(
    "HQ013219 Phytophthora arenaria [taxid=NA] [taxid=456]"
) == (456, "")
assert parse_ncbi_taxid_entry(
    "HQ013219 Phytophthora arenaria [key=value;taxid=456;key=value]"
) == (456, "")
assert parse_ncbi_taxid_entry("HQ013219:Phytophthora_arenaria:taxid=456") == (456, "")


def parse_curated_fasta_entry(
    text: str, known_species: list[str] | None = None
) -> tuple[int, str]:
    """Split an entry of "Accession genus species etc" into fields.

    Does not use the optional known_species argument.

    Returns a two-tuple of taxid (0 unless taxid=... entry found), genus-species.

    >>> parse_curated_fasta_entry("HQ013219 Phytophthora arenaria")
    (0, 'Phytophthora arenaria')

    Will look for an NCBI taxid after the species name (and ignore anything
    following that, such as other key=value entries):

    >>> parse_curated_fasta_entry("P13660 Phytophthora aff infestans taxid=907744 etc")
    (907744, 'Phytophthora aff infestans')

    In this example we expect the NCBI taxid will be matched to a pre-loaded
    species name to be used in preference (i.e. 'Phytophthora aff. infestans'
    with a dot in it).
    """
    acc, sp = text.split(None, 1)
    taxid = 0
    match = taxid_regex.search(sp)
    if match:
        taxid = int(match.group().split("=", 1)[1])
        sp = sp[: match.start()]
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
assert parse_curated_fasta_entry(
    "P13660 Phytophthora aff infestans taxid=907744 ignored text"
) == (
    907744,
    "Phytophthora aff infestans",
)
assert parse_curated_fasta_entry(
    "P13660 Phytophthora aff infestans ncbitaxid=907744 ignored text"
) == (
    907744,
    "Phytophthora aff infestans",
)
assert parse_curated_fasta_entry(
    "P13660 Phytophthora aff infestans [taxid=907744] ignored text"
) == (
    907744,
    "Phytophthora aff infestans",
)


def parse_sintax_fasta_entry(
    text: str, known_species: list[str] | None = None
) -> tuple[int, str]:
    """Extract the species from SINTAX taxonomy annotation.

    See https://drive5.com/usearch/manual/tax_annot.html which defines
    this taxonomy annotation convention as used in USEARCH and VSEARCH.
    The tax=names field is separated from other fields in the FASTA
    description line by semi-colons, for example:

    >>> entry = "X80725_S000004313;tax=d:...,g:Escherichia/Shigella,s:Escherichia_coli"
    >>> parse_sintax_fasta_entry(entry)
    (0, 'Escherichia coli')

    If there is no species entry (prefix ``s:``) then the genus is returned
    (prefix ``g:``), else the empty string:

    >>> parse_sintax_fasta_entry("AB008314;tax=d:...,g:Streptococcus;")
    (0, 'Streptococcus')

    If the species entry is missing the genus information (which may happen
    depending how the file was generated), that is inferred heuristically:

    >>> entry = "X80725_S000004313;tax=d:...,g:Escherichia,s:coli"
    >>> parse_sintax_fasta_entry(entry)
    (0, 'Escherichia coli')

    This can be unclear:

    >>> entry = ">X80725_S000004313;tax=d:...,g:Escherichia/Shigella,s:Escherichia_coli"
    >>> parse_sintax_fasta_entry(entry)
    (0, 'Escherichia coli')
    """
    genus = species = ""
    for part in text.split(";"):
        part = part.replace("_", " ").strip()
        while "  " in part:
            part = part.replace("  ", " ")
        if not part.startswith("tax="):
            continue
        fields = [_.strip() for _ in part.split(",")]
        for field in fields:
            if field.startswith("g:"):
                genus = field[2:]
            elif field.startswith("s:"):
                species = field[2:]
    if genus and species:
        if species.startswith(genus) or species.split(None, 1)[0] in genus:
            # Good, genus is present
            return 0, species
        elif " " not in species or species[0].islower():
            # Looks like the genus is missing, perform ad hoc fix
            return 0, genus + " " + species
        else:
            # Unclear... assume the genus is present as per SINTAX expectations
            sys.stderr.write(f"WARNING: Genus vs species unclear in {text}\n")
            return 0, species
    elif genus:
        return 0, genus
    elif species:
        # Missing genus?
        msg = f"Failed to parse SINTAX entry: {text!r}"
        raise ValueError(msg)
    else:
        # Neither found
        return 0, ""


assert parse_sintax_fasta_entry(
    "AB008314;tax=d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,g:Streptococcus;"  # noqa: E501
) == (0, "Streptococcus")

assert parse_sintax_fasta_entry(
    ">X80725_S000004313;tax=d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,o:Enterobacteriales,f:Enterobacteriaceae,g:Escherichia/Shigella,s:Escherichia_coli"  # noqa: E501
) == (0, "Escherichia coli")


def parse_obitools_fasta_entry(
    text: str, known_species: list[str] | None = None
) -> tuple[int, str]:
    """Parse species from the OBITools extended FASTA header.

    See https://pythonhosted.org/OBITools/attributes.html which explains that
    OBITools splits the FASTA line into identifier, zero or more key=value;
    entries, and a free text description.

    We are specifically interested in the species_name, genus_name (used if
    species_name is missing), and taxid.

    >>> entry = "AP009202 species_name=Abalistes stellaris; taxid=392897; ..."
    >>> parse_obitools_fasta_entry(entry)
    (392897, 'Abalistes stellaris')

    Note this will *not* try to parse any key=value entries embedded in the
    first word (which taken as the identifier).
    """
    taxid = 0
    sp = ""
    identifier, description = text.split(None, 1)
    for part in description.split(";"):
        part = part.strip()  # We may be more lienent that OBITools here
        if part.startswith("taxid="):
            taxid = int(part[6:].strip())
        elif part.startswith("species_name="):
            sp = part[13:].strip()
        elif not sp and part.startswith("genus_name="):
            sp = part[11:].strip()
    return taxid, sp


assert parse_obitools_fasta_entry(
    "AP009202 species_name=Abalistes stellaris; taxid=392897; genus_name=Abalistes; rank=species; Abalistes stellaris mitochondrial DNA, complete genome"  # noqa: E501
) == (392897, "Abalistes stellaris")

assert parse_obitools_fasta_entry(
    "MF101792 family_name=Acipenseridae; species_name=Scaphirhynchus suttkusi; family=7900; reverse_match=CTTCCGGTACACTTACCATG; taxid=36179; rank=species; forward_error=0; forward_tm=60.26; genus_name=Scaphirhynchus; seq_length_ori=16495; forward_match=ACACCGCCCGTCACTCT; reverse_tm=54.79; genus=7909; reverse_error=0; species=36179; strand=D; Scaphirhynchus suttkusi isolate NFWFLH10433 mitochondrion, complete genome"  # noqa: E501
) == (36179, "Scaphirhynchus suttkusi")

fasta_parsing_function = {
    "simple": parse_curated_fasta_entry,
    "ncbi": parse_ncbi_fasta_entry,
    "sintax": parse_sintax_fasta_entry,
    "taxid": parse_ncbi_taxid_entry,
    "obitools": parse_obitools_fasta_entry,
}


def load_taxonomy(session) -> set[str]:
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


def lookup_species(session, name: str):
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


def lookup_genus(session, name: str):
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
    marker,
    left_primer=None,
    right_primer=None,
    min_length=None,
    max_length=None,
    name=None,
    trim=True,
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
    session = connect_to_db(db_url, echo=False)  # echo=debug

    preloaded_taxonomy = load_taxonomy(session)
    if validate_species and not preloaded_taxonomy:
        sys.exit("ERROR: Taxonomy table empty, cannot validate species.\n")
    if debug:
        sys.stderr.write(
            "Taxonomy/synonym tables contains"
            f" {len(preloaded_taxonomy)} distinct species names\n"
        )

    if marker:
        reference_marker = (
            session.query(MarkerDef).filter(MarkerDef.name == marker).one_or_none()
        )
    else:
        # Can omit if DB has one and only one entry...
        reference_marker = session.query(MarkerDef).one_or_none()
        if reference_marker:
            if debug:
                sys.stderr.write(
                    f"DEBUG: Defaulting to only marker in the DB: {marker}\n"
                )
            marker = reference_marker.name
        else:
            # Zero or more than one in DB
            sys.exit(
                "ERROR: Marker name is required (unless there is only one in the DB)."
            )
    if reference_marker:
        # Load primers etc
        if left_primer and reference_marker.left_primer != left_primer:
            sys.exit(
                "ERROR: Given left primer "
                f"{left_primer} does not match DB {reference_marker.left_primer}"
            )
        if right_primer and reference_marker.right_primer != right_primer:
            sys.exit(
                "ERROR: Given right primer "
                f"{right_primer} does not match DB {reference_marker.right_primer}"
            )
        left_primer = reference_marker.left_primer
        right_primer = reference_marker.right_primer
        if min_length is None:
            min_length = reference_marker.min_length
        elif min_length < reference_marker.min_length:
            sys.exit(
                "ERROR: Requested minimum length "
                f"{min_length} lower than that in DB {reference_marker.min_length}"
            )
        if max_length is None:
            max_length = reference_marker.max_length
        elif reference_marker.max_length < max_length:
            sys.exit(
                "ERROR: Requested maximum length "
                f"{max_length} exceeds that in DB {reference_marker.max_length}"
            )
    elif not left_primer or not right_primer:
        sys.exit("ERROR: Both primers must be supplied when defining a new marker.")
    else:
        # New marker!
        if not valid_marker_name(marker):
            sys.exit(
                "ERROR: Inappropriate marker name. Please use only letters,"
                "and if you wish numbers and/or the minus sign."
            )
        if min_length is None:
            min_length = DEF_MIN_LENGTH
        if max_length is None:
            max_length = DEF_MAX_LENGTH
        if debug:
            sys.stderr.write(
                f"DEBUG: New marker {marker} primers {left_primer} & {right_primer}\n"
            )
        reference_marker = MarkerDef()
        reference_marker.name = marker
        reference_marker.left_primer = left_primer
        reference_marker.right_primer = right_primer
        reference_marker.min_length = min_length
        reference_marker.max_length = max_length
        session.add(reference_marker)

    md5 = md5_hexdigest(fasta_file)

    if not name:
        name = "Import of " + os.path.basename(fasta_file)

    # TODO - explicit check for reusing name, and/or unique in schema
    # TODO - explicit check for reusing MD5 (not just DB schema check)
    db_source = DataSource(
        name=name,
        uri=fasta_file,
        md5=md5,
        notes=f"{marker} imported with thapbi_pict v{__version__}",
    )
    session.add(db_source)

    seq_count = 0
    good_seq_count = 0

    bad_species = set()

    entry_count = 0
    bad_entries = 0
    bad_sp_entries = 0
    downgraded_entries = 0  # unknown species --> genus only
    good_entries = 0
    idn_set = set()

    valid_letters = set("GATCRYWSMKHBVDN")

    existing_taxonomy: dict[str, str] = {}
    existing_sequences: dict[str, str] = {}
    additional_taxonomy: dict[str, str] = {}
    additional_sequences: dict[str, str] = {}
    record_entries = []
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
                except ValueError as e:
                    bad_entries += 1
                    sys.stderr.write(f"WARNING: Could not parse entry - {e}\n")
                    continue

                assert isinstance(name, str), name

                if reject_species_name(name):
                    bad_sp_entries += 1
                    sys.stderr.write(
                        "WARNING: Ignoring %r\n"
                        % (entry if len(entry) < 60 else entry[:67] + "...")
                    )
                    continue

                if taxid:
                    # Attempt to lookup the taxid to get the species name
                    taxonomy = (
                        session.query(Taxonomy)
                        .filter_by(ncbi_taxid=taxid)
                        .one_or_none()
                    )
                    if not taxonomy:
                        # Might be in merged.dmp, try our synonym entries
                        taxonomy = (
                            session.query(Taxonomy)
                            .join(Synonym)
                            .filter(Synonym.name == f"NCBI:taxid{taxid}")
                            .one_or_none()
                        )
                    if taxonomy:
                        name = genus_species_name(taxonomy.genus, taxonomy.species)
                    elif not name:
                        sys.stderr.write(
                            f"WARNING: No species information from NCBI:taxid{taxid}\n"
                        )
                elif not name:
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
                assert "  " not in name, title

                if name in existing_taxonomy:
                    taxonomy = existing_taxonomy[name]
                elif name in additional_taxonomy:
                    # Appeared earlier in this import
                    taxonomy = additional_taxonomy[name]
                elif genus_only:
                    taxonomy = lookup_genus(session, name)
                else:
                    taxonomy = lookup_species(session, name)
                    if not taxonomy and validate_species:
                        # In validate mode when have unknown species,
                        # will still take the genus if matches.
                        taxonomy = lookup_genus(session, name.split(None, 1)[0])
                        if taxonomy:
                            # This branch is not expected to be triggered by
                            # the NCBI input (as would have already done this
                            # as part of breaking up the FASTA line)
                            downgraded_entries += 1
                            if debug and name not in bad_species:
                                sys.stderr.write(
                                    f"WARNING: Taking genus only from {name!r}\n"
                                )
                            bad_species.add(name)  # To avoid repeat warnings
                            name = name.split(None, 1)[0]
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
                    assert name not in existing_taxonomy, name
                    assert name not in additional_taxonomy, name
                    # Must add this now
                    genus, species = genus_species_split(name)
                    taxonomy = Taxonomy(
                        genus=genus,
                        species="" if genus_only else species,
                        ncbi_taxid=0,
                    )
                    additional_taxonomy[name] = taxonomy

                assert taxonomy is not None

                if seq in existing_sequences:
                    marker_seq = existing_sequences[seq]
                elif seq in additional_sequences:
                    marker_seq = additional_sequences[seq]
                else:
                    # Is sequence already there? e.g. duplicate sequences in FASTA file
                    marker_seq = (
                        session.query(MarkerSeq).filter_by(sequence=seq).one_or_none()
                    )
                    if marker_seq:
                        existing_sequences[seq] = marker_seq
                    else:
                        marker_seq = MarkerSeq(md5=md5seq(seq), sequence=seq)
                        additional_sequences[seq] = marker_seq
                record_entries.append((entry.split(None, 1)[0], marker_seq, taxonomy))
                good_entries += 1  # count once?
                accepted = True
            if accepted:
                good_seq_count += 1

    del existing_taxonomy, existing_sequences
    # First import Taxonomy and MarkerSeq, will need the new entries' IDs:
    session.bulk_save_objects(additional_taxonomy.values(), return_defaults=True)
    del additional_taxonomy
    session.bulk_save_objects(additional_sequences.values(), return_defaults=True)
    del additional_sequences
    session.flush()
    # Should now be able to access marker_seq.id and taxonomy.id properties
    record_entries = [
        {
            "source_accession": a,
            "source_id": db_source.id,
            "marker_seq_id": s.id,
            "marker_definition_id": reference_marker.id,
            "taxonomy_id": t.id,
        }
        for a, s, t in record_entries
    ]
    session.bulk_insert_mappings(SeqSource, record_entries, return_defaults=False)
    del record_entries
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
        if downgraded_entries:
            # Won't happen on the NCBI import due to how name validation done
            sys.stderr.write(
                f"(Includes {downgraded_entries} downgraded to genus only).\n"
            )
        assert entry_count == good_entries + bad_entries + bad_sp_entries
    else:
        sys.stderr.write(
            f"Of {entry_count} potential entries, loaded {good_entries} entries,"
            f" {bad_entries} failed parsing.\n"
        )
        assert entry_count == good_entries + bad_entries + bad_sp_entries, (
            f"{entry_count} != {good_entries} + {bad_entries} + {bad_sp_entries}"
        )
        assert downgraded_entries == 0, downgraded_entries

    if bad_species and (validate_species or debug):
        sys.stderr.write(
            f"Could not validate {len(bad_species)} different species names\n"
        )


def main(
    fasta,
    db_url,
    marker,
    left_primer=None,
    right_primer=None,
    min_length=0,
    max_length=sys.maxsize,
    name=None,
    convention="simple",
    sep=None,
    validate_species=False,
    genus_only=False,
    ignore_prefixes=None,
    tmp_dir=None,
    debug=False,
):
    r"""Import FASTA file(s) into the database.

    For curated FASTA files, use convention "simple" (default here and at the
    command line), and specify any multi-entry separator you are using.

    For NCBI files, convention "ncbi" and for the separator use Ctrl+A (type
    ``-s $'\001'`` at the command line) if appropriate, or "" or None
    (function default) if single entries are expected.
    """
    if sep:
        if convention in ["sintax", "obitools"]:
            sys.exit(f"ERROR: Can't use separator with {convention} naming.")

        if debug:
            sys.stderr.write(f"DEBUG: Splitting each FASTA entry using {sep!r}.\n")

        def fasta_entry_fn(text):
            """Split FASTA entries on the separator character."""
            return [_.strip() for _ in text.split(sep)]

    else:
        if debug:
            sys.stderr.write("DEBUG: Treating each FASTA entry as a singleton.\n")

        def fasta_entry_fn(text):
            """Treat all FASTA entries as singletons."""
            return [text]

    fasta_files = find_requested_files(
        fasta, (".fasta", ".fa"), ignore_prefixes, debug=debug
    )
    if debug:
        sys.stderr.write(
            f"Classifying {len(fasta_files)} input {convention} FASTA files\n"
        )

    for fasta_file in fasta_files:
        import_fasta_file(
            fasta_file,
            db_url,
            fasta_entry_fn,
            fasta_parsing_function[convention],
            marker,
            left_primer,
            right_primer,
            min_length=min_length,
            max_length=max_length,
            name=name,
            validate_species=validate_species,
            genus_only=genus_only,
            debug=debug,
        )
