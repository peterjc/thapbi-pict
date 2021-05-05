# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
r"""Code for THAPBI PICT to deal with importing NCBI FASTA files.

For example, you might perform this search via the NCBI Entrez website,
against the nucleotide database, then send it to file in FASTA format::

   its1 AND Phytophthora[Organism] AND 150:800[Sequence Length]

Or, using the NCBI Entrez Direct command line tools::

    $ esearch -db nucleotide \
        -query "its1 AND Phytophthora[Organism] AND 150:800[Sequence Length]" \
        | efetch -format fasta > example.fasta

Then, import this into your ITS1 DB using::

    $ thapbi_pict ncbi-import -d example.sqlite -i example.fasta

"""
import sys

from .db_import import import_fasta_file


def parse_fasta_entry(text, known_species=None):
    """Split an entry of Accession_Genus_Species_name_Description.

    Returns a two-tuple: taxid (always zero), presumed genus-species
    (taken as two words by default if cannot be matched to a provided
    known species) which may be the empty string.

    >>> parse_fasta_entry('LC159493.1 Phytophthora drechsleri genes ...')
    (0, 'Phytophthora drechsleri')
    >>> parse_fasta_entry('A57915.1 Sequence 20 from Patent EP0751227')
    (0, '')
    >>> parse_fasta_entry('Y08654.1 P.cambivora ribosomal internal ...')
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


assert parse_fasta_entry("LC159493.1 Phytophthora drechsleri genes ...") == (
    0,
    "Phytophthora drechsleri",
)

assert parse_fasta_entry("A57915.1 Sequence 20 from Patent EP0751227") == (0, "")

assert parse_fasta_entry("Y08654.1 P.cambivora ribosomal internal ...") == (0, "")

assert parse_fasta_entry(
    "MG707849.1 Phytophthora humicola x Phytophthora inundata isolate SCVWD597 internal transcribed spacer 1, ..."  # noqa: E501
) == (0, "Phytophthora humicola x inundata")

assert parse_fasta_entry(
    "MG707849.1 Phytophthora humicola x inundata isolate SCVWD597 internal transcribed spacer 1, ..."  # noqa: E501
) == (0, "Phytophthora humicola x inundata")


def main(
    fasta_file,
    db_url,
    min_length=0,
    max_length=sys.maxsize,
    name=None,
    validate_species=False,
    genus_only=False,
    left_primer=None,
    right_primer=None,
    tmp_dir=None,
    sep=";",
    debug=True,
):
    """Implement the ``thapbi_pict ncbi-import`` command."""
    return import_fasta_file(
        fasta_file,
        db_url,
        min_length=min_length,
        max_length=max_length,
        name=name,
        debug=debug,
        fasta_entry_fn=lambda descr: [_.strip() for _ in descr.split(sep)],
        entry_taxonomy_fn=parse_fasta_entry,
        validate_species=validate_species,
        genus_only=genus_only,
        left_primer=left_primer,
        right_primer=right_primer,
        tmp_dir=tmp_dir,
    )
