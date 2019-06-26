# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

r"""Code for THAPBI PICT to deal with importing NCBI FASTA files.

For example, you might perform this search via the NCBI Entrez
website, against the nucleotide database, then send it to file
in FASTA format::

   its1 AND Phytophthora[Organism] AND 150:800[Sequence Length]

Or, using the NCBI Entrez Direct command line tools::

    esearch -db nucleotide \
            -query "its1 AND Phytophthora[Organism] AND 150:800[Sequence Length]" \
            | efetch -format fasta > example.fasta

Then, import this into our ITS DB using::

    thapbi_pict ncbi-import -d example.sqlite example.fasta

"""  # noqa: E501

import sys

from .db_import import import_fasta_file


def want_fasta_entry(text):
    """Determine if an NCBI entry should be imported or not.

    Note unlike the legacy FASTA import code, each FASTA record
    is a single entity, so if wanted returns a list of one.
    """
    if text.split(None, 1)[0].lower == "uncultured":
        # We can't trust this, even at genus level - reject
        return []
    else:
        return [text]


def parse_fasta_entry(text):
    """Split an entry of Accession_Genus_Species_name_Description.

    Returns a tuple: taxid (always zero), clade (always empty),
    presumed genus-species (here taken as two words by default),
    and spare text which might be more of the species (for use
    with species name validation).

    Note we can't infer the clade without looking up the species,
    so for now this returns an empty clade.

    >>> parse_fasta_entry('LC159493.1 Phytophthora drechsleri genes ...')
    ('', 'Phytophthora drechsleri', 'genes ...')

    Dividing the species name into genus, species, strain etc
    is not handled here.
    """  # noqa: E501
    parts = text.rstrip().split()
    taxid = 0
    clade = ""
    # acc = parts[0]
    name = parts[1:3]  # assumes "Genus species" only (2 words)
    rest = parts[3:]
    if len(name[0]) > 2 and name[0].startswith("P."):
        # Special case, but can we assume these are Phytophthora?
        # e.g. Y08654.1 P.cambivora ribosomal internal transcribed spacer, ITS1
        sys.stderr.write(
            "WARNING: Assuming %s from %s is Phytophthora\n" % (name[0], parts[0])
        )
        name = ["Phytophthora", name[0][2:]]
    if name[0] == "Sequence":
        # Another special case
        # e.g. A57915.1 Sequence 20 from Patent EP0751227
        name = []
        rest = []
    return (taxid, clade, " ".join(name), " ".join(rest))


assert parse_fasta_entry("LC159493.1 Phytophthora drechsleri genes ...") == (
    0,
    "",
    "Phytophthora drechsleri",
    "genes ...",
)


def main(
    fasta_file, db_url, name=None, validate_species=False, genus_only=False, debug=True
):
    """Implement the thapbi_pict ncbi-import command."""
    return import_fasta_file(
        fasta_file,
        db_url,
        name=name,
        debug=debug,
        fasta_entry_fn=want_fasta_entry,
        entry_taxonomy_fn=parse_fasta_entry,
        validate_species=validate_species,
        genus_only=genus_only,
    )
