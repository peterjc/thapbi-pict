# Copyright 2019-2020 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Code for THAPBI PICT to import curated FASTA marker files."""
import sys

from .db_import import import_fasta_file
from .utils import find_requested_files


def parse_fasta_entry(text, known_species=None):
    """Split an entry of "Acession genus species etc" into fields.

    Returns a two-tuple of taxid (always zero), genus-species.

    >>> parse_fasta_entry('HQ013219 Phytophthora arenaria')
    (0, 'Phytophthora arenaria')

    >>> parse_fasta_entry('P13660 Phytophthora aff infestans')
    (0, 'Phytophthora aff infestans')
    """
    acc, sp = text.split(None, 1)
    taxid = 0
    # if sp not in known_species:
    #     sys.stderr.write(f"WARNING: Unexpected species name {sp}\n")
    while "  " in sp:
        sp = sp.replace("  ", " ")
    return taxid, sp.strip()


assert parse_fasta_entry("HQ013219 Phytophthora arenaria") == (
    0,
    "Phytophthora arenaria",
)
assert parse_fasta_entry("P13660 Phytophthora aff infestans") == (
    0,
    "Phytophthora aff infestans",
)


def main(
    fasta,
    db_url,
    min_length=0,
    max_length=sys.maxsize,
    name=None,
    validate_species=False,
    genus_only=False,
    left_primer=None,
    right_primer=None,
    sep=";",
    ignore_prefixes=None,
    debug=True,
):
    """Implement the ``thapbi_pict curated-import`` command."""
    fasta_files = find_requested_files(fasta, ".fasta", ignore_prefixes, debug=debug)
    if debug:
        sys.stderr.write(f"Classifying {len(fasta_files)} input FASTA files\n")

    for fasta_file in fasta_files:
        import_fasta_file(
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
        )
