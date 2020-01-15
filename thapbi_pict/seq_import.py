# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

"""Code for THAPBI PICT to deal with importing classified sequences.

This could be used to import classified sequences using one DB into
another DB, but that is not the motivating use case.

The ``thapbi_pict assess`` command can compare two sets of classifier
predictions, but defaults to comparing to set of a "known" values which
can be created for positive controls (e.g. sequencing a plate where
the samples are all from single isolates).

You can also give these positive control "known" classifications to
``thapbi_pict seq_import`` to import into a database. The idea here
is to capture experimentally real biological variants beyond the
single canonical ITS1 sequence typically available for each species.
In order to guard against importing PCR artefacts or cross-sample
contamination, you can set a minimum abundance for importing.
"""

import sys

from .db_import import import_fasta_file
from .utils import abundance_from_read_name
from .utils import find_paired_files
from .utils import parse_species_tsv


def main(
    inputs,
    method,
    db_url,
    min_abundance=1000,
    min_length=0,
    max_length=sys.maxsize,
    name=None,
    validate_species=False,
    genus_only=False,
    ignore_prefixes=None,
    debug=True,
):
    """Implement the ``thapbi_pict seq-import`` command."""
    assert isinstance(inputs, list)

    input_list = find_paired_files(
        inputs, ".fasta", f".{method}.tsv", ignore_prefixes, debug=debug
    )

    if not input_list:
        sys.exit(
            f"ERROR: Need *.fasta files with matching *.{method}.tsv classification\n"
        )

    sys.stderr.write(
        f"Importing {len(input_list):d} FASTA files with {method} classifications\n"
    )

    for fasta_file, tsv_file in input_list:
        if debug:
            sys.stderr.write(f"DEBUG: Loading meta-data from {tsv_file}\n")
        meta_data = {}
        # Apply minimum abundance threshold during FASTA loading
        try:
            for idn, taxid, genus_species in parse_species_tsv(tsv_file):
                if idn in meta_data:
                    sys.exit(f"ERROR: Duplicated identifier {idn!r} in {tsv_file!r}")
                meta_data[idn] = (int(taxid), "", genus_species)
            if not meta_data:
                sys.stderr.write(
                    f"File {tsv_file} has no sequences, ignoring {fasta_file}\n"
                )
                continue

            def meta_fn(text, known_species=None):
                return meta_data[text]

            def sequence_wanted(title):
                """Check if identifier in TSV file, and passess abundance level."""
                idn = title.split(None, 1)[0]
                if idn not in meta_data:
                    return []
                elif abundance_from_read_name(idn) < min_abundance:
                    return []
                else:
                    return [idn]

        except ValueError as e:
            if str(e) != "Wildcard species name found":
                raise
            genus_species = None
            with open(tsv_file) as handle:
                for line in handle:
                    if line.startswith("*\t"):
                        _, taxid, genus_species, _ = line.split("\t", 3)
            assert genus_species, "Didn't find expected wildcard species line"

            def meta_fn(text, known_species=None):
                return int(taxid), genus_species

            def sequence_wanted(title):
                """Check if identifier passess abundance level."""
                idn = title.split(None, 1)[0]
                if abundance_from_read_name(idn) < min_abundance:
                    return []
                else:
                    return [idn]

            sys.stderr.write(f"File {tsv_file} wild card for {genus_species}\n")

        import_fasta_file(
            fasta_file,
            db_url,
            min_length=min_length,
            max_length=max_length,
            name=name,
            fasta_entry_fn=sequence_wanted,
            entry_taxonomy_fn=meta_fn,
            debug=debug,
            validate_species=validate_species,
            genus_only=genus_only,
        )

    sys.stderr.write(f"Imported {len(input_list)} FASTA files\n")
    return 0
