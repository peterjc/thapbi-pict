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
    name=None,
    validate_species=False,
    genus_only=False,
    debug=True,
):
    """Implement the thapbi_pict seq-import command."""
    assert isinstance(inputs, list)

    input_list = find_paired_files(inputs, ".fasta", ".%s.tsv" % method, debug=debug)

    if not input_list:
        sys.exit(
            "ERROR: Need *.fasta files with matching *.%s.tsv classification\n" % method
        )

    sys.stderr.write(
        "Importing %i FASTA files with %s classifications\n" % (len(input_list), method)
    )

    for fasta_file, tsv_file in input_list:
        if debug:
            sys.stderr.write("DEBUG: Loading meta-data from %s\n" % tsv_file)
        meta_data = {}
        # Apply minimum abundance threshold during FASTA loading
        try:
            for idn, taxid, genus_species in parse_species_tsv(tsv_file):
                if idn in meta_data:
                    sys.exit("ERROR: Duplicated identifier %r in %r" % (idn, tsv_file))
                meta_data[idn] = (int(taxid), "", genus_species, "")
            if not meta_data:
                sys.stderr.write(
                    "File %s has no sequences, ignoring %s\n" % (tsv_file, fasta_file)
                )
                continue
            meta_fn = meta_data.get

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

            def meta_fn(text):
                return (int(taxid), "", genus_species, "")

            def sequence_wanted(title):
                """Check if identifier passess abundance level."""
                idn = title.split(None, 1)[0]
                if abundance_from_read_name(idn) < min_abundance:
                    return []
                else:
                    return [idn]

            sys.stderr.write("File %s wild card for %s\n" % (tsv_file, genus_species))

        import_fasta_file(
            fasta_file,
            db_url,
            name,
            fasta_entry_fn=sequence_wanted,
            entry_taxonomy_fn=meta_fn,
            debug=debug,
            validate_species=validate_species,
            genus_only=genus_only,
        )

    sys.stderr.write("Imported %i FASTA files\n" % len(input_list))
    return 0
