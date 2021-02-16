# Copyright 2019-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Code for THAPBI PICT to deal with importing classified sequences.

Assumes you have Illumina sequenced single isolates. Takes the intermediate
FASTA files from the prepared read step, and their known classifications, and
outputs FASTA files with the species names in the record description - ready
for importing into a database with the ``thapbi_pict curated-import`` command.

This could be used to import classified sequences using one DB into
another DB, but that is not the motivating use case.

The ``thapbi_pict assess`` command can compare two sets of classifier
predictions, but defaults to comparing to set of a "known" values which
can be created for positive controls (e.g. sequencing a plate where
the samples are all from single isolates).

You can also give these positive control "known" classifications to
``thapbi_pict curated-seq`` to then import into a database. The idea here
is to capture experimentally real biological variants beyond the
single canonical ITS1 sequence typically available for each species.
In order to guard against importing PCR artefacts or cross-sample
contamination, you can set a minimum abundance for importing.
"""
import os
import shutil
import sys
import tempfile

from Bio.SeqIO.FastaIO import SimpleFastaParser

from .utils import abundance_from_read_name
from .utils import find_paired_files
from .utils import parse_species_tsv


def main(
    inputs,
    out_dir,
    method,
    min_abundance=1000,
    ignore_prefixes=None,
    tmp_dir=None,
    debug=True,
):
    """Implement the ``thapbi_pict curated-seq`` command."""
    assert isinstance(inputs, list)

    input_list = find_paired_files(
        inputs, ".fasta", f".{method}.tsv", ignore_prefixes, debug, strict=False
    )

    if not input_list:
        sys.exit(
            f"ERROR: Need *.fasta files with matching *.{method}.tsv classification\n"
        )

    sys.stderr.write(
        f"Importing {len(input_list)} FASTA files with {method} classifications\n"
    )

    if out_dir and out_dir != "-" and not os.path.isdir(out_dir):
        sys.stderr.write(f"Making output directory {out_dir!r}\n")
        os.mkdir(out_dir)

    if tmp_dir:
        # Up to the user to remove the files
        tmp = tmp_dir
    else:
        tmp_obj = tempfile.TemporaryDirectory()
        tmp = tmp_obj.name

    for fasta_file, tsv_file in input_list:
        if debug:
            sys.stderr.write(f"DEBUG: Loading meta-data from {tsv_file}\n")

        wild_genus_species = None  # get from TSV
        dict_genus_species = {}  # get from TSV
        try:
            for idn, _, genus_species in parse_species_tsv(tsv_file):
                if idn in dict_genus_species:
                    sys.exit(f"ERROR: Duplicated identifier {idn!r} in {tsv_file!r}")
                dict_genus_species[idn] = genus_species
            if not dict_genus_species:
                sys.stderr.write(
                    f"File {tsv_file} has no sequences, ignoring {fasta_file}\n"
                )
                continue

        except ValueError as e:
            if str(e) != "Wildcard species name found":
                raise
            with open(tsv_file) as handle:
                for line in handle:
                    if line.startswith("*\t"):
                        _, taxid, wild_genus_species, _ = line.split("\t", 3)
            assert wild_genus_species, "Didn't find expected wildcard species line"

        with open(fasta_file) as handle:
            if out_dir == "-":
                output_name = None
                tmp_output = None
                output_handle = sys.stdout
            else:
                # Use the FASTA input name, but in different folder:
                output_name = os.path.join(out_dir, os.path.split(fasta_file)[1])
                tmp_output = os.path.join(tmp, os.path.split(fasta_file)[1])
                output_handle = open(tmp_output, "w")

            for title, seq in SimpleFastaParser(handle):
                idn = title.split(None, 1)[0]
                if abundance_from_read_name(idn) < min_abundance:
                    continue
                genus_species = dict_genus_species.get(idn, wild_genus_species)
                if not genus_species:
                    sys.exit(f"ERROR: No species for {idn} in {fasta_file}")
                output_handle.write(f">{idn} {genus_species}\n{seq}\n")

        if out_dir != "-":
            output_handle.close()
            # Move our temp file into position...
            shutil.move(tmp_output, output_name)

    if tmp_dir:
        # Should be empty as moved files out into place...
        pass
        # sys.stderr.write(
        #     f"WARNING: Please remove temporary files written to {tmp_dir}\n"
        # )
    else:
        tmp_obj.cleanup()

    sys.stderr.write(f"Processed {len(input_list)} FASTA files\n")
    return 0
