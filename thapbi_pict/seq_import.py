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

from .utils import find_paired_files


def main(
    inputs,
    method,
    db_url,
    min_abundance=1000,
    name=None,
    validate_species=False,
    debug=True,
):
    """Implement the thapbi_pict seq-import command."""
    assert isinstance(inputs, list)

    input_list = find_paired_files(inputs, ".fasta", ".%s.tsv" % method, debug=debug)

    sys.stderr.write(
        "Importing %i FASTA files with %s classifications\n" % (len(input_list), method)
    )

    return 0
