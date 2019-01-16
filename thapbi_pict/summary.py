"""Summarise ITS1 classification results at folder (plate) level.

This implements the ``thapbi_pict plate-summary ...`` command.
"""

import sys
import tempfile


from .utils import find_paired_files


def main(inputs, output, method, min_abundance=1, debug=False):
    """Implement the thapbi_pict plate-summary command.

    The expectation is that the inputs represent all the samples
    from one (96 well) plate, or some other meaningful batch.
    """
    assert isinstance(inputs, list)

    input_list = find_paired_files(inputs, ".fasta", ".%s.tsv" % method, debug=False)

    if debug:
        sys.stderr.write(
            "Reporting on %i samples using %s classifier\n" % (len(input_list), method)
        )

    # Context manager should remove the temp dir:
    with tempfile.TemporaryDirectory() as shared_tmp:
        if debug:
            sys.stderr.write("DEBUG: Shared temp folder %s\n" % shared_tmp)
        for fasta_file, predicted_file in input_list:
            if debug:
                sys.stderr.write(
                    "DEBUG: Reading %s %s\n" % (fasta_file, predicted_file)
                )
