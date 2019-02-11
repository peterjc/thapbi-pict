"""Summarise ITS1 classification results at folder (plate) level.

This implements the ``thapbi_pict plate-summary ...`` command.
"""

import os
import sys
import tempfile

from collections import Counter

from Bio.SeqIO.FastaIO import SimpleFastaParser

from .utils import find_paired_files
from .utils import parse_species_tsv
from .utils import split_read_name_abundance
from .utils import untangle_species


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

        samples = set()
        md5_abundance = Counter()
        abundance_by_samples = dict()
        md5_species = dict()
        md5_to_seq = dict()

        for fasta_file, predicted_file in input_list:
            if debug:
                sys.stderr.write(
                    "DEBUG: Reading %s %s\n" % (fasta_file, predicted_file)
                )
            sample = os.path.basename(fasta_file).rsplit(".", 1)[0]
            samples.add(sample)
            # Load the TSV
            for name, taxid, genus, species in parse_species_tsv(
                predicted_file, min_abundance
            ):
                md5, abundance = split_read_name_abundance(name)
                assert min_abundance < 1 or min_abundance <= abundance, name
                sp_list = untangle_species(taxid, genus, species).split(";")
                try:
                    md5_species[md5].union(sp_list)
                except KeyError:
                    md5_species[md5] = set(sp_list)
                abundance_by_samples[md5, sample] = abundance
                md5_abundance[md5] += abundance
            # Load the FASTA
            with open(fasta_file) as handle:
                for title, seq in SimpleFastaParser(handle):
                    md5, abundance = split_read_name_abundance(title.split(None, 1)[0])
                    if min_abundance > 1 and abundance < min_abundance:
                        continue
                    assert abundance_by_samples[md5, sample] == abundance
                    md5_to_seq[md5] = seq
        samples = sorted(samples)

        if output == "-":
            if debug:
                sys.stderr.write("DEBUG: Output to stdout...\n")
            handle = sys.stdout
        else:
            handle = open(output, "w")

        handle.write(
            "#ITS1-MD5\t%s-predictions\tSequence\tSample-count\tTotal-abundance\t%s\n"
            % (method, "\t".join(samples))
        )
        for total_abundance, md5 in reversed(
            sorted((v, k) for (k, v) in md5_abundance.items())
        ):
            md5_in_xxx_samples = sum(
                1 for _ in samples if (md5, _) in abundance_by_samples
            )
            handle.write(
                "%s\t%s\t%s\t%i\t%i\t%s\n"
                % (
                    md5,
                    ";".join(sorted(md5_species[md5])),
                    md5_to_seq[md5],
                    md5_in_xxx_samples,
                    total_abundance,
                    "\t".join(
                        str(abundance_by_samples.get((md5, _), 0)) for _ in samples
                    ),
                )
            )

        if output != "-":
            handle.close()

        try:
            sys.stdout.flush()
        except BrokenPipeError:
            pass
        try:
            sys.stderr.flush()
        except BrokenPipeError:
            pass

        return 0
