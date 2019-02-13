"""Summarise ITS1 classification results at folder (plate) level.

This implements the ``thapbi_pict plate-summary ...`` command.
"""

import os
import sys
import tempfile

from collections import Counter

from Bio.SeqIO.FastaIO import SimpleFastaParser

from .utils import find_requested_files
from .utils import parse_species_tsv
from .utils import split_read_name_abundance


def main(inputs, output, method, min_abundance=1, debug=False):
    """Implement the thapbi_pict plate-summary command.

    The expectation is that the inputs represent all the samples
    from one (96 well) plate, or some other meaningful batch.
    """
    assert isinstance(inputs, list)

    # Context manager should remove the temp dir:
    with tempfile.TemporaryDirectory() as shared_tmp:
        if debug:
            sys.stderr.write("DEBUG: Shared temp folder %s\n" % shared_tmp)

        samples = set()
        md5_abundance = Counter()
        abundance_by_samples = dict()
        md5_species = dict()
        md5_to_seq = dict()

        if debug:
            sys.stderr.write("Loading FASTA sequences and abundances\n")
        for fasta_file in find_requested_files(
            [_ for _ in inputs if not _.endswith(".tsv")], ".fasta", debug
        ):
            sample = os.path.basename(fasta_file).rsplit(".", 1)[0]
            samples.add(sample)
            with open(fasta_file) as handle:
                for title, seq in SimpleFastaParser(handle):
                    md5, abundance = split_read_name_abundance(title.split(None, 1)[0])
                    if min_abundance > 1 and abundance < min_abundance:
                        continue
                    abundance_by_samples[md5, sample] = abundance
                    md5_abundance[md5] += abundance
                    md5_to_seq[md5] = seq
                    md5_species[md5] = set()
        samples = sorted(samples)

        methods = method.split(",")
        for method in methods:
            if debug:
                sys.stderr.write("Loading predictions for %s\n" % method)
            tsv_files = find_requested_files(
                [_ for _ in inputs if not _.endswith(".fasta")],
                ".%s.tsv" % method,
                debug,
            )
            if len(samples) != len(tsv_files):
                sys.exit(
                    "Identified %i samples from FASTA files, but %i TSV files for %s"
                    % (len(samples), len(tsv_files), method)
                )
            for predicted_file in tsv_files:
                sample = os.path.basename(predicted_file).rsplit(".", 2)[0]
                assert sample in samples, predicted_file
                # TODO: Look at taxid here?
                for name, _, sp in parse_species_tsv(predicted_file, min_abundance):
                    md5, abundance = split_read_name_abundance(name)
                    if min_abundance > 1 and abundance < min_abundance:
                        continue
                    assert abundance_by_samples[md5, sample] == abundance, name
                    # Combining over all methods!
                    if sp:
                        md5_species[md5].update(sp.split(";"))

        if output == "-":
            if debug:
                sys.stderr.write("DEBUG: Output to stdout...\n")
            handle = sys.stdout
        else:
            handle = open(output, "w")

        handle.write(
            "#ITS1-MD5\t%s-predictions\tSequence\tSample-count\tTotal-abundance\t%s\n"
            % (",".join(methods), "\t".join(samples))
        )
        handle.write(
            "TOTAL\t-\t-\t%i\t%i\t%s\n"
            % (
                sum(
                    1
                    for md5 in md5_to_seq
                    for sample in samples
                    if (md5, sample) in abundance_by_samples
                ),
                sum(md5_abundance.values()),
                "\t".join(
                    str(
                        sum(
                            abundance_by_samples.get((md5, sample), 0)
                            for md5 in md5_to_seq
                        )
                    )
                    for sample in samples
                ),
            )
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
