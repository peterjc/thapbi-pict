"""Extract HMM matches from sequences (using hmmer3).

This is part of the THAPBI *Phytophthora* ITS1 Classifier Tool (PICT).
We have a Hidden Markov Model (HMM) based on an alignment of manually
curated *Phytophthora* ITS1 sequences, and we will use this to find
ITS1 regions within sequences of potential interest.

The first use case is to trim down ITS1 sequences during import to our
database (for example, sequences submitted to GenBank will often be
longer than the region we are focused on). Second, we can apply this
during our classification pipeline - both to sort sequencing reads into
ITS1 matching or not, but also to trim to a standard region for simpler
comparison.

The code uses Biopython's SearchIO module to parse the output from the
command line tool hmmscan from the hmmer3 suite.
"""

import os
import shutil
import subprocess
import sys
import tempfile
import warnings

from Bio import BiopythonExperimentalWarning
from Bio import SeqIO

with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    # This experimenal code (beta) warning was dropped in Biopython 1.73
    from Bio import SearchIO


def run_and_parse_hmmscan(hmm_file, fasta_input_file, hmmscan='hmmscan',
                          bitscore_threshold=None,
                          debug=False):
    """Run hmmscan and parse the output with Bio.SearchIO.

    This will run the command line tool ``hmmscan`` (assumed to be on the
    ``$PATH`` but an explicit path can be given), writing the output to a
    temporary file which will be deleted after parsing it.

    This will yield 2-tuples of Biopython objects via the FASTA parser in
    Bio.SeqIO, and the HMMER3 text parser in Bio.SearchIO.
    """
    if not os.path.isfile(fasta_input_file):
        raise ValueError("Missing FASTA input %r" % fasta_input_file)
    if not os.path.isfile(hmm_file):
        raise ValueError("Missing HMM input %r" % hmm_file)

    # TODO - Once this is working, refactor to parse stdout
    tmp_dir = tempfile.mkdtemp()
    hmm_out = os.path.join(tmp_dir, "hmmscan.txt")
    cmd = [hmmscan, "--noali"]
    if bitscore_threshold is not None:
        cmd += ["-T", bitscore_threshold, "--domT", bitscore_threshold]
    cmd += ["-o", hmm_out, hmm_file, fasta_input_file]
    # cmd = "'%s' --noali -o '%s' '%s' '%s'" % (hmmscan, hmm_out, hmm_file,
    #                                           fasta_input_file)
    if debug:
        sys.stderr.write("DEBUG: Executing command: %s\n" % " ".join(cmd))
    subprocess.check_call(cmd)

    if debug:
        sys.stderr.write(
            "DEBUG: hmmscan finished, about to parse the output\n")

    for record, result in zip(SeqIO.parse(fasta_input_file, "fasta"),
                              SearchIO.parse(hmm_out, "hmmer3-text")):
        if record.id != result.id:
            raise ValueError("FASTA %s vs HMMER %s" % (record.id, result.id))
        yield record, result

    os.remove(hmm_out)
    shutil.rmtree(tmp_dir)


def filter_for_ITS1(input_fasta, bitscore_threshold="6", debug=False):
    """Search for the expected single ITS1 sequence within FASTA entries.

    The arbitrary low bitscore_threshold default is based on ensuring
    all the sequences in our legacy FASTA files pass without false
    positive additional weak hits complicating things.
    """
    hmm = os.path.join(os.path.split(__file__)[0], "phytophthora_its1.hmm")
    for record, result in run_and_parse_hmmscan(
                hmm, input_fasta,
                bitscore_threshold=bitscore_threshold,
                debug=debug):
        title = record.description
        seq = str(record.seq).upper()
        if not result:
            yield title, seq, None
            continue

        assert len(result) == 1
        hit = result[0]
        assert hit.id == "phytophthora_its1"

        if len(hit) == 0:
            yield title, seq, None
            continue

        yield (title, seq,
               [seq[hsp.query_start:hsp.query_end] for hsp in hit])
