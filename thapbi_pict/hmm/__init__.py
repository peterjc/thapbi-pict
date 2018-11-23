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

from Bio import SearchIO
from Bio import SeqIO


def run_and_parse_hmmscan(hmm_file, fasta_input_file, hmmscan='hmmscan',
                          bitscore_threshold=None,
                          debug=False):
    """Run hmmscan and parse the output with Bio.SearchIO.

    This will run the command line tool ``hmmscan`` (assumed to be on the
    ``$PATH`` but an explicit path can be given), writing the output to a
    temporary file which will be deleted after parsing it.

    This will yield Biopython objects via the HMMER3 parser in Bio.SearchIO.
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

    for record in SearchIO.parse(hmm_out, "hmmer3-text"):
        yield record

    os.remove(hmm_out)
    shutil.rmtree(tmp_dir)


def filter_fasta_for_ITS1(input_fasta, output_fasta,
                          bitscore_threshold="50",
                          debug=False):
    """Search for ITS1 sequences and produce an output of any matches.

    Returns a tuple of counts:
     - number of inputs
     - number of outputs

    Initally enforces at most one ITS1 match per sequence.

    Initially the HMM is hard coded.
    """
    hmm = os.path.join(os.path.split(__file__)[0], "phytophthora_its1.hmm")
    input_count = 0
    output_count = 0

    with open(output_fasta, "w") as out_handle:
        for record, result in zip(
                SeqIO.parse(input_fasta, "fasta"),
                run_and_parse_hmmscan(hmm, input_fasta,
                                      bitscore_threshold=bitscore_threshold,
                                      debug=debug)):
            assert record.id == result.id, "%s vs %s" % (record.id, result.id)
            input_count += 1
            if not result:
                if debug:
                    sys.stderr.write("DEBUG: Ignoring %s\n" % result.id)
                continue
            assert len(result) == 1  # Searching against only one HMM profile
            hit = result[0]
            assert hit.id == "phytophthora_its1"
            if len(hit) == 0:
                if debug:
                    sys.stderr.write("DEBUG: Ignoring %s\n" % result.id)
                continue
            if len(hit) != 1:
                print(hit)
                for dom in hit:
                    print(dom)
                assert False
            assert len(hit) == 1
            hsp = hit[0]
            print(hsp.query_id, hsp.bitscore, hsp.query_start, hsp.query_end)
            assert record.id == hsp.query_id
            old_len = len(record)
            record = record[hsp.query_start:hsp.query_end]
            if len(record) == old_len:
                record.description = "HMM bitscore %s, untrimmed at %ibp" % (
                    hsp.bitscore, old_len)
            else:
                record.description = "HMM bitscore %s, cut %ibp to %ibp" % (
                    hsp.bitscore, old_len, len(record))
            output_count += SeqIO.write(record, out_handle, "fasta")
    return input_count, output_count


if __name__ == "__main__":
    # For testing
    if len(sys.argv) == 3:
        in_fasta = sys.argv[1]
        out_fasta = sys.argv[2]
        in_count, out_count = filter_fasta_for_ITS1(in_fasta, out_fasta,
                                                    debug=True)
        print("Extracted %i from %i inputs" % (out_count, in_count))
