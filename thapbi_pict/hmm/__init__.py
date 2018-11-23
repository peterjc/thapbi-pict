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
import subprocess
import sys
import tempfile


def find_motif(hmm_file, fasta_input_file, hmmscan='hmmscan', debug=False):
    """Run hmmscan and parse the output.

    This will run the command line tool ``hmmscan`` (assumed to be on the
    ``$PATH`` but an explicit path can be given), writing the output to a
    temporary file which will be deleted after parsing it.

    It will yield tuples for each HMM match,
     - sequence identifier (string)
     - motif sequence (string)
     - number of matches (integer)
     - ...

    The match details are potentially useful if you need to generate unique
    names for each HMM match.
    """
    if not os.path.isfile(fasta_input_file):
        raise ValueError("Missing FASTA input %r" % fasta_input_file)
    if not os.path.isfile(hmm_file):
        raise ValueError("Missing HMM input %r" % hmm_file)

    tmp_dir = tempfile.mkdtemp()
    hmm_out = os.path.join(tmp_dir, "hmmscan.txt")
    cmd = [hmmscan,
           "--noali",
           "-o", hmm_out,
           hmm_file,
           fasta_input_file]
    # cmd = "'%s' --noali -o '%s' '%s' '%s'" % (hmmscan, hmm_out, hmm_file,
    #                                           fasta_input_file)
    if debug:
        sys.stderr.write("DEBUG: Executing command: %s\n" % " ".join(cmd))
    subprocess.check_call(cmd)


if __name__ == "__main__":
    # For testing
    if len(sys.argv) > 2:
        hmm = sys.argv[1]
        for fasta in sys.argv[2:]:
            print("Checking %s against %s" % (fasta, hmm))
            for motif in find_motif(hmm, fasta, debug=True):
                print(motif)
