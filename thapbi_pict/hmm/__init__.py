# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

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
import sys
import tempfile
import warnings

from ..utils import run

from Bio import BiopythonExperimentalWarning
from Bio import SeqIO

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    # This experimenal code (beta) warning was dropped in Biopython 1.73
    from Bio import SearchIO


def run_and_parse_hmmscan(hmm_file, fasta_input_file, hmmscan="hmmscan", debug=False):
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
    cmd = [hmmscan, "--noali", "--cut_ga"]
    cmd += ["-o", hmm_out, hmm_file, fasta_input_file]
    # cmd = "'%s' --noali -o '%s' '%s' '%s'" % (hmmscan, hmm_out, hmm_file,
    #                                           fasta_input_file)
    run(cmd, debug=debug, attempts=3)

    if debug:
        sys.stderr.write("DEBUG: hmmscan finished, about to parse the output\n")

    for record, result in zip(
        SeqIO.parse(fasta_input_file, "fasta"), SearchIO.parse(hmm_out, "hmmer3-text")
    ):
        if record.id != result.id:
            raise ValueError("FASTA %s vs HMMER %s" % (record.id, result.id))
        yield record, result

    os.remove(hmm_out)
    shutil.rmtree(tmp_dir)


def hmm_cache(hmm_file, cache_dir, debug=False):
    """Cache HMM files for calling hmmscan."""
    old_dir, stem = os.path.split(hmm_file)

    new = os.path.join(cache_dir, stem)
    if os.path.isfile(new):
        # Already cached
        return new

    if debug:
        sys.stderr.write("DEBUG: Caching HMM at %s\n" % new)

    # Copy the extra files first:
    for f in os.listdir(old_dir):
        if f.startswith(stem + ".h"):
            shutil.copyfile(os.path.join(old_dir, f), os.path.join(cache_dir, f))

    # Do this last in case interupted (as later only check this exists):
    shutil.copyfile(hmm_file, new)

    return new


def filter_for_ITS1(
    input_fasta, cache_dir, min_length=100, max_length=250, debug=False
):
    """Search for the ITS1 sequence(s) within FASTA entries.

    Expect one ITS1 match in the vast majority of cases.
    """
    hmm = os.path.join(os.path.split(__file__)[0], "combined.hmm")
    if cache_dir:
        hmm = hmm_cache(hmm, cache_dir, debug=debug)

    for record, result in run_and_parse_hmmscan(hmm, input_fasta, debug=debug):
        title = record.description
        seq = str(record.seq).upper()
        if len(record) < min_length or not result:
            yield title, seq, None, []
            continue

        # Not interested in cases like this with no actual hits:
        # [No individual domains that satisfy reporting thresholds
        # (although complete target did)]
        if debug and len([_ for _ in result if _]) > 1:
            sys.stderr.write(
                "DEBUG: %s matched HMM for %s\n"
                % (record.id, ";".join(hit.id for hit in result))
            )

        its1_seqs = [
            (hit.id, seq[hsp.query_start : hsp.query_end])
            for hit in result
            for hsp in hit
        ]

        # Apply length filter to all the matches
        if min_length:
            its1_seqs = [_ for _ in its1_seqs if min_length <= len(_[1])]
        if max_length:
            its1_seqs = [_ for _ in its1_seqs if len(_[1]) <= max_length]

        if not its1_seqs:
            yield title, seq, None, []
            continue

        if len({_[0] for _ in its1_seqs}) > 1:
            # Depending on the orthogonality of the HMM set, this could be fine
            # (e.g. HMM for close sister genera), or a potential problen
            # (e.g. two synthetic controls)
            sys.stderr.write(
                "ERROR: Conflicting HMM matches for %s: %s\n"
                % (record.id, ";".join(sorted({_[0] for _ in its1_seqs})))
            )
            sys.stderr.write("%s length %s:\n" % (_[0], len(_[1])) for _ in its1_seqs)
            sys.exit(1)
        name = its1_seqs[0][0]  # Just checked all the same name

        # Discard HMM names, keep just acceptable length sub-sequences
        its1_seqs = [_[1] for _ in its1_seqs]

        if len(set(its1_seqs)) < len(its1_seqs):
            # e.g. DQ641247.1 Phytophthora ramorum isolate 2195
            sys.stderr.write(
                "WARNING: Discarding exactly duplicated ITS1 matches in %s\n"
                % title.split(None, 1)[0]
            )
            new = []
            for _ in its1_seqs:
                if _ not in new:
                    new.append(_)
            its1_seqs = new
        yield title, seq, name, its1_seqs
