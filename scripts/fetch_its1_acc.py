#!/usr/bin/env python3
"""Python 3 script to fetch accessions and trim to ITS1 region.

While it could be generalised (e.g. configurable primers and trimming
behaviour), this is initially very specific to the Phytophthora
targeting ITS1 primers used in THAPBI PICT's default database. Usage:

    $ python fetch_its1_acc.py ACC1 ... > candidates.fasta

It takes as command line arguments one or more NCBI accessions, and
prints out trimmed FASTA sequences hopefully of their ITS1 amplicon.
"""

import re
import sys

from Bio import Entrez
from Bio.SeqIO.FastaIO import SimpleFastaParser

from thapbi_pict.db_import import parse_ncbi_fasta_entry

leader32bp = "TTTCCGTAGGTGAACCTGCGGAAGGATCATTA"

left_primer = "GAAGGTGAAGTCGTAACAAGG"  # no ambiguities
# The reverse complement of the right primer is GYRGGGACGAAAGTCYYTGC
# which should appear at the end of the region of interest
# Turn into a regular expression, Y=C or T, and R=A or G.
re_right = re.compile("G[CTY][AGR]GGGACGAAAGTC[CTY][CTY]TGC")


def fetch_acc(identifier, out_handle):
    """Download sequence for the accession give, writing FASTA to output handle."""
    title, seq = next(
        SimpleFastaParser(
            Entrez.efetch(db="nucleotide", id=identifier, rettype="fasta")
        )
    )
    acc = title.split(None, 1)[0]
    assert acc.startswith(acc)

    # Trim right primer if present
    try:
        seq = seq[: re_right.search(seq).start()]
    except AttributeError:
        sys.stderr.write(f"ERROR: Right primer not found in {acc}\n")
        return

    # Trim left primer if present in full
    if left_primer in seq:
        seq = seq[seq.index[left_primer] + len(left_primer) :]
    if leader32bp in seq:
        index = seq.index(leader32bp)
        leader = seq[:index]
        if index == 0:
            pass  # no primer (anymore)
        elif index < len(left_primer):
            if not left_primer.endswith(leader):
                sys.stderr.write(
                    f"WARNING: {acc}, got {leader} not primer {left_primer}\n"
                )
        else:
            sys.stderr.write(f"ERROR: {acc}, got {leader} not primer {left_primer}\n")
            return
        seq = seq[index:]
    # TODO - Handle starting with partial 32bp leader automatically?

    # Add assumed 32bp leader if missing
    if seq.startswith("CCACACCTAAAAA"):
        # missing leading 32bp
        seq = leader32bp.lower() + seq
        acc = acc.split(".", 1)[0] + ".X"  # our convention for extended

    species = parse_ncbi_fasta_entry(title)[1]
    out_handle.write(f">{acc} {species}\n{seq}\n")


if len(sys.argv) < 2:
    sys.exit("ERROR: Expects one or more accessions as arguments.")

for acc in sys.argv[1:]:
    fetch_acc(acc, sys.stdout)
