#!/usr/bin/env python
"""Reformat a legacy curated FASTA file, and sort by sequence.

Mimics a cut-down NCBI style::

    >accession genus species etc
    sequence

Follows NCBI BLAST non-redundant (NR) FASTA style of using
Ctrl+A separated multiple entries on the FASTA title, e.g.::

    >accession1 genus species etc(Ctrl+A)accession2 genus species etc
    sequence

Sadly the Ctrl+A is not a printed character, and it may appear
as invisible in some displays.

Note this discards any Phytophthora clade information.

Note this discards the four control sequences (if present).
"""

import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser

from thapbi_pict.legacy import split_composite_entry, parse_fasta_entry

CTRL_A = chr(1)  # NCBI use for multiple entries in BLAST NR FASTA

sequences = {}
for title, seq in SimpleFastaParser(sys.stdin):
    if title.startswith("Control_"):
        continue
    if seq not in sequences:
        sequences[seq] = []
    for entry in split_composite_entry(title):
        # print(entry)
        acc = entry.split(" ", 1)[0].strip("_")
        if "_" in acc:
            acc = acc.rsplit("_", 1)[1]
        sp = parse_fasta_entry(entry)[1]
        if not sp:
            sp = "Phytophthora"
        else:
            assert sp.startswith("Phytophthora"), entry
        sequences[seq].append("%s %s" % (acc, sp))

for seq in sorted(sequences):
    names = sequences[seq]
    print(">%s\n%s" % (CTRL_A.join(names), seq))
