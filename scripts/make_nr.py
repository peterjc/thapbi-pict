#!/usr/bin/env python
"""Make FASTA files non-redundant using semi-colon separator.

Expected usage::

    $ python make_nr.py example1.fasta example2.fasta > nr.fasta

This does not normalise the sequence case, only entries with
the same case-sensitive sequence are merged.

The output is sorted by sequence (case in-sensitive sort).

See also the thapbi_pict fasta-nr and various import commands.
"""
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser

if len(sys.argv) == 1:
    sys.exit("ERROR: Requires one or more FASTA filenames")

sep = ";"
seq_dict = {}

for filename in sys.argv[1:]:
    with open(filename) as handle:
        for title, seq in SimpleFastaParser(handle):
            for entry in title.split(sep):
                if seq in seq_dict:
                    seq_dict[seq].add(entry)
                else:
                    seq_dict[seq] = {entry}

for seq in sorted(seq_dict, key=lambda s: (s.upper(), s)):
    print(">%s\n%s" % (sep.join(sorted(seq_dict[seq])), seq))
