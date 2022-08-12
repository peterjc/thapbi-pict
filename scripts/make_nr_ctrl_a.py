"""Make FASTA files non-redundant using Ctrl+A separator."""
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser

CTRL_A = chr(1)

seq_dict = {}
for filename in sys.argv[1:]:
    with open(filename) as handle:
        for title, seq in SimpleFastaParser(handle):
            for entry in title.split(CTRL_A):
                if seq in seq_dict:
                    seq_dict[seq].add(entry)
                else:
                    seq_dict[seq] = {entry}

for seq in sorted(seq_dict):
    print(">%s\n%s" % (CTRL_A.join(sorted(seq_dict[seq])), seq))
