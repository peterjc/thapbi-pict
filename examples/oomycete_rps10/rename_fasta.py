#!/usr/bin/env python
"""Rename entries in OomyceteDB style FASTA file (via stdin/stdout).

NOTE: Does not preserve sequence case or any line breaks, gets converted to
upper case as a single line.
"""
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser

for title, seq in SimpleFastaParser(sys.stdin):
    species = title.rstrip().rsplit(";")[-1].replace("_", " ")
    idn = title[title.find("|oodb_id=") :].split("|", 2)[1]
    sys.stdout.write(f">{idn} {species}\n{seq}\n")
