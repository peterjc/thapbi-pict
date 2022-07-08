#!/usr/bin/env python
"""Rename entries in OomyceteDB style FASTA file (via stdin/stdout).

Outputs a FASTA file in the ObiTools style to capture both the
species name and the NCBI taxid.

NOTE: Does not preserve sequence case or any line breaks, gets converted to
upper case as a single line.
"""
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser

for title, seq in SimpleFastaParser(sys.stdin):
    species = title.rstrip().rsplit(";")[-1].replace("_", " ")
    species = species.replace(".", ". ").replace("  ", " ").strip()
    idn = title[title.index("|oodb_id=") :].split("|", 2)[1]
    taxid = title[title.index("|ncbi_taxid=") + 12 :].split("|", 1)[0]
    sys.stdout.write(f">{idn} species_name={species}; taxid={taxid};\n{seq}\n")
