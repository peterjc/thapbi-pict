#!/usr/bin/env python
"""Python 3 script to reformat abundance annotation in FASTA files.

Expected input is FASTA records like ">identifier_abundance [optional text]"
with the output being ">identifier;size=abundance [optional text]" instead.

    $ python swarm2usearch.py swarm.fasta > usearch.fasta

"""
import os
import sys

if len(sys.argv) != 2 or not os.path.isfile(sys.argv[1]):
    sys.exit("ERROR: Expects a single argument, a FASTA filename.")


def fix(header):
    """Switch header style for embedding read abundance."""
    try:
        identifier_abundance, rest = header.split(None, 1)
    except ValueError:
        identifier_abundance = header
        rest = None
    identifier, abundance = identifier_abundance.rsplit("_", 1)
    abundance = int(abundance)
    if rest:
        return f"{identifier};size={abundance} {rest}"
    else:
        return f"{identifier};size={abundance}"


with open(sys.argv[1]) as handle:
    for line in handle:
        if line.startswith(">"):
            sys.stdout.write(f">{fix(line[1:].strip())}\n")
        else:
            sys.stdout.write(line)
