#!/usr/bin/env python
"""Turn BLAST TSV output into a FASTA file for curated import.

Example usage to create an ad-hoc DB based on close NCBI NT BLAST matches:

    $ blastn -qcov_hsp_perc 100 -perc_identity 99 \
      -query example.all_reads.fasta -db nt \
      -out example.all_reads.blast.tsv \
      -outfmt "6 qseqid saccver ssciname sseq" \
      -num_threads 8

    $ python blast_to_fasta.py example.all_reads.blast.tsv > example_db.fasta

    $ thapbi_pict curated-import -d example_db.sqlite -i example_db.fasta -x

For brevity this example loads the species names in lax mode. You may prefer
to first load the NCBI taxonomy, and drop the -x import argument.

You are likely to require manual curation for best results - for example
removing problematic entries with suspect species assignments.

WARNING: This may result in duplicated identifiers where the end points
of a matched sequence differ slightly from the BLAST matches.
"""
import os
import sys

if len(sys.argv) != 2 or not os.path.isfile(sys.argv[1]):
    sys.exit("ERROR: Require a BLAST TSV filename as input.")

seq_to_entry = {}
for line in open(sys.argv[1]):
    parts = line.rstrip("\n").split("\t")
    if len(parts) != 4:
        sys.exit(
            "ERROR: Wrong column count, please use:\n"
            'blastn -outfmt "6 qseqid saccver ssciname sseq" ...'
        )

    qseqid, saccver, ssciname, sseq = parts

    if "," in saccver or ";" in saccver:
        sys.exit(
            "ERROR: Unexpected punctuation in saccver field (column 2):\n" + saccver
        )
    if "," in ssciname or ";" in ssciname:
        sys.exit(
            "ERROR: Unexpected punctuation in ssciname field (column 3):\n" + ssciname
        )
    if "synthetic construct" in ssciname:
        continue
    sseq = sseq.upper().replace("-", "")
    entry = f"{saccver} {ssciname}"
    assert "," not in entry, line
    assert ";" not in entry, line
    try:
        seq_to_entry[sseq].add(entry)
    except KeyError:
        seq_to_entry[sseq] = {entry}

for sseq, entry in sorted(seq_to_entry.items()):
    print(">" + ";".join(sorted(entry)))
    print(sseq)
