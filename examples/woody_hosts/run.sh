#!/bin/bash
set -euo pipefail

echo "NOTE: Expected first time run time is about a minute,"
echo "repeat runs a few seconds just to regenerate reports."
echo
echo ==================
echo Woody hosts - ITS1
echo ==================

echo "Running analysis"
mkdir -p intermediate/ summary/

echo "Pipeline without metadata..."
thapbi_pict pipeline -i raw_data/ -s intermediate/ \
        -o summary/thapbi-pict -n raw_data/NEGATIVE*.fastq.gz

echo "Pipeline with metadata & assess classifier..."
# Reuses the intermediate files (prepared FASTA files)
# Giving different report name stem (so not to over-write reports without metadata)
thapbi_pict pipeline -i raw_data/ expected/ -s intermediate/ \
        -o summary/with-metadata -n raw_data/NEGATIVE*.fastq.gz \
        -t metadata.tsv -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -x 16

echo ====
echo Done
echo ====
