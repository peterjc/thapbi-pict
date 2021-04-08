#!/bin/bash
set -eup pipeline

echo NOTE: Expected first time run time is under 5 minues,
echo repeat runs about 1 minute just to regenerate reports.
echo
echo ==================
echo Woody hosts - ITS1
echo ==================

echo "Running analysis"
mkdir -p intermediate/ summary/

echo "Pipeline without metadata..."
# Using default report naming
thapbi_pict pipeline -i raw_data/ -s intermediate/ -o summary/ \
        -n raw_data/NEGATIVE*.fastq.gz

echo "Pipeline with metadata & assess classifier..."
# Reuses the intermediate files (prepared FASTA and classifer output)
# Giving report name stem (so not to over-write reports without metadata)
thapbi_pict pipeline -i raw_data/ expected/ -s intermediate/ -o summary/ \
        -n raw_data/NEGATIVE*.fastq.gz -r with-metadata \
        -t metadata.tsv -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -x 16 -f 20

echo ====
echo Done
echo ====
