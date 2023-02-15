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

echo =============================
echo Exploring threshold for mocks
echo =============================

# Not bothering with the metadata here,
# using same -s intermediate/ folder as that's
# before the abundance threshold is applied
thapbi_pict pipeline -a 2 -f 0 -i expected \
        raw_data/DNA15MIX_R?.fastq.gz \
        raw_data/DNA10MIX_undiluted_R?.fastq.gz \
        -s intermediate -o summary/mocks_a2

# Crude parameter sweep of the abundance threshold (-a) with -f 0 assumed.
echo -e "#Threshold\tTP\tFP\tFN\tTN\tsensitivity\tspecificity\tprecision\tF1\tHamming-loss\tAd-hoc-loss" > summary/mocks_a2.assess-vs-abundance.tsv
for A in 2 10 20 30 40 50 60 70 80 90 100; do
    # Name files explicitly to avoid warnings about unused entries
    thapbi_pict assess -i \
            summary/mocks_a2.ITS1.tally.tsv \
            summary/mocks_a2.ITS1.onebp.tsv \
            expected/DNA15MIX.known.tsv \
            expected/DNA10MIX_undiluted.known.tsv \
            -a $A | grep OVERALL | sed "s/OVERALL/A=$A/g" \
            >> summary/mocks_a2.assess-vs-abundance.tsv
done
cut -f 1-5,9,11 summary/mocks_a2.assess-vs-abundance.tsv

echo ====
echo Done
echo ====
