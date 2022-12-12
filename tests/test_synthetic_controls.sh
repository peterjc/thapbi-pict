#!/bin/bash

# Copyright 2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -euo pipefail
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp/thapbi_pict}/synthetic_controls
rm -rf $TMP
mkdir -p $TMP

echo "=============================================="
echo "Checking prepare-reads with synthetic controls"
echo "=============================================="

mkdir $TMP/warning
thapbi_pict prepare-reads -a 2 -f 0 \
    -i tests/reads/DNAMIX_S95_L001_R?_001.fastq.gz \
    -o $TMP/warning

thapbi_pict sample-tally \
    -i $TMP/warning/ITS1/DNAMIX_S95_L001.fasta \
    -y $TMP/warning/ITS1/DNAMIX_S95_L001.fasta \
    -o $TMP/warning/tally.tsv \
    2>&1 | grep "WARNING: Control DNAMIX_S95_L001 suggests overly high fractional abundance threshold 23.6"

echo "------------------"
echo "Four plate example"
echo "------------------"

rm -rf $TMP/mock_plates/
mkdir -p $TMP/mock_plates/merged $TMP/mock_plates/prepared $TMP/mock_plates/intermediate_a2

for PLATE in A B C D; do
    # Making mock plates, each with a sample pair and a control pair
    mkdir -p $TMP/mock_plates/plate-${PLATE}
    # Setup the biological sample pair
    cp tests/reads/DNAMIX_S95_L001_R1_001.fastq.gz \
       $TMP/mock_plates/plate-${PLATE}/sample-${PLATE}_R1.fastq.gz
    cp tests/reads/DNAMIX_S95_L001_R2_001.fastq.gz \
       $TMP/mock_plates/plate-${PLATE}/sample-${PLATE}_R2.fastq.gz
    # Create empty FASTQ pair to setup mock control input
    mkdir -p -p $TMP/mock_plates/plate-${PLATE}
    echo | gzip > $TMP/mock_plates/plate-${PLATE}/spike-in-${PLATE}_R1.fastq.gz
    echo | gzip > $TMP/mock_plates/plate-${PLATE}/spike-in-${PLATE}_R2.fastq.gz
    # The merged cache uses gzipped deduplicated FASTA files:
    cat tests/synthetic_controls/spike-in-${PLATE}.fasta \
        | gzip > $TMP/mock_plates/merged/spike-in-${PLATE}.fasta.gz
done

echo "Checking spike-in controls used via pipeline:"

rm -rf $TMP/mock_plates/prepared/*
thapbi_pict pipeline -d - -a 75 -f 0.001 \
            -i $TMP/mock_plates/plate-* \
            -n $TMP/mock_plates/plate-*/spike-in-* \
            --merged-cache $TMP/mock_plates/merged/ \
            -s $TMP/mock_plates/prepared/ \
            -m 1s3g -o $TMP/mock_plates/report

for f in tests/synthetic_controls/report.*.tsv; do
    echo diff $TMP/mock_plates/${f##*/} $f
    diff $TMP/mock_plates/${f##*/} $f
done

echo "Checking spike-in controls used via sample-tally:"

rm -rf $TMP/mock_plates/intermediate_a2/*
thapbi_pict prepare-reads -d - -a 2 -f 0 \
            -i $TMP/mock_plates/plate-* \
            --merged-cache $TMP/mock_plates/merged/ \
            -o $TMP/mock_plates/intermediate_a2/
thapbi_pict sample-tally -d - -a 75 -f 0.001 \
            -i $TMP/mock_plates/intermediate_a2/ITS1/*.fasta \
            -n $TMP/mock_plates/intermediate_a2/ITS1/spike-in-*.fasta \
            -o $TMP/mock_plates/tally.tsv
echo diff $TMP/mock_plates/tally.tsv tests/synthetic_controls/report.ITS1.tally.tsv
diff $TMP/mock_plates/tally.tsv tests/synthetic_controls/report.ITS1.tally.tsv
# Effectively just did this:
echo diff $TMP/mock_plates/tally.tsv $TMP/mock_plates/report.ITS1.tally.tsv
diff $TMP/mock_plates/tally.tsv $TMP/mock_plates/report.ITS1.tally.tsv

echo "--------------------"
echo "Single plate example"
echo "--------------------"

rm -rf $TMP/single_plate/
mkdir -p $TMP/single_plate/raw_data/ $TMP/single_plate/merged $TMP/single_plate/prepared

cp tests/reads/DNAMIX_S95_L001_R1_001.fastq.gz \
   $TMP/single_plate/raw_data/sample_R1.fastq.gz
cp tests/reads/DNAMIX_S95_L001_R2_001.fastq.gz \
   $TMP/single_plate/raw_data/sample_R2.fastq.gz
for PLATE in A B C D; do
    # Create empty FASTQ pair to setup mock control input
    echo | gzip > $TMP/single_plate/raw_data/spike-in-${PLATE}_R1.fastq.gz
    echo | gzip > $TMP/single_plate/raw_data/spike-in-${PLATE}_R2.fastq.gz
    # The merged cache uses gzipped deduplicated FASTA files:
    cat tests/synthetic_controls/spike-in-${PLATE}.fasta \
        | gzip > $TMP/single_plate/merged/spike-in-${PLATE}.fasta.gz
done

thapbi_pict prepare-reads -d - -a 2 -f 0 \
            -i $TMP/single_plate/raw_data/ \
            --merged-cache $TMP/single_plate/merged/ \
            -o $TMP/single_plate/prepared/

echo "Checking spike-in controls..."

# Should all be same as above (except for the #threshold_pool: path)
for PLATE in A B C D; do
    diff \
      <(grep -v "^#threshold_pool:" $TMP/single_plate/prepared/ITS1/spike-in-${PLATE}.fasta) \
      <(grep -v "^#threshold_pool:" $TMP/mock_plates/prepared/ITS1/spike-in-${PLATE}.fasta)
done

echo "Checking the mock sample and threshold used..."

# Should be same as plate D above since that had the highest threshold:
diff \
  <(grep -v "^#threshold_pool:" $TMP/single_plate/prepared/ITS1/sample.fasta) \
  <(grep -v "^#threshold_pool:" $TMP/mock_plates/prepared/ITS1/sample-D.fasta)

echo "Checking spike-in controls used via sample-tally:"

rm -rf $TMP/single_plate/prepared/*
thapbi_pict prepare-reads -d - -a 2 -f 0 \
            -i $TMP/single_plate/raw_data/ \
            --merged-cache $TMP/single_plate/merged/ \
            -o $TMP/single_plate/prepared/
thapbi_pict sample-tally -d - -a 75 -f 0.001 \
            -i $TMP/single_plate/prepared/ITS1/*.fasta \
            -n $TMP/single_plate/prepared/ITS1/spike-in-*.fasta \
            -o $TMP/single_plate/tally.tsv
# Very similar to four-plate tests/synthetic_controls/report.ITS1.tally.tsv
# if ignores sample-A, sample-B and sample-C and compare sample-D with the
# lone sample (i.e. compare the samples with the higher 107 not 75 threshold)
echo diff $TMP/single_plate/tally.tsv tests/synthetic_controls/single-plate.ITS1.tally.tsv
diff $TMP/single_plate/tally.tsv tests/synthetic_controls/single-plate.ITS1.tally.tsv

echo "===="
echo "Done"
echo "===="

echo "$0 - test_synthetic_controls.sh passed"
