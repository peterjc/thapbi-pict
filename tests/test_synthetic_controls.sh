#!/bin/bash

# Copyright 2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -euo pipefail

export TMP=${TMP:-/tmp/thapbi_pict}/synthetic_controls
rm -rf $TMP
mkdir -p $TMP

echo "=============================================="
echo "Checking prepare-reads with synthetic controls"
echo "=============================================="

echo "------------------"
echo "Four plate example"
echo "------------------"

rm -rf $TMP/mock_plates/
mkdir -p $TMP/mock_plates/merged $TMP/mock_plates/prepared

for PLATE in A B C D; do
    # Making mock plates, each with a sample pair and a control pair
    mkdir -p -p $TMP/mock_plates/plate-${PLATE}
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

thapbi_pict prepare-reads -d - -a 75 \
            -i $TMP/mock_plates/plate-* \
            -n $TMP/mock_plates/plate-*/spike-in-* \
            --merged-cache $TMP/mock_plates/merged/ \
            -o $TMP/mock_plates/prepared/

echo "Checking spike-in controls..."

# A:
if [ `grep -c "^>" $TMP/mock_plates/prepared/ITS1/spike-in-A.fasta` -ne "12" ]; then
    echo "Wrong unique count after abundance threshold in spike-in-A.fasta"; false
fi
if [ `grep "^#abundance:" $TMP/mock_plates/prepared/ITS1/spike-in-A.fasta` != "#abundance:38473" ]; then
    echo "Wrong count accepted after abundance threshold in spike-in-A.fasta"; false
fi
# B:
if [ `grep -c "^>" $TMP/mock_plates/prepared/ITS1/spike-in-B.fasta` -ne "8" ]; then
    echo "Wrong unique count after abundance threshold in spike-in-B.fasta"; false
fi
if [ `grep "^#abundance:" $TMP/mock_plates/prepared/ITS1/spike-in-B.fasta` != "#abundance:84648" ]; then
    echo "Wrong count accepted after abundance threshold in spike-in-B.fasta"; false
fi
# C:
if [ `grep -c "^>" $TMP/mock_plates/prepared/ITS1/spike-in-C.fasta` -ne "7" ]; then
    echo "Wrong unique count after abundance threshold in spike-in-C.fasta"; false
fi
if [ `grep "^#abundance:" $TMP/mock_plates/prepared/ITS1/spike-in-C.fasta` != "#abundance:44501" ]; then
    echo "Wrong count accepted after abundance threshold in spike-in-C.fasta"; false
fi
# D:
if [ `grep -c "^>" $TMP/mock_plates/prepared/ITS1/spike-in-D.fasta` -ne "6" ]; then
    echo "Wrong unique count after abundance threshold in spike-in-D.fasta"; false
fi
if [ `grep "^#abundance:" $TMP/mock_plates/prepared/ITS1/spike-in-D.fasta` != "#abundance:25102" ]; then
    echo "Wrong count accepted after abundance threshold in spike-in-D.fasta"; false
fi

echo "Checking the mock samples and thresholds used..."

# A, threshold kept at 75:
if [ `grep "^#threshold:" $TMP/mock_plates/prepared/ITS1/sample-A.fasta` != "#threshold:75" ]; then
    echo "Wrong abundance threshold in sample-A.fasta"; false
fi
if [ `grep -c "^>" $TMP/mock_plates/prepared/ITS1/sample-A.fasta` -ne "8" ]; then
    echo "Wrong unique count after abundance threshold in sample-A.fasta"; false
fi
if [ `grep "^#abundance:" $TMP/mock_plates/prepared/ITS1/sample-A.fasta` != "#abundance:3683" ]; then
    echo "Wrong count accepted after abundance threshold in sample-A.fasta"; false
fi
# B, threshold kept at 75:
if [ `grep "^#threshold:" $TMP/mock_plates/prepared/ITS1/sample-B.fasta` != "#threshold:75" ]; then
    echo "Wrong abundance threshold in sample-B.fasta"; false
fi
if [ `grep -c "^>" $TMP/mock_plates/prepared/ITS1/sample-B.fasta` -ne "8" ]; then
    echo "Wrong unique count after abundance threshold in sample-B.fasta"; false
fi
if [ `grep "^#abundance:" $TMP/mock_plates/prepared/ITS1/sample-B.fasta` != "#abundance:3683" ]; then
    echo "Wrong count accepted after abundance threshold in sample-B.fasta"; false
fi
# C, threshold kept at 75:
if [ `grep "^#threshold:" $TMP/mock_plates/prepared/ITS1/sample-C.fasta` != "#threshold:75" ]; then
    echo "Wrong abundance threshold in sample-C.fasta"; false
fi
if [ `grep -c "^>" $TMP/mock_plates/prepared/ITS1/sample-C.fasta` -ne "8" ]; then
    echo "Wrong unique count after abundance threshold in sample-C.fasta"; false
fi
if [ `grep "^#abundance:" $TMP/mock_plates/prepared/ITS1/sample-C.fasta` != "#abundance:3683" ]; then
    echo "Wrong count accepted after abundance threshold in sample-C.fasta"; false
fi
# D, threshold raised to 107:
if [ `grep "^#threshold:" $TMP/mock_plates/prepared/ITS1/sample-D.fasta` != "#threshold:107" ]; then
    echo "Wrong abundance threshold in sample-D.fasta"; false
fi
if [ `grep -c "^>" $TMP/mock_plates/prepared/ITS1/sample-D.fasta` -ne "7" ]; then
    echo "Wrong unique count after abundance threshold in sample-D.fasta"; false
fi
if [ `grep "^#abundance:" $TMP/mock_plates/prepared/ITS1/sample-D.fasta` != "#abundance:3585" ]; then
    echo "Wrong count accepted after abundance threshold in sample-D.fasta"; false
fi

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

thapbi_pict prepare-reads -d - -a 75 \
            -i $TMP/single_plate/raw_data/ \
            -n $TMP/single_plate/raw_data/spike-in-* \
            --merged-cache $TMP/single_plate/merged/ \
            -o $TMP/single_plate/prepared/

echo "Checking spike-in controls..."

# Should all be same as above:
for PLATE in A B C D; do
    diff $TMP/single_plate/prepared/ITS1/spike-in-${PLATE}.fasta $TMP/mock_plates/prepared/ITS1/spike-in-${PLATE}.fasta
done

echo "Checking the mock sample and threshold used..."

# Should be same as plate D above since that had the highest threshold:
diff $TMP/single_plate/prepared/ITS1/sample.fasta $TMP/mock_plates/prepared/ITS1/sample-D.fasta

echo "===="
echo "Done"
echo "===="

echo "$0 - test_synthetic_controls.sh passed"
