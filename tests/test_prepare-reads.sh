#!/bin/bash

# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "======================"
echo "Checking prepare-reads"
echo "======================"
set -x
thapbi_pict prepare-reads 2>&1 | grep "the following arguments are required"
set -o pipefail

# Try a real example
rm -rf $TMP/DNAMIX_S95_L001.fasta
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 0
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.fasta` -ne "735" ]; then echo "Wrong FASTA output count"; false; fi

# In this case, --flip makes no difference
rm -rf $TMP/DNAMIX_S95_L001.fasta
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 0 --flip
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.fasta` -ne "735" ]; then echo "Wrong FASTA output count"; false; fi

rm -rf $TMP/DNAMIX_S95_L001.fasta
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 5
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.fasta` -ne "24" ]; then echo "Wrong FASTA output count"; false; fi

rm -rf $TMP/DNAMIX_S95_L001.fasta
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 5 --hmm thapbi_pict/hmm/combined.hmm
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.fasta` -ne "24" ]; then echo "Wrong FASTA output count"; false; fi

rm -rf $TMP/DNAMIX_S95_L001.fasta
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 5 --hmm ''
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.fasta` -ne "24" ]; then echo "Wrong FASTA output count"; false; fi

echo "Generating mock control file"
# Using just 50 real reads (50 * 4 = 200 lines)
set +o pipefail
cat tests/reads/DNAMIX_S95_L001_R1_001.fastq.gz | gunzip | head -n 200 > $TMP/MOCK_CONTROL_R1.fastq
cat tests/reads/DNAMIX_S95_L001_R2_001.fastq.gz | gunzip | head -n 200 > $TMP/MOCK_CONTROL_R2.fastq
set -o pipefail

rm -rf $TMP/DNAMIX_S95_L001.fasta
rm -rf $TMP/MOCK_CONTROL.fasta
# Starting low threshold, should be increased to 19, so get new output count...
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 5 -n $TMP/MOCK_CONTROL_R?.fastq
if [ `grep -c "^>" $TMP/MOCK_CONTROL.fasta` -ne "1" ]; then echo "Wrong FASTA control output count"; false; fi
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.fasta` -ne "9" ]; then echo "Wrong FASTA output count"; false; fi

rm -rf $TMP/DNAMIX_S95_L001.fasta
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 100
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.fasta` -ne "7" ]; then echo "Wrong FASTA output count"; false; fi
diff $TMP/DNAMIX_S95_L001.fasta tests/prepare-reads/DNAMIX_S95_L001.fasta

rm -rf $TMP/DNAMIX_S95_L001.fasta
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 100 --hmm ''
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.fasta` -ne "7" ]; then echo "Wrong FASTA output count"; false; fi
# Should be identical but without the HMM names in the FASTA title lines:
diff <(grep -v ">" $TMP/DNAMIX_S95_L001.fasta) <(grep -v "^>" tests/prepare-reads/DNAMIX_S95_L001.fasta)

echo "$0 - test_prepare-reads.sh passed"
