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
rm -rf $TMP/DNAMIX_S95_L001.fasta $TMP/DNAMIX_S95_L001.failed-primers.fasta
rm -rf $TMP/merged_cache/
mkdir $TMP/merged_cache/
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz \
    -a 0 --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.fasta` -ne "827" ]; then echo "Wrong FASTA output count"; false; fi

# In this case, --flip makes no difference, as does -n ""
# Also using -p / --primers too.
# Using merged cache also should make no difference
rm -rf $TMP/DNAMIX_S95_L001.fasta $TMP/DNAMIX_S95_L001.failed-primers.fasta
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz \
    --flip -n "" --merged-cache $TMP/merged_cache/ -p $TMP/ \
    -a 0 --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.fasta` -ne "827" ]; then echo "Wrong FASTA output count"; false; fi
diff $TMP/DNAMIX_S95_L001.failed-primers.fasta tests/prepare-reads/DNAMIX_S95_L001.failed-primers.fasta

rm -rf $TMP/DNAMIX_S95_L001.fasta
# Reusing the pre-primer time pre-abundance cache here:
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz \
    -n "-" --merged-cache $TMP/merged_cache/ \
    -a 5 --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.fasta` -ne "27" ]; then echo "Wrong FASTA output count"; false; fi

rm -rf $TMP/DNAMIX_S95_L001.fasta
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz \
    -a 5 --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA --hmm thapbi_pict/hmm/controls.hmm
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.fasta` -ne "27" ]; then echo "Wrong FASTA output count"; false; fi

rm -rf $TMP/DNAMIX_S95_L001.fasta
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz \
    -a 5 --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA --hmm ''
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.fasta` -ne "27" ]; then echo "Wrong FASTA output count"; false; fi

echo "Generating mock control file"
# Using just 50 real reads (50 * 4 = 200 lines)
set +o pipefail
cat tests/reads/DNAMIX_S95_L001_R1_001.fastq.gz | gunzip | head -n 200 > $TMP/MOCK_CONTROL_R1.fastq
cat tests/reads/DNAMIX_S95_L001_R2_001.fastq.gz | gunzip | head -n 200 > $TMP/MOCK_CONTROL_R2.fastq
set -o pipefail

rm -rf $TMP/DNAMIX_S95_L001.fasta
rm -rf $TMP/MOCK_CONTROL.fasta
rm -rf $TMP/pool/
mkdir  $TMP/pool/
cp tests/reads/DNAMIX_S95_L001_*.fastq.gz $TMP/MOCK_CONTROL_R?.fastq $TMP/pool/
# Starting low threshold, should be increased to 19, so get new output count...
thapbi_pict prepare-reads -o $TMP -i $TMP/pool/ -a 5 -n $TMP/pool/MOCK_CONTROL_R?.fastq
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

# Testing primers (default)
rm -rf $TMP/SRR6303948*.fasta
thapbi_pict prepare-reads -o $TMP -i tests/reads/SRR6303948_sample_*.fastq -a 2
diff $TMP/SRR6303948_sample.fasta tests/prepare-reads/SRR6303948_sample_default.fasta

# Testing primers (actual)
rm -rf $TMP/SRR6303948*.fasta
thapbi_pict prepare-reads -o $TMP -i tests/reads/SRR6303948_sample_*.fastq -a 2 \
        --left GAAGGTGAAGTCGTAACAAGG --right AGCGTTCTTCATCGATGTGC
diff $TMP/SRR6303948_sample.fasta tests/prepare-reads/SRR6303948_sample_primers.fasta

echo "Testing --flip works"
# Took only 50 reads from original file (about 6%),
# dropping abudnance threshold from 100 to only 10.

rm -rf $TMP/sample.fasta
thapbi_pict prepare-reads --hmm '' --left CTGCTGCTGGATCATTACCC --right CGCCAGCACAGCCGTTAG --minlen 150 --maxlen 350 \
    -i tests/nematodes/sample_R*.fastq.gz -a 10 -o $TMP/
diff $TMP/sample.fasta tests/nematodes/sample_noflip_a10.fasta  # empty!

rm -rf $TMP/sample.fasta
thapbi_pict prepare-reads --hmm '' --left CTGCTGCTGGATCATTACCC --right CGCCAGCACAGCCGTTAG --minlen 150 --maxlen 350 \
    -i tests/nematodes/sample_R*.fastq.gz -a 10 --flip -o $TMP/
diff $TMP/sample.fasta tests/nematodes/sample_flip_a10.fasta

echo "Testing pathological trimming example"
rm -rf $TMP/6e847180a4da6eed316e1fb98b21218f.fasta
thapbi_pict prepare-reads -i tests/reads/6e847180a4da6eed316e1fb98b21218f_R?.fastq \
    -o $TMP/ -a 1
diff $TMP/6e847180a4da6eed316e1fb98b21218f.fasta tests/prepare-reads/6e847180a4da6eed316e1fb98b21218f.fasta

echo "$0 - test_prepare-reads.sh passed"
