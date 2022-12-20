#!/bin/bash

# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp/thapbi_pict}/prepare_reads
rm -rf $TMP
mkdir -p $TMP

echo "======================"
echo "Checking prepare-reads"
echo "======================"
set -x
thapbi_pict prepare-reads 2>&1 | grep "the following arguments are required"
set -o pipefail

export DB=$TMP/w32_on_primer.sqlite
rm -rf $DB
thapbi_pict import -d $DB -x -k ITS1 --minlen 100 --maxlen 700 \
    --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA \
    --right GCARRGACTTTCGTCCCYRC -i database/Phytophthora_ITS1_curated.fasta

# Try a real example with the above primers
rm -rf $TMP/ITS1
rm -rf $TMP/merged_cache/
mkdir $TMP/merged_cache/
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz \
    -a 0 -f 0 -d $DB
if [ `grep -c "^>" $TMP/ITS1/DNAMIX_S95_L001.fasta` -ne "822" ]; then echo "Wrong FASTA output count"; false; fi

# In this case, --flip makes no difference
# Using merged cache also should make no difference
rm -rf $TMP/ITS1
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz \
    --flip --merged-cache $TMP/merged_cache/ -a 0 -f 0 -d $DB
if [ `grep -c "^>" $TMP/ITS1/DNAMIX_S95_L001.fasta` -ne "822" ]; then echo "Wrong FASTA output count"; false; fi

rm -rf $TMP/ITS1
# Reusing the pre-primer-trim pre-abundance cache here:
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz \
    --merged-cache $TMP/merged_cache/ -a 5 -d $DB
if [ `grep -c "^>" $TMP/ITS1/DNAMIX_S95_L001.fasta` -ne "27" ]; then echo "Wrong FASTA output count"; false; fi

rm -rf $TMP/ITS1
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz \
    -a 5 -d $DB --database '-'
if [ `grep -c "^>" $TMP/ITS1/DNAMIX_S95_L001.fasta` -ne "27" ]; then echo "Wrong FASTA output count"; false; fi

rm -rf $TMP/ITS1
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz \
    -a 5 -d $DB
if [ `grep -c "^>" $TMP/ITS1/DNAMIX_S95_L001.fasta` -ne "27" ]; then echo "Wrong FASTA output count"; false; fi

echo "---------------------"
echo "With negative control"
echo "---------------------"
echo "Generating mock control file"
rm -rf $TMP/pool/
mkdir  $TMP/pool/
cp tests/reads/DNAMIX_S95_L001_*.fastq.gz $TMP/pool/
# Using just 50 real reads (50 * 4 = 200 lines), and alternative extension
set +o pipefail
cat tests/reads/DNAMIX_S95_L001_R1_001.fastq.gz | gunzip | head -n 200 > $TMP/pool/MOCK_CONTROL_R1.fq
cat tests/reads/DNAMIX_S95_L001_R2_001.fastq.gz | gunzip | head -n 200 > $TMP/pool/MOCK_CONTROL_R2.fq
set -o pipefail

rm -rf $TMP/ITS1
thapbi_pict prepare-reads -o $TMP -i tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 100
if [ `grep -c "^>" $TMP/ITS1/DNAMIX_S95_L001.fasta` -ne "7" ]; then echo "Wrong FASTA output count"; false; fi
diff $TMP/ITS1/DNAMIX_S95_L001.fasta tests/prepare-reads/DNAMIX_S95_L001.fasta

# Testing primers (default)
rm -rf $TMP/ITS1
thapbi_pict prepare-reads -o $TMP -i tests/reads/SRR6303948_sample_*.fastq -a 2
diff $TMP/ITS1/SRR6303948_sample.fasta tests/prepare-reads/SRR6303948_sample_default.fasta

# Testing primers (actual)
rm -rf $TMP/ITS1
export DB=$TMP/alt_primers.sqlite
rm -rf $DB
thapbi_pict import -d $DB -x -k ITS1 --minlen 100 --maxlen 700 \
    --left GAAGGTGAAGTCGTAACAAGG --right AGCGTTCTTCATCGATGTGC \
    -i database/Phytophthora_ITS1_curated.fasta
thapbi_pict prepare-reads -o $TMP -i tests/reads/SRR6303948_sample_*.fastq -a 2 -d $DB
diff $TMP/ITS1/SRR6303948_sample.fasta tests/prepare-reads/SRR6303948_sample_primers.fasta

echo "Testing --flip works"
# Took only 50 reads from original file (about 6%),
# dropping abudnance threshold from 100 to only 10.

export DB=$TMP/nematode_primers.sqlite
rm -rf $DB
thapbi_pict import -d $DB -x -k Nema --minlen 150 --maxlen 350 \
    --left CTGCTGCTGGATCATTACCC --right CGCCAGCACAGCCGTTAG \
    -i database/Phytophthora_ITS1_curated.fasta # empty file would work

rm -rf $TMP/Nema
thapbi_pict prepare-reads -d $DB \
    -i tests/nematodes/sample_R*.fastq.gz -a 10 -o $TMP/
diff $TMP/Nema/sample.fasta tests/nematodes/sample_noflip_a10.fasta  # empty!

rm -rf $TMP/Nema
thapbi_pict prepare-reads --database $DB \
    -i tests/nematodes/sample_R*.fastq.gz -a 10 --flip -o $TMP/
diff $TMP/Nema/sample.fasta tests/nematodes/sample_flip_a10.fasta

echo "Testing pathological trimming example"
rm -rf $TMP/ITS1
thapbi_pict prepare-reads -i tests/reads/6e847180a4da6eed316e1fb98b21218f_R?.fastq \
    -o $TMP/ -a 1
diff $TMP/ITS1/6e847180a4da6eed316e1fb98b21218f.fasta tests/prepare-reads/6e847180a4da6eed316e1fb98b21218f.fasta

echo "$0 - test_prepare-reads.sh passed"
