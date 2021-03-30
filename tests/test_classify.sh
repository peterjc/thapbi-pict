#!/bin/bash

# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "================="
echo "Checking classify"
echo "================="
set -x
thapbi_pict classify 2>&1 | grep "the following arguments are required"
thapbi_pict classify -d "sqlite:///:memory:" -i hypothetical_example.fasta 2>&1 | grep "cannot classify anything"
set -o pipefail

export DB=$TMP/curated.sqlite
if [ ! -f $DB ]; then echo "Run test_curated-import.sh to setup test DB"; false; fi

# Passing one filename; default output dir:
rm -rf $TMP/classify/
mkdir -p $TMP/classify/
cp database/Phytophthora_ITS1_curated.fasta $TMP/classify/
thapbi_pict classify -m identity -d $DB -i $TMP/classify/Phytophthora_ITS1_curated.fasta
if [ "`grep -c -v '^#' $TMP/classify/Phytophthora_ITS1_curated.identity.tsv`" -ne "`grep -c '^>' $TMP/classify/Phytophthora_ITS1_curated.fasta`" ]; then echo "Expected one line per input seq"; false; fi

rm -rf $TMP/DNAMIX_S95_L001.identity.tsv
rm -rf $TMP/thapbi_onebp
rm -rf $TMP/thapbi_swarm
rm -rf $TMP/thapbi_blast
rm -rf $TMP/thapbi_1s3g
mkdir -p $TMP/thapbi_onebp
mkdir -p $TMP/thapbi_swarm
mkdir -p $TMP/thapbi_blast
mkdir -p $TMP/thapbi_1s3g

# Explicitly setting output directory, would be here anyway:
thapbi_pict classify -m identity -d $DB -i tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/
thapbi_pict classify -m onebp -d $DB -i tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/thapbi_onebp
thapbi_pict classify -m blast -d $DB -i tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/thapbi_blast
thapbi_pict classify -m swarm -d $DB -i tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/thapbi_swarm
thapbi_pict classify -m swarmid -d $DB -i tests/prepare-reads/DNAMIX_S95_L001.fasta
thapbi_pict classify -m 1s3g -d $DB -i tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/1s3g

# Passing one directory name (should get all 2 FASTA files):
rm -rf $TMP/duo
mkdir -p $TMP/duo
cp database/Phytophthora_ITS1_curated.fasta $TMP/duo/
cp database/controls.fasta $TMP/duo/
thapbi_pict classify -m identity -d $DB -i $TMP/duo -o $TMP/duo
if [ "`ls -1 $TMP/duo/*.identity.tsv | wc -l`" -ne "2" ]; then echo "Expected 4 files"; false; fi

# Test using sequences from a single isolate control,
rm -rf $TMP/P-infestans-T30-4.*.tsv
for M in identity onebp substr blast swarm swarmid 1s3g; do
    echo "Checking single isolate control with $M"
    thapbi_pict classify -d $DB -i tests/classify/P-infestans-T30-4.fasta -o $TMP/ -m $M
    diff $TMP/P-infestans-T30-4.$M.tsv tests/classify/P-infestans-T30-4.$M.tsv
done

rm -rf $TMP/multiple_its1.*.tsv
# Have not handled this in swarm classifier....
for M in identity onebp blast 1s3g; do
    echo "Checking multiple HMM containing sequences with $M"
    thapbi_pict classify -i tests/classify/multiple_its1.fasta -o $TMP/ -m $M
    diff $TMP/multiple_its1.$M.tsv tests/classify/multiple_its1.$M.tsv
done

rm -rf $TMP/hmm_trim.*.tsv
# Swarm classifier can't cope with multiple HMM hits...
for M in identity onebp blast 1s3g; do
    echo "Checking HMM trim corner cases with $M"
    thapbi_pict classify -d $DB -i tests/classify/hmm_trim.fasta -o $TMP/ -m $M -a 50
    diff $TMP/hmm_trim.$M.tsv tests/classify/hmm_trim.$M.tsv
done

echo "$0 - test_classify.sh passed"
