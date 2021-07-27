#!/bin/bash

# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp/thapbi_pict}/classify
rm -rf $TMP
mkdir -p $TMP

echo "================="
echo "Checking classify"
echo "================="
set -x
thapbi_pict classify 2>&1 | grep "the following arguments are required"
thapbi_pict classify -d "sqlite:///:memory:" -i hypothetical_example.fasta 2>&1 | grep "cannot classify anything"
set -o pipefail

# The same DB is also created and used in tests/test_curated-import.sh
export DB=$TMP/curated.sqlite
if [ ! -f $DB ]; then
    thapbi_pict load-tax -d $DB -t new_taxdump_2019-09-01
    thapbi_pict import -d $DB -i database/Phytophthora_ITS1_curated.fasta -s ";" \
                -k ITS1 -l GAAGGTGAAGTCGTAACAAGG -r GCARRGACTTTCGTCCCYRC
fi

# Passing one filename; default output dir:
rm -rf $TMP/input/
mkdir -p $TMP/input/
cp database/Phytophthora_ITS1_curated.fasta $TMP/input/
thapbi_pict classify -m identity -d $DB -i $TMP/input/Phytophthora_ITS1_curated.fasta
if [ "`grep -c -v '^#' $TMP/input/Phytophthora_ITS1_curated.identity.tsv`" -ne "`grep -c '^>' $TMP/input/Phytophthora_ITS1_curated.fasta`" ]; then echo "Expected one line per input seq"; false; fi

rm -rf $TMP/DNAMIX_S95_L001.identity.tsv
rm -rf $TMP/thapbi_onebp
rm -rf $TMP/thapbi_blast
rm -rf $TMP/thapbi_1s3g
mkdir -p $TMP/thapbi_onebp
mkdir -p $TMP/thapbi_blast
mkdir -p $TMP/thapbi_1s3g

# Explicitly setting output directory, would be here anyway:
thapbi_pict classify -m identity -d $DB -i tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/
thapbi_pict classify -m onebp -d $DB -i tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/thapbi_onebp
thapbi_pict classify -m blast -d $DB -i tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/thapbi_blast
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
for M in identity onebp substr blast 1s3g; do
    echo "Checking single isolate control with $M"
    thapbi_pict classify -d $DB -i tests/classify/P-infestans-T30-4.fasta -o $TMP/ -m $M
    diff $TMP/P-infestans-T30-4.$M.tsv tests/classify/P-infestans-T30-4.$M.tsv
done

rm -rf $TMP/hmm_trim.*.tsv
for M in identity onebp blast 1s3g; do
    echo "Checking HMM trim corner cases with $M"
    thapbi_pict classify -d $DB -i tests/classify/hmm_trim.fasta -o $TMP/ -m $M -a 50
    diff $TMP/hmm_trim.$M.tsv tests/classify/hmm_trim.$M.tsv
done

for M in identity onebp blast 1s2g 1s3g 1s4g 1s5g; do
    # Using default DB
    echo "Checking genus corner cases with $M"
    thapbi_pict classify -i tests/classifier/corner_cases_query.fasta -o $TMP/ -m $M -k ITS1
    diff $TMP/corner_cases_query.$M.tsv tests/classifier/corner_cases_query.$M.tsv
done

for EXAMPLE in P_bilorbang P_vulcanica genus_boundary; do

    DB=$TMP/${EXAMPLE}.sqlite
    rm -rf $DB $TMP/${EXAMPLE}_query.*
    set -x
    thapbi_pict import -d $DB -i tests/classifier/${EXAMPLE}.fasta -x -k ITS1 -l NNN -r NNN
    diff tests/classifier/${EXAMPLE}.fasta <(thapbi_pict dump -m -f fasta -d $DB)
    set +x
    for M in identity onebp 1s2g 1s3g 1s4g 1s5g; do
        echo
        echo "Checking ${EXAMPLE} example with $M classifier"
        set -x
        thapbi_pict classify -d $DB -i tests/classifier/${EXAMPLE}_query.fasta -m $M -o $TMP
        diff $TMP/${EXAMPLE}_query.$M.tsv tests/classifier/${EXAMPLE}_query.$M.tsv
        set +x
    done
done

echo "$0 - test_classify.sh passed"
