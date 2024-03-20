#!/bin/bash

# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp/thapbi_pict}/sintax_import
rm -rf $TMP
mkdir -p $TMP

export LEFT=GAAGGTGAAGTCGTAACAAGG
export RIGHT=GCARRGACTTTCGTCCCYRC
export RIGHT_RC=GYRGGGACGAAAGTCYYTGC

echo "======================"
echo "Checking sintax-import"
echo "======================"
set -x
thapbi_pict import 2>&1 | grep "the following arguments are required"
# Cannot use validation without having some taxonomy entries
thapbi_pict import -c sintax -d sqlite:///:memory: --input tests/sintax-import/refs.fasta 2>&1 | grep "Taxonomy table empty"
set -o pipefail

# Check hybrid like "Phytophthora humicola x Phytophthora inundata"
# imports as genus="Phytophthora", species="humicola x inundata"
export DB=$TMP/refs.sqlite
rm -rf $DB
thapbi_pict import -x -k ITS1 -l $LEFT -r $RIGHT -c sintax -d $DB -i tests/sintax-import/refs.fasta

if [ "$(sqlite3 "$DB" "SELECT COUNT(id) FROM data_source;")" -ne "1" ]; then
    echo "Wrong data_source count"
    false
fi
if [ "$(sqlite3 "$DB" "SELECT COUNT(id) FROM sequence_source;")" -ne "2" ]; then
    echo "Wrong sequence_source count"
    false
fi
if [ "$(sqlite3 "$DB" "SELECT COUNT(id) FROM marker_sequence;")" -ne "2" ]; then
    echo "Wrong marker_sequence count"
    false
fi
if [ "$(sqlite3 "$DB" "SELECT COUNT(id) FROM taxonomy;")" -ne "2" ]; then
    echo "Wrong taxonomy count"
    false
fi
if [ "$(sqlite3 "$DB" "SELECT MAX(LENGTH(sequence)) FROM marker_sequence;")" -ne "226" ]; then
    echo "Wrong max ITS1 sequence length"
    false
fi

thapbi_pict dump -d $DB -m -o $TMP/refs.tsv
diff $TMP/refs.tsv tests/sintax-import/refs.tsv

echo "$0 - test_sintax-import.sh passed"
