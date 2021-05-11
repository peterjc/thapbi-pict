#!/bin/bash

# Copyright 2019-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
set -o pipefail

# Note all tests here (initially) using default database:

export TMP=${TMP:-/tmp/thapbi_pict}/conflicts
rm -rf $TMP
mkdir -p $TMP

echo "=================="
echo "Checking conflicts"
echo "=================="
set -x
set -o pipefail

thapbi_pict conflicts -o $TMP/conflicts.tsv
diff $TMP/conflicts.tsv tests/conflicts/default.tsv

# Same test is used in tests/test_curated-import.sh
# Only 6 FASTA records, but two are double entries so want 8 here
export DB=$TMP/dup_seqs.sqlite
rm -rf $DB
thapbi_pict import -x -d $DB -i tests/curated-import/dup_seqs.fasta -c ncbi -s $'\001' -k ITS1 -l "N" -r "N"
thapbi_pict conflicts -d $TMP/dup_seqs.sqlite -o $TMP/dup_seqs.tsv
diff $TMP/dup_seqs.tsv tests/conflicts/dup_seqs.tsv

echo "$0 - test_dump.sh passed"
