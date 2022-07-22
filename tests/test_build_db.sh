#!/bin/bash

# Copyright 2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -euo pipefail

export TMP=${TMP:-/tmp/thapbi_pict}/build_db
rm -rf $TMP
mkdir -p $TMP

echo "==================================="
echo "Dumping current DB as TSV and FASTA"
echo "==================================="
set -x

thapbi_pict dump -f fasta -o "$TMP/old.fasta"
thapbi_pict dump -m -o "$TMP/old.txt"
thapbi_pict dump -o "$TMP/old.tsv"

set +x
echo "========================="
echo "Rebuilding the default DB"
echo "========================="
set -x

cd database/
rm -rf ITS1_DB.sqlite ITS1_DB.sql ITS1_DB.txt ITS1_DB.fasta
./build_ITS1_DB.sh
cd ..

set +x
echo "========================================"
echo "Checking rebuilt DB matches expectations"
echo "========================================"
set -x

diff "$TMP/old.fasta" database/ITS1_DB.fasta
diff "$TMP/old.txt" database/ITS1_DB.txt
diff "$TMP/old.tsv" database/ITS1_DB.tsv
# git diff database/ITS1_DB.sql

echo "$0 - test_build_db.sh passed"
