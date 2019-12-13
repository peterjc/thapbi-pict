#!/bin/bash

# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
set -o pipefail

export TMP=${TMP:-/tmp}

# Note all tests here (initially) using default database:

echo "=================="
echo "Checking conflicts"
echo "=================="
set -x
set -o pipefail

thapbi_pict conflicts -o $TMP/conflicts.tsv
diff $TMP/conflicts.tsv tests/conflicts/default.tsv

if [ ! -f $TMP/dup_seqs.sqlite ]; then echo "Run tests/test_curated-import.sh to setup test DB"; false; fi
thapbi_pict conflicts -d $TMP/dup_seqs.sqlite -o $TMP/dup_seqs.tsv
diff $TMP/dup_seqs.tsv tests/conflicts/dup_seqs.tsv

echo "$0 - test_dump.sh passed"
