#!/bin/bash

# Copyright 2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -euo pipefail

export TMP=${TMP:-/tmp/thapbi_pict}/classify
mkdir -p $TMP

echo "==================="
echo "Checking pooling.py"
echo "==================="

scripts/pooling.py -i tests/pooling/example.samples.onebp.tsv -c 2,3,4,5 --pcr -z -p 6 -o $TMP/example
diff $TMP/example.tsv tests/pooling/example.pooled.tsv

echo "$0 - test_pooling.sh passed"
