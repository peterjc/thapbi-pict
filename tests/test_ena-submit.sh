#!/bin/bash

# Copyright 2019-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp/thapbi_pict}/ena_submit
rm -rf $TMP
mkdir -p $TMP

echo "==================="
echo "Checking ena-submit"
echo "==================="
set -x

echo "Checking ena-submit"
thapbi_pict ena-submit 2>&1 | grep "the following arguments are required"
set -o pipefail

thapbi_pict ena-submit --study PRJEB00000 -i tests/reads/ \
    -t tests/reads/metadata.tsv -c 5 -x 1 \
    -o $TMP/ena_submit.tsv -e UTF-8 --flat
diff $TMP/ena_submit.tsv tests/reads/ena_submit.tsv

thapbi_pict ena-submit --study PRJEB00000 -i tests/reads/ \
    -t tests/reads/metadata.tsv -c 5 -x 1 \
    -o $TMP/ena_submit_custom.tsv --metaencoding latin1 \
    --library "Test set" --instrument "Illumina Widget"
diff $TMP/ena_submit_custom.tsv tests/reads/ena_submit_custom.tsv

echo "$0 - test_ena-submit.sh passed"
