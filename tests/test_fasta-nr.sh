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

echo "================="
echo "Checking fasta-nr"
echo "================="
set -x
thapbi_pict fasta-nr 2>&1 | grep "Require at least one of"
set -o pipefail

rm -rf $TMP/all.fasta
thapbi_pict fasta-nr -i tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/all.fasta
# Should be identical bar no header lines, and loss of HMM names in FASTA descriptions
diff $TMP/all.fasta <(grep -v "^#" tests/prepare-reads/DNAMIX_S95_L001.fasta | cut -f 1 -d " ")

echo "$0 - test_fasta-nr.sh passed"
