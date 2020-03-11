#!/bin/bash

# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "====================="
echo "Checking read-summary"
echo "====================="
set -x
thapbi_pict summary 2>&1 | grep "the following arguments are required"
thapbi_pict summary -o '' -i tests/classify 2>&1 | grep "Output directory name blank"
set -o pipefail

thapbi_pict summary -i tests/prepare-reads/DNAMIX_S95_L001.fasta \
	    $TMP/DNAMIX_S95_L001.identity.tsv \
	    -m identity -o $TMP/ -r summary

# With metadata, using default method, -m onebp
thapbi_pict summary --input tests/classify/P-infestans-T30-4.fasta tests/classify/P-infestans-T30-4.onebp.tsv -o $TMP/ -r summary -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5
diff $TMP/summary.reads.onebp.tsv tests/classify/P-infestans-T30-4.summary.tsv

# Now require metadata, but give entire folder as input
thapbi_pict summary --input tests/classify/ -o $TMP/ -r summary -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5 -r summary -q
diff $TMP/summary.reads.onebp.tsv tests/classify/P-infestans-T30-4.summary.tsv


echo "$0 - test_read-summary.sh passed"
