#!/bin/bash

# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eux
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking read-summary"
thapbi_pict read-summary 2>&1 | grep "the following arguments are required"
thapbi_pict read-summary -o '' -i tests/classify 2>&1 | grep "No output file specified"
set -o pipefail

# One method:
thapbi_pict read-summary -i tests/prepare-reads/DNAMIX_S95_L001.fasta $TMP/DNAMIX_S95_L001.identity.tsv -m identity -o $TMP/read-summary_identity.tsv -e $TMP/read-summary_swarm_and_blast.xlsx

# Two methods:
thapbi_pict read-summary -i tests/prepare-reads/DNAMIX_S95_L001.fasta  $TMP/thapbi_swarm/ $TMP/thapbi_blast/ -m blast,swarm -o $TMP/read-summary_swarm_and_blast.tsv

# With metadata, using default method, -m onebp
thapbi_pict read-summary --input tests/classify/P-infestans-T30-4.fasta tests/classify/P-infestans-T30-4.onebp.tsv -o $TMP/read-summary_onebp.tsv -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5 -e $TMP/read-summary_onebp.xlsx
diff $TMP/read-summary_onebp.tsv tests/classify/P-infestans-T30-4.summary.tsv

echo "$0 - test_read-summary.sh passed"
