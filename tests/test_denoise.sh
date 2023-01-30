#!/bin/bash

# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp/thapbi_pict}/sample-tally
    rm -rf $TMP
mkdir -p $TMP

echo "================"
echo "Checking denoise"
echo "================"
set -x
thapbi_pict denoise 2>&1 | grep "the following arguments are required"
set -o pipefail

# These tests are also used with the sample-tally command:
for AFTER in tests/read-correction/*.unoise.fasta; do
    BEFORE=${AFTER%%.*}.before.fasta
    echo "Checking denoising $BEFORE --> $AFTER"
    thapbi_pict denoise -i $BEFORE -o $TMP/after.fasta --minlen 60 -t 0 -α 2.0 -γ 4
    echo diff $TMP/after.fasta $AFTER
    diff $TMP/after.fasta $AFTER

    # TODO - compare with usearch or vsearch?
    # Something like this with a suitable example...
    # $ vsearch --unoise_alpha 2 --minsize 4 --cluster_unoise \
    #     <(python scripts/swarm2usearch.py $BEFORE) \
    #     --centroids $TMP/vsearch.fasta --sizein --sizeout --sizeorder
    # $ diff $TMP/vsearch.fasta \
    #     <(python scripts/swarm2usearch.py $AFTER)
done

echo "$0 - test_denoise.sh passed"
