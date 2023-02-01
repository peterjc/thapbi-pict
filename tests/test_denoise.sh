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
    echo "Checking denoising $BEFORE --> $AFTER (UNOISE)"
    thapbi_pict denoise -i $BEFORE -o $TMP/after.fasta \
                --denoise unoise --minlen 60 -t 0 -α 2.0 -γ 4
    echo diff $TMP/after.fasta $AFTER
    diff $TMP/after.fasta $AFTER
done

if ! [ -x "$(command -v vsearch)" ]; then
    echo "Skipping testing using VSEARCH"
else
    for AFTER in tests/read-correction/*.vsearch.fasta; do
        BEFORE=${AFTER%%.*}.before.fasta
        echo "Checking denoising $BEFORE --> $AFTER (VSEARCH)"
        thapbi_pict denoise -i $BEFORE -o $TMP/after.fasta \
                --denoise vsearch --minlen 60 -t 0 -α 2.0 -γ 4
        echo diff $TMP/after.fasta $AFTER
        diff $TMP/after.fasta $AFTER

        echo "Checking versus direct use of vsearch"
        python scripts/swarm2usearch.py $BEFORE > $TMP/before.fasta
        vsearch --unoise_alpha 2.0 --minsize 4 --sizein --sizeout \
            --cluster_unoise $TMP/before.fasta \
            --centroids $TMP/vsearch.fasta
        # The THAPBI-PICT output is upper case and not line wrapped,
        # the VSEARCH output is mixed case (masking?) and wrapped.
        # Also while both are sorted by abundance, tie breaking differs.
        # So, can just compare the title lines (MD5 and abundance)
        diff <(grep "^>" $TMP/vsearch.fasta | sort) \
             <(python scripts/swarm2usearch.py $AFTER | grep "^>" | sort)
    done
fi

echo "$0 - test_denoise.sh passed"
