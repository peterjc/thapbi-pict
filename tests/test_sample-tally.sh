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

echo "====================="
echo "Checking sample-tally"
echo "====================="
set -x
thapbi_pict sample-tally 2>&1 | grep "the following arguments are required"
set -o pipefail

rm -rf $TMP/x.tsv
thapbi_pict sample-tally -i tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/x.tsv
diff $TMP/x.tsv tests/sample-tally/DNAMIX_S95_L001.tally.tsv

if ! [ -x "$(command -v biom)" ]; then
    echo 'WARNING: biom-format not installed, skipping some tests'
else
    for X in $TMP/x.tsv; do
        echo "Checking biom convert works with $X"
        biom convert -i $X -o $TMP/biom.tsv --table-type="OTU table" --to-tsv \
             --header-key Sequence --tsv-metadata-formatter naive
        biom convert -i $X -o $TMP/biom.hdf5 --table-type="OTU table" --to-hdf5
        biom convert -i $X -o $TMP/biom.json --table-type="OTU table" --to-json
        # Don't know how stable the output is, so don't use diff
        grep -c "Biological Observation Matrix 1.0.0" $TMP/biom.json

    done
fi

echo "----------------"
echo "Checking denoise"
echo "----------------"

thapbi_pict sample-tally -i tests/read-correction/Rhabditis.before.fasta \
    -o $TMP/Rhabditis.tally.tsv --fasta $TMP/Rhabditis.after.fasta \
    -a 0 -f 0 # No denoise, should be unchangeed
diff $TMP/Rhabditis.after.fasta tests/read-correction/Rhabditis.before.fasta

thapbi_pict sample-tally -i tests/read-correction/Rhabditis.before.fasta \
    -o $TMP/Rhabditis.tally.tsv --fasta $TMP/Rhabditis.after.fasta \
    -a 0 -f 0 --denoise
diff $TMP/Rhabditis.after.fasta tests/read-correction/Rhabditis.after.fasta
thapbi_pict edit-graph -i $TMP/Rhabditis.tally.tsv -f matrix -d '' | cut -f 1-3

echo "$0 - test_sample-tally.sh passed"
