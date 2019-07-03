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

echo "Checking sample-summary"
thapbi_pict sample-summary 2>&1 | grep "the following arguments are required"
thapbi_pict sample-summary -o '' -i tests/classify 2>&1 | grep "No output file specified"
set -o pipefail

# Passing filename, default method, explicit min abundance
rm -rf $TMP/human.txt $TMP/computer.tsv
thapbi_pict sample-summary -m identity -a 1 -r $TMP/human.txt -o $TMP/computer.tsv -i tests/classify/*.identity.tsv
diff $TMP/human.txt tests/sample-summary/classify.identity.txt
diff $TMP/computer.tsv tests/sample-summary/classify.identity.tsv

# Passing a folder, trying different methods
# Skipping swarm based classifiers and they don't work on one of the test cases
for M in identity onebp blast; do
    rm -rf $TMP/human.txt $TMP/computer.tsv
    thapbi_pict sample-summary -m $M -r $TMP/human.txt --output $TMP/computer.tsv --input tests/classify/
    diff $TMP/human.txt tests/sample-summary/classify.$M.txt
    diff $TMP/computer.tsv tests/sample-summary/classify.$M.tsv
    # And again, but with metadata
    rm -rf $TMP/human.txt $TMP/computer.tsv
    thapbi_pict sample-summary -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5 -m $M -r $TMP/human.txt -o $TMP/computer.tsv -i tests/classify/
    diff $TMP/human.txt tests/sample-summary/classify-meta.$M.txt
    # This currently does not include metadata, and sorting does not change as only one sample
    diff $TMP/computer.tsv tests/sample-summary/classify.$M.tsv
done

# More complicated metadata testing
for M in identity; do
    rm -rf $TMP/human.txt $TMP/computer.tsv
    thapbi_pict sample-summary -m $M -r $TMP/human.txt -o $TMP/computer.tsv -i tests/assess/
    diff $TMP/human.txt tests/sample-summary/assess.$M.txt
    diff $TMP/computer.tsv tests/sample-summary/assess.$M.tsv
    # And again, but with metadata
    rm -rf $TMP/human.txt $TMP/computer.tsv
    thapbi_pict sample-summary -m $M -t tests/assess/meta.tsv -x 2 -c 1 -r $TMP/human.txt -o $TMP/computer.tsv -i tests/assess/
    diff $TMP/human.txt tests/sample-summary/assess-meta.$M.txt
    diff $TMP/computer.tsv tests/sample-summary/assess-meta.$M.tsv
done

echo "$0 - test_sample-summary.sh passed"
