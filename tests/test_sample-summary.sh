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

echo "======================="
echo "Checking sample-summary"
echo "======================="
set -x
thapbi_pict sample-summary 2>&1 | grep "the following arguments are required"
thapbi_pict sample-summary -o '' -i tests/classify 2>&1 | grep "Output directory does not exist"
set -o pipefail

# Passing filename, default method, explicit min abundance
rm -rf $TMP/human.txt $TMP/test-case.tsv $TMP/test-case.xlsx
thapbi_pict sample-summary -m identity -a 99 -o $TMP/ -r test-case -i tests/classify/*.identity.tsv
diff $TMP/test-case.samples.identity.txt tests/sample-summary/classify.identity.txt
diff $TMP/test-case.samples.identity.tsv tests/sample-summary/classify.identity.tsv

# Passing a folder, trying different methods
# Skipping swarm based classifiers and they don't work on one of the test cases
for M in identity onebp blast swarmid swarm; do
    rm -rf $TMP/test-case.samples.$M.txt $TMP/test-case.samples.$M.tsv
    thapbi_pict sample-summary -m $M -r test-case --output $TMP/ --input tests/classify/
    diff $TMP/test-case.samples.$M.txt tests/sample-summary/classify.$M.txt
    diff $TMP/test-case.samples.$M.tsv tests/sample-summary/classify.$M.tsv
    # And again, but with metadata
    rm -rf $TMP/test-case.samples.$M.txt $TMP/test-case.samples.$M.tsv
    thapbi_pict sample-summary -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5 -m $M -o $TMP -r test-case -i tests/classify/
    diff $TMP/test-case.samples.$M.txt tests/sample-summary/classify-meta.$M.txt
    diff $TMP/test-case.samples.$M.tsv tests/sample-summary/classify-meta.$M.tsv
    # Now require metadata...
    rm -rf $TMP/test-case.samples.$M.txt $TMP/test-case.samples.$M.tsv
    thapbi_pict sample-summary -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5 -m $M -o $TMP/ -r test-case -i tests/classify/ -q
    diff $TMP/test-case.samples.$M.txt tests/sample-summary/classify-meta-req.$M.txt
    diff $TMP/test-case.samples.$M.tsv tests/sample-summary/classify-meta-req.$M.tsv
done

# More complicated metadata testing
for M in identity; do
    rm -rf $TMP/test-case.samples.$M.txt $TMP/test-case.samples.$M.tsv $TMP/test-case.xlsx
    thapbi_pict sample-summary -m $M -o $TMP/ -r test-case -i tests/assess/
    diff $TMP/test-case.samples.$M.txt tests/sample-summary/assess.$M.txt
    diff $TMP/test-case.samples.$M.tsv tests/sample-summary/assess.$M.tsv
    # And again, but with metadata
    rm -rf $TMP/test-case.samples.$M.txt $TMP/test-case.samples.$M.tsv $TMP/test-case.xlsx
    thapbi_pict sample-summary -m $M -t tests/assess/meta.tsv -x 2 -c 1 -o $TMP/ -r test-case -i tests/assess/
    diff $TMP/test-case.samples.$M.txt tests/sample-summary/assess-meta.$M.txt
    diff $TMP/test-case.samples.$M.tsv tests/sample-summary/assess-meta.$M.tsv
    # And now requiring the metadata (no change):
    rm -rf $TMP/test-case.samples.$M.txt $TMP/test-case.samples.$M.tsv $TMP/test-case.xlsx
    thapbi_pict sample-summary -m $M -t tests/assess/meta.tsv -x 2 -c 1 -o $TMP/ -r test-case -i tests/assess/ -q
    diff $TMP/test-case.samples.$M.txt tests/sample-summary/assess-meta.$M.txt
    diff $TMP/test-case.samples.$M.tsv tests/sample-summary/assess-meta.$M.tsv
done

echo "$0 - test_sample-summary.sh passed"
