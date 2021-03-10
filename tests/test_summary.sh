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

echo "================"
echo "Checking summary"
echo "================"
set -x
thapbi_pict summary 2>&1 | grep "the following arguments are required"
thapbi_pict summary -o '' -i tests/classify 2>&1 | grep "Output directory name blank"
set -o pipefail

thapbi_pict summary -i tests/summary_meta/ -m 1s3g -o $TMP/ \
    -t tests/summary_meta/metadata.tsv -x 4 -c 1,2 -r summary
diff $TMP/summary.reads.1s3g.tsv tests/summary_meta/summary.reads.1s3g.tsv
diff $TMP/summary.samples.1s3g.tsv tests/summary_meta/summary.samples.1s3g.tsv
diff $TMP/summary.samples.1s3g.txt tests/summary_meta/summary.samples.1s3g.txt

thapbi_pict summary -i tests/summary_meta/ -m 1s3g -o $TMP/ \
    -t tests/summary_meta/metadata.tsv -x 4 -c 1,2 -u -r summary-u
diff $TMP/summary-u.reads.1s3g.tsv tests/summary_meta/summary.reads.1s3g.tsv  # no change
diff $TMP/summary-u.samples.1s3g.tsv tests/summary_meta/summary-u.samples.1s3g.tsv
diff $TMP/summary-u.samples.1s3g.txt tests/summary_meta/summary-u.samples.1s3g.txt

thapbi_pict summary -i tests/summary_meta/ -m 1s3g -o $TMP/ \
    -t tests/summary_meta/metadata.tsv -x 4 -c 1,2 -q -r summary-q
diff $TMP/summary-q.reads.1s3g.tsv tests/summary_meta/summary-q.reads.1s3g.tsv
diff $TMP/summary-q.samples.1s3g.tsv tests/summary_meta/summary-q.samples.1s3g.tsv
diff $TMP/summary-q.samples.1s3g.txt tests/summary_meta/summary-q.samples.1s3g.txt

thapbi_pict summary -i tests/summary_meta/ -m 1s3g -o $TMP/ \
    -t tests/summary_meta/metadata.tsv -x 4 -c 1,2 -q -u -r summary-qu
diff $TMP/summary-qu.reads.1s3g.tsv tests/summary_meta/summary-q.reads.1s3g.tsv  # no change
diff $TMP/summary-qu.samples.1s3g.tsv tests/summary_meta/summary-qu.samples.1s3g.tsv
diff $TMP/summary-qu.samples.1s3g.txt tests/summary_meta/summary-qu.samples.1s3g.txt


thapbi_pict summary -i tests/prepare-reads/DNAMIX_S95_L001.fasta \
    $TMP/DNAMIX_S95_L001.identity.tsv \
    -m identity -o $TMP/ -r summary

# With metadata, using default method, -m onebp
thapbi_pict summary --input tests/classify/P-infestans-T30-4.fasta tests/classify/P-infestans-T30-4.onebp.tsv -o $TMP/ -r summary -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5
diff $TMP/summary.reads.onebp.tsv tests/classify/P-infestans-T30-4.summary.tsv

# Now require metadata, but give entire folder as input
thapbi_pict summary --input tests/classify/ -o $TMP/ -r summary -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5 -r summary -q
diff $TMP/summary.reads.onebp.tsv tests/classify/P-infestans-T30-4.summary.tsv


# Passing filename, default method, explicit min abundance
rm -rf $TMP/human.txt $TMP/test-case.tsv $TMP/test-case.xlsx
thapbi_pict summary -m identity -a 99 -o $TMP/ -r test-case \
    -i tests/classify/*.fasta tests/classify/*.identity.tsv
diff $TMP/test-case.samples.identity.txt tests/summary/classify.identity.txt
diff $TMP/test-case.samples.identity.tsv tests/summary/classify.identity.tsv

# Passing a folder, trying different methods
for M in identity onebp blast; do
    rm -rf $TMP/test-case.samples.$M.txt $TMP/test-case.samples.$M.tsv
    thapbi_pict summary -m $M -r test-case --output $TMP/ --input tests/classify/
    diff $TMP/test-case.samples.$M.txt tests/summary/classify.$M.txt
    diff $TMP/test-case.samples.$M.tsv tests/summary/classify.$M.tsv
    # And again, but with metadata
    rm -rf $TMP/test-case.samples.$M.txt $TMP/test-case.samples.$M.tsv
    thapbi_pict summary -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5 -m $M -o $TMP -r test-case -i tests/classify/
    diff $TMP/test-case.samples.$M.txt tests/summary/classify-meta.$M.txt
    diff $TMP/test-case.samples.$M.tsv tests/summary/classify-meta.$M.tsv
    # Now require metadata...
    rm -rf $TMP/test-case.samples.$M.txt $TMP/test-case.samples.$M.tsv
    thapbi_pict summary -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5 -m $M -o $TMP/ -r test-case -i tests/classify/ -q
    diff $TMP/test-case.samples.$M.txt tests/summary/classify-meta-req.$M.txt
    diff $TMP/test-case.samples.$M.tsv tests/summary/classify-meta-req.$M.tsv
done

echo "$0 - test_summary.sh passed"
