#!/bin/bash

# Copyright 2019-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp/thapbi_pict}/summary
rm -rf $TMP
mkdir -p $TMP

echo "================"
echo "Checking summary"
echo "================"
set -x
thapbi_pict summary 2>&1 | grep "the following arguments are required"
thapbi_pict summary -o '' -i tests/classify 2>&1 | grep "Output stem is blank"
set -o pipefail

thapbi_pict summary -i tests/summary_meta/ -m 1s3g -o $TMP/summary \
    -t tests/summary_meta/metadata.tsv -x 4 -c 1,2
diff $TMP/summary.reads.1s3g.tsv tests/summary_meta/summary.reads.1s3g.tsv
diff $TMP/summary.samples.1s3g.tsv tests/summary_meta/summary.samples.1s3g.tsv
diff $TMP/summary.samples.1s3g.txt tests/summary_meta/summary.samples.1s3g.txt

thapbi_pict summary -i tests/summary_meta/ -m 1s3g -o $TMP/summary-u \
    -t tests/summary_meta/metadata.tsv -x 4 -c 1,2 -u
diff $TMP/summary-u.reads.1s3g.tsv tests/summary_meta/summary.reads.1s3g.tsv  # no change
diff $TMP/summary-u.samples.1s3g.tsv tests/summary_meta/summary-u.samples.1s3g.tsv
diff $TMP/summary-u.samples.1s3g.txt tests/summary_meta/summary-u.samples.1s3g.txt

thapbi_pict summary -i tests/summary_meta/ -m 1s3g -o $TMP/summary-q \
    -t tests/summary_meta/metadata.tsv -x 4 -c 1,2 -q
diff $TMP/summary-q.reads.1s3g.tsv tests/summary_meta/summary-q.reads.1s3g.tsv
diff $TMP/summary-q.samples.1s3g.tsv tests/summary_meta/summary-q.samples.1s3g.tsv
diff $TMP/summary-q.samples.1s3g.txt tests/summary_meta/summary-q.samples.1s3g.txt

thapbi_pict summary -i tests/summary_meta/ -m 1s3g -o $TMP/summary-qu \
    -t tests/summary_meta/metadata.tsv -x 4 -c 1,2 -q -u
diff $TMP/summary-qu.reads.1s3g.tsv tests/summary_meta/summary-q.reads.1s3g.tsv  # no change
diff $TMP/summary-qu.samples.1s3g.tsv tests/summary_meta/summary-qu.samples.1s3g.tsv
diff $TMP/summary-qu.samples.1s3g.txt tests/summary_meta/summary-qu.samples.1s3g.txt

# This was originally created in test_classify.sh
if [ ! -f $TMP/DNAMIX_S95_L001.identity.tsv ]; then
    thapbi_pict classify -m identity -i tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/
fi
thapbi_pict summary -i tests/prepare-reads/DNAMIX_S95_L001.fasta \
    $TMP/DNAMIX_S95_L001.identity.tsv -m identity -o $TMP/summary

# With metadata, using default method, -m onebp
thapbi_pict summary \
    --input tests/classify/P-infestans-T30-4.fasta tests/classify/P-infestans-T30-4.onebp.tsv \
    -o $TMP/summary -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5
diff $TMP/summary.reads.onebp.tsv tests/classify/P-infestans-T30-4.summary.tsv

# Now require metadata, but give entire folder as input
thapbi_pict summary --input tests/classify/ -o $TMP/summary -q \
    -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5 -e latin1
diff $TMP/summary.reads.onebp.tsv tests/classify/P-infestans-T30-4.summary.tsv


# Passing filename, default method, explicit min abundance
# Note:
#
#    $ grep threshold tests/classify/*.fasta
#    tests/classify/P-infestans-T30-4.fasta:#threshold:500
#    tests/classify/hmm_trim.fasta:#threshold:50
#
# i.e. the default of 100 or an explicit value more than 50
# could exclude some reads - the file is not abundance sorted
# but the lowest and only values under 100 are:
#
#    >62e52211a6661fb05ae292808d79a4a3_75
#    >330d9e67a26944344219464449fed619_68
#
rm -rf $TMP/human.txt $TMP/test-case.tsv $TMP/test-case.xlsx
thapbi_pict summary -m identity -a 99 -o $TMP/test-case \
    -i tests/classify/*.fasta tests/classify/*.identity.tsv
diff $TMP/test-case.samples.identity.txt tests/summary/classify.identity.txt
diff $TMP/test-case.samples.identity.tsv tests/summary/classify.identity.tsv

# Passing a folder, trying different methods
for M in identity onebp blast; do
    rm -rf $TMP/test-case.samples.$M.txt $TMP/test-case.samples.$M.tsv
    thapbi_pict summary -m $M -a 99 --output $TMP/test-case --input tests/classify/
    diff $TMP/test-case.samples.$M.txt tests/summary/classify.$M.txt
    diff $TMP/test-case.samples.$M.tsv tests/summary/classify.$M.tsv
    # And again, but with metadata
    rm -rf $TMP/test-case.samples.$M.txt $TMP/test-case.samples.$M.tsv
    thapbi_pict summary -m $M -a 99 -o $TMP/test-case -i tests/classify/ \
        -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5
    diff $TMP/test-case.samples.$M.txt tests/summary/classify-meta.$M.txt
    diff $TMP/test-case.samples.$M.tsv tests/summary/classify-meta.$M.tsv
    # Now require metadata...
    rm -rf $TMP/test-case.samples.$M.txt $TMP/test-case.samples.$M.tsv
    thapbi_pict summary -m $M -a 99 -o $TMP/test-case -i tests/classify/ -q \
        -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5
    diff $TMP/test-case.samples.$M.txt tests/summary/classify-meta-req.$M.txt
    diff $TMP/test-case.samples.$M.tsv tests/summary/classify-meta-req.$M.tsv
done

echo "$0 - test_summary.sh passed"
