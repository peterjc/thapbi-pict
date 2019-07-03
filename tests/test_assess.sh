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

echo "Checking assess"
thapbi_pict assess 2>&1 | grep "the following arguments are required"
set -o pipefail

# Simple examples with expected output to compare against
# Because these are single samples, sseq and useq should agree
for LEVEL in sseq useq; do
    diff tests/assess/ex1.assess.tsv <(thapbi_pict assess -m identity -l $LEVEL -i tests/assess/ex1.known.tsv tests/assess/ex1.identity.tsv -c /dev/null)
    diff tests/assess/ex2.assess.tsv <(thapbi_pict assess -m identity -l $LEVEL -i tests/assess/ex2.known.tsv tests/assess/ex2.identity.tsv -c /dev/null)
    diff tests/assess/ex3.assess.tsv <(thapbi_pict assess -m identity -l $LEVEL -i tests/assess/ex3.known.tsv tests/assess/ex3.identity.tsv -c /dev/null)
    diff tests/assess/ex4.assess.tsv <(thapbi_pict assess -m identity -l $LEVEL -i tests/assess/ex4.known.tsv tests/assess/ex4.identity.tsv -c /dev/null)
    diff tests/assess/unclassified.assess.tsv <(thapbi_pict assess -m identity -l $LEVEL -i tests/assess/unclassified.known.tsv tests/assess/unclassified.identity.tsv -c /dev/null)
    diff tests/assess/fp.assess.tsv <(thapbi_pict assess -m identity -l $LEVEL -i tests/assess/fp.known.tsv tests/assess/fp.identity.tsv -c /dev/null)
done

echo "Checking classifier assessment at different levels"
for LEVEL in sample sseq useq; do
    thapbi_pict assess -i tests/assess/ -o $TMP/assess.tsv -t $TMP/tally.tsv -c $TMP/confusion.tsv -l $LEVEL -m identity
    diff tests/assess/tally_$LEVEL.tsv $TMP/tally.tsv
    diff tests/assess/assess_$LEVEL.tsv $TMP/assess.tsv
    diff tests/assess/confusion_$LEVEL.tsv $TMP/confusion.tsv
done

echo "Checking warning for unexpected species"
set +o pipefail
thapbi_pict assess -m identity --input tests/assess/*.identity.tsv tests/assess/*.known.tsv 2>&1 | grep "Expected species Phytophthora fallax was not a possible prediction"
set -o pipefail


if [ ! -f $TMP/thapbi_swarm/DNAMIX_S95_L001.swarm.tsv ]; then echo "Run test_classify.sh to setup test input"; false; fi
if [ ! -f $TMP/DNAMIX_S95_L001.identity.tsv ]; then echo "Run test_classify.sh to setup test input"; false; fi

rm -rf $TMP/assess_swarm_vs_identity.tsv
rm -rf $TMP/confusion_swarm_vs_identity.tsv

# Don't have a gold standard known truth to test this against, so test swarm vs identity
thapbi_pict assess -i $TMP/thapbi_swarm/DNAMIX_S95_L001.swarm.tsv $TMP/DNAMIX_S95_L001.identity.tsv -m swarm -k identity -o $TMP/assess_swarm_vs_identity.tsv -c $TMP/confusion_swarm_vs_identity.tsv

# Check assessment output to stdout works (default):
thapbi_pict assess -i $TMP/thapbi_swarm/DNAMIX_S95_L001.swarm.tsv $TMP/DNAMIX_S95_L001.identity.tsv -m swarm -k identity > $TMP/stdout.txt
diff $TMP/stdout.txt $TMP/assess_swarm_vs_identity.tsv

# Check confusion matrix output to stdout works, also give one input as a directory name
thapbi_pict assess -i $TMP/thapbi_swarm $TMP/DNAMIX_S95_L001.identity.tsv -m swarm -k identity -o $TMP/assess_swarm_vs_identity.tsv -c - > $TMP/stdout.txt
diff $TMP/stdout.txt $TMP/confusion_swarm_vs_identity.tsv

echo "$0 - test_assess.sh passed"
