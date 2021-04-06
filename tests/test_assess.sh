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

echo "==============="
echo "Checking assess"
echo "==============="
set -x
thapbi_pict assess 2>&1 | grep "the following arguments are required"
set -o pipefail

export DB=$TMP/seven.sqlite
rm -rf $DB
thapbi_pict curated-import -d $DB -i tests/assess/seven.fasta -x

# Simple examples with expected output to compare against
for SAMPLE in ex1 ex2 ex3 ex4 unclassified fp; do
    thapbi_pict assess -m identity -i tests/assess/$SAMPLE.known.tsv tests/assess/$SAMPLE.identity.tsv -o $TMP/$SAMPLE.assess.tsv -d $DB
    diff $TMP/$SAMPLE.assess.tsv tests/assess/$SAMPLE.assess.tsv
done

echo "Checking all classifier assessment outputs"
thapbi_pict assess -i tests/assess/ -o $TMP/assess.tsv -t $TMP/tally.tsv -c $TMP/confusion.tsv -m identity -d $DB
diff $TMP/tally.tsv tests/assess/samples.tally.tsv
diff $TMP/assess.tsv tests/assess/samples.assess.tsv
diff $TMP/confusion.tsv tests/assess/samples.confusion.tsv

echo "Checking warning for unexpected species"
set +o pipefail
thapbi_pict assess -d $DB -m identity --input tests/assess/*.identity.tsv tests/assess/*.known.tsv 2>&1 | grep "WARNING: 1 expected species were not a possible prediction: Phytophthora fallax"
set -o pipefail

if [ ! -f $TMP/thapbi_swarm/DNAMIX_S95_L001.swarm.tsv ]; then echo "Run test_classify.sh to setup test input"; false; fi
if [ ! -f $TMP/DNAMIX_S95_L001.identity.tsv ]; then echo "Run test_classify.sh to setup test input"; false; fi

rm -rf $TMP/assess_swarm_vs_identity.tsv
rm -rf $TMP/confusion_swarm_vs_identity.tsv

echo "Testing swarm vs identity"
# Note using the default DB here, as would be used in test_classify.sh to generate inputs
# Don't have a gold standard known truth to test this against, so test swarm vs identity
thapbi_pict assess -i $TMP/thapbi_swarm/DNAMIX_S95_L001.swarm.tsv $TMP/DNAMIX_S95_L001.identity.tsv -m swarm -k identity -o $TMP/assess_swarm_vs_identity.tsv -c $TMP/confusion_swarm_vs_identity.tsv

# Check assessment output to stdout works (default):
thapbi_pict assess -i $TMP/thapbi_swarm/DNAMIX_S95_L001.swarm.tsv $TMP/DNAMIX_S95_L001.identity.tsv -m swarm -k identity > $TMP/stdout.txt
diff $TMP/stdout.txt $TMP/assess_swarm_vs_identity.tsv

# Check confusion matrix output to stdout works, also give one input as a directory name
thapbi_pict assess -i $TMP/thapbi_swarm $TMP/DNAMIX_S95_L001.identity.tsv -m swarm -k identity -o $TMP/assess_swarm_vs_identity.tsv -c - > $TMP/stdout.txt
diff $TMP/stdout.txt $TMP/confusion_swarm_vs_identity.tsv

echo "$0 - test_assess.sh passed"
