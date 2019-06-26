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

echo "Checking pipeline"
thapbi_pict pipeline 2>&1 | grep "the following arguments are required"
set -o pipefail

rm -rf $TMP/intermediate $TMP/output
mkdir $TMP/intermediate $TMP/output
thapbi_pict pipeline -s $TMP/intermediate -o $TMP/output -i tests/reads/
diff $TMP/intermediate/DNAMIX_S95_L001.fasta tests/prepare-reads/DNAMIX_S95_L001.fasta
diff $TMP/output/thapbi-pict.samples.txt tests/pipeline/thapbi-pict.samples.onebp.txt
diff $TMP/output/thapbi-pict.samples.tsv tests/pipeline/thapbi-pict.samples.onebp.tsv
diff $TMP/output/thapbi-pict.reads.tsv tests/pipeline/thapbi-pict.reads.onebp.tsv
diff $TMP/output/thapbi-pict.edit-graph.xgmml tests/pipeline/thapbi-pict.edit-graph.xgmml

# Leaving the intermediate files in place... plus some stray files:
touch $TMP/intermediate/unwanted.fasta
touch $TMP/intermediate/unwanted.onebp.tsv
touch $TMP/intermediate/unwanted.identity.tsv
touch $TMP/intermediate/distraction.fasta
touch $TMP/intermediate/ignore-me.onebp.tsv
# Run again with some explicit options set (shouldn't change output)
rm -rf $TMP/output/*
thapbi_pict pipeline -i tests/reads/ -s $TMP/intermediate -o $TMP/output -m onebp -a 250 -r report
diff $TMP/intermediate/DNAMIX_S95_L001.fasta tests/prepare-reads/DNAMIX_S95_L001.fasta
diff $TMP/output/report.samples.txt tests/pipeline/thapbi-pict.samples.onebp.txt
diff $TMP/output/report.samples.tsv tests/pipeline/thapbi-pict.samples.onebp.tsv
diff $TMP/output/report.reads.tsv tests/pipeline/thapbi-pict.reads.onebp.tsv
diff $TMP/output/report.edit-graph.xgmml tests/pipeline/thapbi-pict.edit-graph.xgmml

echo "$0 - test_pipeline.sh passed"
