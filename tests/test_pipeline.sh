#!/bin/bash

# Copyright 2019-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "================="
echo "Checking pipeline"
echo "================="
set -x
thapbi_pict pipeline 2>&1 | grep "the following arguments are required"
set -o pipefail

rm -rf $TMP/intermediate $TMP/output
mkdir $TMP/intermediate $TMP/output
thapbi_pict pipeline -s $TMP/intermediate -o $TMP/output -i tests/reads/
diff $TMP/intermediate/DNAMIX_S95_L001.fasta tests/prepare-reads/DNAMIX_S95_L001.fasta
diff $TMP/output/thapbi-pict.all_reads.fasta tests/pipeline/thapbi-pict.all_reads.fasta
diff $TMP/output/thapbi-pict.samples.onebp.txt tests/pipeline/thapbi-pict.samples.onebp.txt
diff $TMP/output/thapbi-pict.samples.onebp.tsv tests/pipeline/thapbi-pict.samples.onebp.tsv
diff $TMP/output/thapbi-pict.reads.onebp.tsv tests/pipeline/thapbi-pict.reads.onebp.tsv
diff $TMP/output/thapbi-pict.edit-graph.onebp.xgmml tests/pipeline/thapbi-pict.edit-graph.onebp.xgmml

# Leaving the intermediate files in place... plus some stray files:
touch $TMP/intermediate/unwanted.fasta
touch $TMP/intermediate/unwanted.onebp.tsv
touch $TMP/intermediate/unwanted.identity.tsv
touch $TMP/intermediate/distraction.fasta
touch $TMP/intermediate/ignore-me.onebp.tsv
# Run again with some explicit options set (shouldn't change output)
# Using --flip will have no effect as already have the intermediate files
rm -rf $TMP/output/*
thapbi_pict pipeline -i tests/reads/ -s $TMP/intermediate -o $TMP/output -m onebp -a 250 -r report --hmm thapbi_pict/hmm/controls.hmm --flip
diff $TMP/intermediate/DNAMIX_S95_L001.fasta tests/prepare-reads/DNAMIX_S95_L001.fasta
diff $TMP/output/report.samples.onebp.txt tests/pipeline/thapbi-pict.samples.onebp.txt
diff $TMP/output/report.samples.onebp.tsv tests/pipeline/thapbi-pict.samples.onebp.tsv
diff $TMP/output/report.reads.onebp.tsv tests/pipeline/thapbi-pict.reads.onebp.tsv
diff $TMP/output/report.edit-graph.onebp.xgmml tests/pipeline/thapbi-pict.edit-graph.onebp.xgmml


# Clear the intermediate, run again with --merged-cache
rm -rf $TMP/intermediate_with_cache $TMP/output $TMP/merged_cache
mkdir $TMP/intermediate_with_cache $TMP/output $TMP/merged_cache
thapbi_pict pipeline --merged-cache $TMP/merged_cache -s $TMP/intermediate_with_cache -o $TMP/output -i tests/reads/
for F in $TMP/intermediate_with_cache/*.fasta; do
    diff $F $TMP/intermediate/${F##*/}
done
diff $TMP/intermediate/DNAMIX_S95_L001.fasta tests/prepare-reads/DNAMIX_S95_L001.fasta
diff $TMP/output/thapbi-pict.samples.onebp.txt tests/pipeline/thapbi-pict.samples.onebp.txt
diff $TMP/output/thapbi-pict.samples.onebp.tsv tests/pipeline/thapbi-pict.samples.onebp.tsv
diff $TMP/output/thapbi-pict.reads.onebp.tsv tests/pipeline/thapbi-pict.reads.onebp.tsv
diff $TMP/output/thapbi-pict.edit-graph.onebp.xgmml tests/pipeline/thapbi-pict.edit-graph.onebp.xgmml

echo "$0 - test_pipeline.sh passed"
