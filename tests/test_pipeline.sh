#!/bin/bash

# Copyright 2019-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp/thapbi_pict}/pipeline
rm -rf $TMP
mkdir -p $TMP

echo "================="
echo "Checking pipeline"
echo "================="
set -x
thapbi_pict pipeline 2>&1 | grep "the following arguments are required"
set -o pipefail

rm -rf $TMP/intermediate $TMP/output
mkdir $TMP/intermediate $TMP/output
thapbi_pict pipeline -s $TMP/intermediate -o $TMP/output/thapbi-pict -i tests/reads/
diff $TMP/intermediate/ITS1/DNAMIX_S95_L001.fasta tests/prepare-reads/DNAMIX_S95_L001.fasta
diff $TMP/output/thapbi-pict.ITS1.all_reads.fasta tests/pipeline/thapbi-pict.all_reads.fasta
diff $TMP/output/thapbi-pict.ITS1.samples.onebp.txt tests/pipeline/thapbi-pict.samples.onebp.txt
diff $TMP/output/thapbi-pict.ITS1.samples.onebp.tsv tests/pipeline/thapbi-pict.samples.onebp.tsv
diff $TMP/output/thapbi-pict.ITS1.reads.onebp.tsv tests/pipeline/thapbi-pict.reads.onebp.tsv

# Leaving the intermediate files in place... plus some stray files:
touch $TMP/intermediate/ITS1/unwanted.fasta
touch $TMP/intermediate/ITS1/unwanted.onebp.tsv
touch $TMP/intermediate/ITS1/unwanted.identity.tsv
touch $TMP/intermediate/ITS1/distraction.fasta
touch $TMP/intermediate/ITS1/ignore-me.onebp.tsv
# Run again with some explicit options set (shouldn't change output)
# Using --flip will have no effect as already have the intermediate files
rm -rf $TMP/output/*
thapbi_pict pipeline -i tests/reads/ -s $TMP/intermediate -o $TMP/output/report -m onebp -a 250 -d - --flip --metaencoding UTF-8
diff $TMP/intermediate/ITS1/DNAMIX_S95_L001.fasta tests/prepare-reads/DNAMIX_S95_L001.fasta
diff $TMP/output/report.ITS1.samples.onebp.txt tests/pipeline/thapbi-pict.samples.onebp.txt
diff $TMP/output/report.ITS1.samples.onebp.tsv tests/pipeline/thapbi-pict.samples.onebp.tsv
diff $TMP/output/report.ITS1.reads.onebp.tsv tests/pipeline/thapbi-pict.reads.onebp.tsv

# Clear the intermediate, run again with --merged-cache
rm -rf $TMP/intermediate_with_cache $TMP/output $TMP/merged_cache
mkdir $TMP/intermediate_with_cache $TMP/output $TMP/merged_cache
thapbi_pict pipeline --merged-cache $TMP/merged_cache -s $TMP/intermediate_with_cache -o $TMP/output/thapbi-pict -i tests/reads/
for F in $TMP/intermediate_with_cache/ITS1/*.fasta; do
    diff $F $TMP/intermediate/ITS1/${F##*/}
done
diff $TMP/intermediate/ITS1/DNAMIX_S95_L001.fasta tests/prepare-reads/DNAMIX_S95_L001.fasta
diff $TMP/output/thapbi-pict.ITS1.samples.onebp.txt tests/pipeline/thapbi-pict.samples.onebp.txt
diff $TMP/output/thapbi-pict.ITS1.samples.onebp.tsv tests/pipeline/thapbi-pict.samples.onebp.tsv
diff $TMP/output/thapbi-pict.ITS1.reads.onebp.tsv tests/pipeline/thapbi-pict.reads.onebp.tsv

echo "$0 - test_pipeline.sh passed"
