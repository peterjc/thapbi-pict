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
diff <(head -n 30 $TMP/intermediate/ITS1/DNAMIX_S95_L001.fasta) tests/prepare-reads/DNAMIX_S95_L001-a2-head.fasta
diff $TMP/output/thapbi-pict.ITS1.tally.tsv tests/pipeline/thapbi-pict.tally.tsv
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
thapbi_pict pipeline -i tests/reads/ -s $TMP/intermediate -o $TMP/output/report -m onebp -f 0 -a 100 -d - --flip --metaencoding UTF-8
diff <(head -n 30 $TMP/intermediate/ITS1/DNAMIX_S95_L001.fasta) tests/prepare-reads/DNAMIX_S95_L001-a2-head.fasta
diff $TMP/output/report.ITS1.samples.onebp.tsv tests/pipeline/thapbi-pict.samples.onebp.tsv
diff $TMP/output/report.ITS1.reads.onebp.tsv tests/pipeline/thapbi-pict.reads.onebp.tsv

# Clear the intermediate, run again with --merged-cache
rm -rf $TMP/intermediate_with_cache $TMP/output $TMP/merged_cache
mkdir $TMP/intermediate_with_cache $TMP/output $TMP/merged_cache
if [ -x "$(command -v biom)" ]; then
    thapbi_pict pipeline --merged-cache $TMP/merged_cache -s $TMP/intermediate_with_cache -o $TMP/output/thapbi-pict -i tests/reads/ --biom
    biom validate-table -i $TMP/output/thapbi-pict.ITS1.onebp.biom
else
    # can't use --biom option:
    thapbi_pict pipeline --merged-cache $TMP/merged_cache -s $TMP/intermediate_with_cache -o $TMP/output/thapbi-pict -i tests/reads/
fi
for F in $TMP/intermediate_with_cache/ITS1/*.fasta; do
    diff $F $TMP/intermediate/ITS1/${F##*/}
done
diff <(head -n 30 $TMP/intermediate/ITS1/DNAMIX_S95_L001.fasta) tests/prepare-reads/DNAMIX_S95_L001-a2-head.fasta
diff $TMP/output/thapbi-pict.ITS1.samples.onebp.tsv tests/pipeline/thapbi-pict.samples.onebp.tsv
diff $TMP/output/thapbi-pict.ITS1.reads.onebp.tsv tests/pipeline/thapbi-pict.reads.onebp.tsv

# Now with denoising, changes the counts but not the 10 sequences themselves
thapbi_pict pipeline --merged-cache $TMP/merged_cache -s $TMP/intermediate_with_cache \
            -o $TMP/output/unoise-l -i tests/reads/ --denoise unoise-l
diff <(cut -f 1-3 $TMP/output/unoise-l.ITS1.samples.onebp.tsv) \
     <(cut -f 1-3 tests/pipeline/thapbi-pict.samples.onebp.tsv)
if ! [ -x "$(command -v vsearch)" ]; then
    echo "Skipping testing using VSEARCH"
else
    thapbi_pict pipeline --merged-cache $TMP/merged_cache -s $TMP/intermediate_with_cache \
                -o $TMP/output/vsearch -i tests/reads/ --denoise vsearch
    diff <(cut -f 1-3 $TMP/output/vsearch.ITS1.samples.onebp.tsv) \
         <(cut -f 1-3 tests/pipeline/thapbi-pict.samples.onebp.tsv)
fi
if ! [ -x "$(command -v usearch)" ]; then
    echo "Skipping testing using USEARCH"
else
    thapbi_pict pipeline --merged-cache $TMP/merged_cache -s $TMP/intermediate_with_cache \
                -o $TMP/output/usearch -i tests/reads/ --denoise usearch
    diff <(cut -f 1-3 $TMP/output/usearch.ITS1.samples.onebp.tsv) \
         <(cut -f 1-3 tests/pipeline/thapbi-pict.samples.onebp.tsv)
fi

echo "$0 - test_pipeline.sh passed"
