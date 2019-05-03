#!/bin/bash
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
diff $TMP/output/sample-summary.txt tests/pipeline/sample-summary.onebp.txt
diff $TMP/output/sample-summary.tsv tests/pipeline/sample-summary.onebp.tsv
diff $TMP/output/plate-summary.tsv tests/pipeline/plate-summary.onebp.tsv

# Leaving the intermediate files in place... plus some stray files:
touch $TMP/intermediate/unwanted.fasta
touch $TMP/intermediate/unwanted.onebp.tsv
touch $TMP/intermediate/unwanted.identity.tsv
touch $TMP/intermediate/distraction.fasta
touch $TMP/intermediate/ignore-me.onebp.tsv
# Run again with some explicit options set (shouldn't change output)
rm -rf $TMP/output/*
thapbi_pict pipeline -i tests/reads/ -s $TMP/intermediate -o $TMP/output -m onebp -a 250
diff $TMP/intermediate/DNAMIX_S95_L001.fasta tests/prepare-reads/DNAMIX_S95_L001.fasta
diff $TMP/output/sample-summary.txt tests/pipeline/sample-summary.onebp.txt
diff $TMP/output/sample-summary.tsv tests/pipeline/sample-summary.onebp.tsv
diff $TMP/output/plate-summary.tsv tests/pipeline/plate-summary.onebp.tsv

echo "$0 - test_pipeline.sh passed"
