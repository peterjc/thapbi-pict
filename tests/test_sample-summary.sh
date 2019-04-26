#!/bin/bash
IFS=$'\n\t'
set -eux
# Note not using "set -o pipefile" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking sample-summary"
thapbi_pict sample-summary 2>&1 | grep "the following arguments are required"
thapbi_pict sample-summary -o '' tests/classify 2>&1 | grep "No output file specified"
set -o pipefail

# Passing filename, default method, explicit min abundance
rm -rf $TMP/human.txt $TMP/computer.tsv
thapbi_pict sample-summary -a 1 -r $TMP/human.txt -o $TMP/computer.tsv tests/classify/*.identity.tsv
diff $TMP/human.txt tests/sample-summary/classify.identity.txt
diff $TMP/computer.tsv tests/sample-summary/classify.identity.tsv

# Passing a folder, trying different methods
for M in identity onebp blast swarmid swarm; do
    rm -rf $TMP/human.txt $TMP/computer.tsv
    thapbi_pict sample-summary -m $M -r $TMP/human.txt -o $TMP/computer.tsv tests/classify/
    diff $TMP/human.txt tests/sample-summary/classify.$M.txt
    diff $TMP/computer.tsv tests/sample-summary/classify.$M.tsv
    # And again, but with metadata
    rm -rf $TMP/human.txt $TMP/computer.tsv
    thapbi_pict sample-summary -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5 -m $M -r $TMP/human.txt -o $TMP/computer.tsv tests/classify/
    diff $TMP/human.txt tests/sample-summary/classify-meta.$M.txt
    # This currently does not include metadata:
    diff $TMP/computer.tsv tests/sample-summary/classify.$M.tsv
done

echo "$0 - test_sample-summary.sh passed"
