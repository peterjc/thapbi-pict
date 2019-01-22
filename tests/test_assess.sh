#!/bin/bash
IFS=$'\n\t'
set -eux

# Note not using "set -o pipefile" as want to use that
# with grep to check error messages

export TMP=${TMP:-/tmp}

echo "Checking assess"

if [ ! -f $TMP/DNAMIX_S95_L001.fasta ]; then echo "Run test_prepare.sh to setup test input"; false; fi
if [ ! -f $TMP/DNAMIX_S95_L001.swarm.tsv ]; then echo "Run test_classify.sh to setup test input"; false; fi
if [ ! -f $TMP/DNAMIX_S95_L001.identity.tsv ]; then echo "Run test_classify.sh to setup test input"; false; fi

rm -rf $TMP/assess_swarm_vs_identity.tsv
rm -rf $TMP/confusion_swarm_vs_identity.tsv

# Don't have a gold standard known truth to test this against, so test swarm vs identity
thapbi_pict assess $TMP/DNAMIX_S95_L001.swarm.tsv -m swarm -k identity -o $TMP/assess_swarm_vs_identity.tsv -c $TMP/confusion_swarm_vs_identity.tsv

# Check assessment output to stdout works (default):
thapbi_pict assess $TMP/DNAMIX_S95_L001.swarm.tsv -m swarm -k identity > $TMP/stdout.txt
diff $TMP/stdout.txt $TMP/assess_swarm_vs_identity.tsv

# Check confusion matrix output to stdout works:
thapbi_pict assess $TMP/DNAMIX_S95_L001.swarm.tsv -m swarm -k identity -o $TMP/assess_swarm_vs_identity.tsv -c - > $TMP/stdout.txt
diff $TMP/stdout.txt $TMP/confusion_swarm_vs_identity.tsv

echo "$0 passed"
