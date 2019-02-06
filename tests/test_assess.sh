#!/bin/bash
IFS=$'\n\t'
set -eux
# Note not using "set -o pipefile" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking assess"
thapbi_pict assess 2>&1 | grep "the following arguments are required"
set -o pipefail

# Simple examples with expected output to compare against:
diff tests/assess/ex1.assess.tsv <(thapbi_pict assess tests/assess/ex1.known.tsv tests/assess/ex1.identity.tsv)
diff tests/assess/ex2.assess.tsv <(thapbi_pict assess tests/assess/ex2.known.tsv tests/assess/ex2.identity.tsv)
diff tests/assess/ex3.assess.tsv <(thapbi_pict assess tests/assess/ex3.known.tsv tests/assess/ex3.identity.tsv)
diff tests/assess/ex4.assess.tsv <(thapbi_pict assess tests/assess/ex4.known.tsv tests/assess/ex4.identity.tsv)
diff tests/assess/unclassified.assess.tsv <(thapbi_pict assess tests/assess/unclassified.known.tsv tests/assess/unclassified.identity.tsv)

if [ ! -f $TMP/thapbi_swarm/DNAMIX_S95_L001.swarm.tsv ]; then echo "Run test_classify.sh to setup test input"; false; fi
if [ ! -f $TMP/DNAMIX_S95_L001.identity.tsv ]; then echo "Run test_classify.sh to setup test input"; false; fi

rm -rf $TMP/assess_swarm_vs_identity.tsv
rm -rf $TMP/confusion_swarm_vs_identity.tsv

# Don't have a gold standard known truth to test this against, so test swarm vs identity
thapbi_pict assess $TMP/thapbi_swarm/DNAMIX_S95_L001.swarm.tsv $TMP/DNAMIX_S95_L001.identity.tsv -m swarm -k identity -o $TMP/assess_swarm_vs_identity.tsv -c $TMP/confusion_swarm_vs_identity.tsv

# Check assessment output to stdout works (default):
thapbi_pict assess $TMP/thapbi_swarm/DNAMIX_S95_L001.swarm.tsv $TMP/DNAMIX_S95_L001.identity.tsv -m swarm -k identity > $TMP/stdout.txt
diff $TMP/stdout.txt $TMP/assess_swarm_vs_identity.tsv

# Check confusion matrix output to stdout works, also give one input as a directory name
thapbi_pict assess $TMP/thapbi_swarm $TMP/DNAMIX_S95_L001.identity.tsv -m swarm -k identity -o $TMP/assess_swarm_vs_identity.tsv -c - > $TMP/stdout.txt
diff $TMP/stdout.txt $TMP/confusion_swarm_vs_identity.tsv

echo "$0 - test_assess.sh passed"
