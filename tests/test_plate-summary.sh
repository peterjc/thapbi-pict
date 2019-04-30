#!/bin/bash
IFS=$'\n\t'
set -eux
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking plate-summary"
thapbi_pict plate-summary 2>&1 | grep "the following arguments are required"
thapbi_pict plate-summary -o '' tests/classify 2>&1 | grep "No output file specified"
set -o pipefail

# One method:
thapbi_pict plate-summary tests/prepare-reads/DNAMIX_S95_L001.fasta $TMP/DNAMIX_S95_L001.identity.tsv -o $TMP/plate-summary_identity.tsv

# Two methods:
thapbi_pict plate-summary tests/prepare-reads/DNAMIX_S95_L001.fasta  $TMP/thapbi_swarm/ $TMP/thapbi_blast/ -m blast,swarm -o $TMP/plate-summary_swarm_and_blast.tsv

# With metadata
thapbi_pict plate-summary tests/classify -m onebp -o $TMP/plate-summary_onebp.tsv -t tests/classify/P-infestans-T30-4.meta.tsv -x 1 -c 2,3,4,5
diff $TMP/plate-summary_onebp.tsv tests/classify/P-infestans-T30-4.summary.tsv

echo "$0 - test_plate-summary.sh passed"
