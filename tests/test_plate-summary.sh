#!/bin/bash
IFS=$'\n\t'
set -eux
# Note not using "set -o pipefile" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking plate-summary"
thapbi_pict plate-summary 2>&1 | grep "the following arguments are required"
set -o pipefail

# One method:
thapbi_pict plate-summary $TMP/DNAMIX_S95_L001.fasta $TMP/DNAMIX_S95_L001.identity.tsv -o $TMP/plate-summary_identity.tsv

# Two methods:
thapbi_pict plate-summary $TMP/DNAMIX_S95_L001.fasta  $TMP/thapbi_swarm/ $TMP/thapbi_blast/ -m blast,swarm -o $TMP/plate-summary_swarm_and_blast.tsv

echo "$0 - test_plate-summary.sh passed"
