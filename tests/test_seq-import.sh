#!/bin/bash
IFS=$'\n\t'
set -eux
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking seq-import"
thapbi_pict seq-import 2>&1 | grep "the following arguments are required"
# Cannot use validation without having some taxonomy entries
#thapbi_pict seq-import -d "sqlite:///:memory:" tests/seq-import/dup_seqs.fasta 2>&1 | grep "Taxonomy table empty"
# No tsv file given:
thapbi_pict seq-import -d "sqlite:///:memory:" tests/legacy-import/dup_seqs.fasta 2>&1 | grep "Need \*.fasta files with matching \*.known.tsv classification"
set -o pipefail

echo "$0 - test_seq-import.sh passed"
