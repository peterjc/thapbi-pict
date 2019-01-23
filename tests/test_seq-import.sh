#!/bin/bash
IFS=$'\n\t'
set -eux

# Note not using "set -o pipefile" as want to use that
# with grep to check error messages

export TMP=${TMP:-/tmp}

echo "Checking seq-import"
thapbi_pict seq-import 2>&1 | grep "the following arguments are required"

# Cannot use validation without having some taxonomy entries
#thapbi_pict seq-import -d "sqlite:///:memory:" tests/seq-import/dup_seqs.fasta 2>&1 | grep "Taxonomy table empty"

# No tsv file given:
thapbi_pict seq-import -d "sqlite:///:memory:" tests/legacy-import/dup_seqs.fasta 2>&1 | grep "but missing"

echo "$0 passed"
