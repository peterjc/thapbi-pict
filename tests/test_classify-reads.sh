#!/bin/bash
IFS=$'\n\t'
set -eux

# Note not using "set -o pipefile" as want to use that
# with grep to check error messages

export TMP=${TMP:-/tmp}

echo "Checking classify-reads"
thapbi_pict classify-reads 2>&1 | grep "the following arguments are required"

# Passing one filename:
thapbi_pict classify-reads -m identity -d $TMP/legacy_004_and_005_validated.sqlite database/legacy/database.fasta

# Passing one directory name (should get all three FASTA files):
thapbi_pict classify-reads -m identity -d $TMP/legacy_004_and_005_validated.sqlite database/legacy/

echo "$0 passed"
