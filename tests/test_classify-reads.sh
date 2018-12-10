#!/bin/bash
IFS=$'\n\t'
set -eux

# Note not using "set -o pipefile" as want to use that
# with grep to check error messages

export TMP=${TMP:-/tmp}

echo "Checking classify-reads"
thapbi_pict classify-reads 2>&1 | grep "the following arguments are required"
thapbi_pict classify-reads -d "sqlite:///:memory:" hypothetical_example.fasta 2>&1 | grep "cannot classify anything"

export DB=$TMP/legacy_004_and_005_validated.sqlite
if [ ! -f $DB ]; then echo "Run test_legacy-import.sh to setup test DB"; false; fi

# Passing one filename:
thapbi_pict classify-reads -m identity -d $DB database/legacy/database.fasta

# Passing one directory name (should get all three FASTA files):
thapbi_pict classify-reads -m identity -d $DB database/legacy/

echo "$0 passed"
