#!/bin/bash
IFS=$'\n\t'
set -eux

# Note not using "set -o pipefile" as want to use that
# with grep to check error messages

export TMP=${TMP:-/tmp}

echo "Checking prepare-reads"
thapbi_pict prepare-reads 2>&1 | grep "the following arguments are required"


echo "$0 passed"
