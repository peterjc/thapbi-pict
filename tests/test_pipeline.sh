#!/bin/bash
IFS=$'\n\t'
set -eux
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking pipeline"
thapbi_pict pipeline 2>&1 | grep "the following arguments are required"
set -o pipefail

# TODO...

echo "$0 - test_pipeline.sh passed"
