#!/bin/bash

# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "==================="
echo "Checking ena-submit"
echo "==================="
set -x

echo "Checking ena-submit"
thapbi_pict ena-submit 2>&1 | grep "the following arguments are required"
set -o pipefail

thapbi_pict ena-submit -i tests/reads/ -t tests/reads/metadata.tsv -c 5 -x 1 -o $TMP/ena_submit.tsv
diff $TMP/ena_submit.tsv tests/reads/ena_submit.tsv

thapbi_pict ena-submit -i tests/reads/ -t tests/reads/metadata.tsv -c 5 -x 1 -o $TMP/ena_submit_custom.tsv --library "Test set" --instrument "Illumina Widget" --design "Ad hoc" --protocol "Making it up" --insert 275
diff $TMP/ena_submit_custom.tsv tests/reads/ena_submit_custom.tsv

echo "$0 - test_ena-submit.sh passed"
