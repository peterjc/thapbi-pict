#!/bin/bash

# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eux
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking ena-submit"
thapbi_pict ena-submit 2>&1 | grep "the following arguments are required"
set -o pipefail

# TODO
# thapbi_pict ena-submit -i tests/reads/ -t tests/reads/metadata.tsv -c 2,3,4 -x 1 -m tests/reads/ena_mapping.tsv

echo "$0 - test_ena-submit.sh passed"
