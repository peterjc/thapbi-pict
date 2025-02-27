#!/bin/bash

# Copyright 2019-2024 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

# Note all tests here (initially) using default database:

echo "============="
echo "Checking dump"
echo "============="
set -x
thapbi_pict dump -s fallax 2>&1 | grep "species requires a single genus"
set -o pipefail

if [ "$(thapbi_pict dump | grep -c -v "^#")" -ne 4688 ]; then
    echo "Wrong source count for table dump"
    false
fi
if [ "$(thapbi_pict dump -f fasta | grep -c "^>")" -ne 1552 ]; then
    echo "Wrong sequence count for fasta dump"
    false
fi

if [ "$(thapbi_pict dump --minimal | grep -c -v "^#")" -ne 1552 ]; then
    echo "Wrong sequence count for minimal table dump"
    false
fi
if [ "$(thapbi_pict dump -m -f fasta | grep -c "^>")" -ne 1552 ]; then
    echo "Wrong sequence count for minimal fasta dump"
    false
fi

# With genus filter,
if [ "$(thapbi_pict dump -f txt -g Phytophthora | grep -v -c "^#")" -ne 3028 ]; then
    echo "Wrong source count for Phytophthora genus"
    false
fi
if [ "$(thapbi_pict dump -f fasta -g Phytophthora -m | grep -c "^>")" -ne 830 ]; then
    echo "Wrong sequence count for Phytophthora genus"
    false
fi

# With genus and species filter,
if [ "$(thapbi_pict dump -f txt -g Phytophthora -s "fallax, andina" | grep -v -c "^#")" -ne 7 ]; then
    echo "Wrong source count for two species"
    false
fi

echo "$0 - test_dump.sh passed"
