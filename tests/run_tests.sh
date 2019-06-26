#!/bin/bash

# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -euo pipefail

python -c "import thapbi_pict; print('Direct import says version ' + thapbi_pict.__version__)"
thapbi_pict -v
python -m thapbi_pict -v

tests/test_dump.sh
tests/test_load-tax.sh
tests/test_legacy-import.sh
tests/test_ncbi-import.sh

# Currently can't easily install these on TravisCI
if ! [ -x "$(command -v flash)" ]; then
    echo 'WARNING: flash not installed, skipping some tests'
elif ! [ -x "$(command -v cutadapt)" ]; then
    echo 'WARNING: cutadapt not installed, skipping some tests'
else
    tests/test_prepare-reads.sh
    tests/test_pipeline.sh
fi

tests/test_classify.sh
tests/test_seq-import.sh
tests/test_assess.sh
tests/test_sample-summary.sh
tests/test_read-summary.sh
tests/test_edit-graph.sh
