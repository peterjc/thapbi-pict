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

if [ ! -f tests/run_tests.sh ]; then
    echo "Please run tests/run_tests.sh from the top level directory"
    false
fi

time tests/test_build_db.sh

time tests/test_woody_hosts.sh

time tests/test_ena-submit.sh
time tests/test_dump.sh
time tests/test_load-tax.sh
time tests/test_curated-import.sh
time tests/test_conflicts.sh

if ! [ -x "$(command -v cutadapt)" ]; then
    echo 'WARNING: cutadapt not installed, skipping some tests'
else
    time tests/test_ncbi-import.sh
fi

# Currently can't easily install these on TravisCI
if ! [ -x "$(command -v flash)" ]; then
    echo 'WARNING: flash not installed, skipping some tests'
elif ! [ -x "$(command -v cutadapt)" ]; then
    echo 'WARNING: cutadapt not installed, skipping some tests'
else
    time tests/test_prepare-reads.sh
    time tests/test_pipeline.sh
fi

time tests/test_fasta-nr.sh
time tests/test_classify.sh
time tests/test_curated-seq.sh
time tests/test_assess.sh
time tests/test_summary.sh
time tests/test_edit-graph.sh

echo "================="
echo "Test suite passed"
echo "================="
