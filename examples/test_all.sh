#!/bin/bash
set -euo pipefail

# Assumes setup_all.sh has been run to download the data etc,
# and probably need to have used run_all.sh too.

# Only look at documented examples, and run them in that order:
for example in `grep "/index" ../docs/examples/index.rst | cut -f 1 -d "/" | cut -c 4-`; do
    echo "================================="
    echo "Checking docs for $example"
    echo "================================="
    cd $example
    ../../scripts/rst_doc_test.py ../../docs/examples/$example/*.rst
    cd ..
done

echo "================================="
echo "Checked all example documentation"
echo "================================="
