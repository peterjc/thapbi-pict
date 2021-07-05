#!/bin/bash
set -euo pipefail

# Assumes setup_all.sh has been run to download the data etc,
# and probably need to have used run_all.sh too.

for example in woody_hosts recycled_water fungal_mock microalgal_mock great_lakes fecal_sequel soil_nematodes endangered_species; do
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
