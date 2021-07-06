#!/bin/bash
set -euo pipefail

# Intended to remove intermediate files, custom DB, but not summary reports.

for example in woody_hosts recycled_water fungal_mock great_lakes fecal_sequel soil_nematodes endangered_species; do
    echo "========================="
    echo "Cleaning $example"
    echo "========================="
    rm -rf $example/*.sqlite $example/*/*.sqlite

    # Do not touch $example/tmp_merged/

    # Do not touch woody_hosts/intermediate.tar.bz2
    rm -rf $example/intermediate $example/intermediate_*

    # The XGMML files are not overwritten by default
    rm -rf $example/summary/*.xgmml
    echo "Done"
done

echo "=============================="
echo "Removed all intermediate files"
echo "=============================="
