#!/bin/bash
set -eup pipeline

# Intended to remove intermediate files, custom DB, but not summary reports.

for example in woody_hosts recycled_water fungal_mock microalgal_mock great_lakes fecal_sequel endangered_species; do
    echo "========================="
    echo "Cleaning $example"
    echo "========================="
    rm -rf $example/*.sqlite $example/*/*.sqlite

    # Do not touch woody_hosts/intermediate.tar.bz2
    rm -rf $example/intermediate $example/intermediate_* $example/tmp_merged

    # The XGMML files are not overwritten by default
    rm -rf $example/summary/*.xgmml
    echo "Done"
done

echo "=============================="
echo "Removed all intermediate files"
echo "=============================="
