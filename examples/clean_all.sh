#!/bin/bash
set -euo pipefail

# Intended to remove intermediate files, custom DB, but not summary reports.

# Only look at documented examples, and run them in that order:
for example in `grep "/index" ../docs/examples/index.rst | cut -f 1 -d "/" | cut -c 4-`; do
    echo "========================="
    echo "Cleaning $example"
    echo "========================="
    rm -rf $example/*.sqlite $example/*/*.sqlite

    # Do not touch $example/tmp_merged/

    # Do not touch woody_hosts/intermediate.tar.bz2
    rm -rf $example/intermediate $example/intermediate_*

    # The classifier TSV is not (currently) overwritten by default:
    rm -rf $example/summary/*.all_reads.*.tsv

    # The XGMML files are not overwritten by default
    rm -rf $example/summary/*.xgmml
    echo "Done"
done

echo "=============================="
echo "Removed all intermediate files"
echo "=============================="
