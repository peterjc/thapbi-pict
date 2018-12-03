#!/bin/bash
IFS=$'\n\t'
set -eux

# Note not using "set -o pipefile" as want to use that
# with grep to check error messages

export TMP=${TMP:-/tmp}

echo "Checking load-tax"
thapbi_pict load-tax 2>&1 | grep "the following arguments are required"

if [ ! -f "new_taxdump_2018-12-01.zip" ]; then curl -L -O "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_2018-12-01.zip"; fi
if [ ! -d "new_taxdump_2018-12-01" ]; then unzip new_taxdump_2018-12-01.zip -d new_taxdump_2018-12-01; fi

thapbi_pict load-tax -d "sqlite:///:memory:" -t new_taxdump_2018-12-01 -v
