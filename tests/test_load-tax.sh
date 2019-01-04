#!/bin/bash
IFS=$'\n\t'
set -eux

# Note not using "set -o pipefile" as want to use that
# with grep to check error messages

export TMP=${TMP:-/tmp}

echo "Checking load-tax"
thapbi_pict load-tax 2>&1 | grep "the following arguments are required"

if [ ! -f "new_taxdump_2018-12-01.zip" ]; then curl -L -O "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_2018-12-01.zip"; fi
if [ ! -d "new_taxdump_2018-12-01" ]; then unzip new_taxdump_2018-12-01.zip -d new_taxdump_2018-12-01; fi

thapbi_pict load-tax -d "sqlite:///:memory:" -t new_taxdump_2018-12-01

if [ ! -f "taxdmp_2014-08-01.zip" ]; then curl -L -O "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2014-08-01.zip"; fi
if [ ! -d "taxdmp_2014-08-01" ]; then unzip taxdmp_2014-08-01.zip -d taxdmp_2014-08-01; fi

# Defaults to all Phytophthora
thapbi_pict load-tax -d "sqlite:///:memory:" -t taxdmp_2014-08-01

# Request Phytophthora and Peronospora
thapbi_pict load-tax -d sqlite:///:memory: -t taxdmp_2014-08-01 -a 4783,70742

# Request all Peronosporales
thapbi_pict load-tax -d sqlite:///:memory: -t taxdmp_2014-08-01 -a 4776

# Check this error condition
thapbi_pict load-tax -d sqlite:///:memory: -t taxdmp_2014-08-01 -a 12908 2>&1 | grep "Could not identify any genus names"

echo "$0 passed"
