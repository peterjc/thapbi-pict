#!/bin/bash

# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp/thapbi_pict}/load_tax
rm -rf $TMP
mkdir -p $TMP

echo "================="
echo "Checking load-tax"
echo "================="
set -x
thapbi_pict load-tax 2>&1 | grep "the following arguments are required"
set -o pipefail

# Same taxonomy as test_ncbi-import.sh (new style, ~100MB)
if [ ! -f "new_taxdump_2019-09-01.zip" ]; then curl -L -O "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_2019-09-01.zip"; fi
if [ ! -d "new_taxdump_2019-09-01" ]; then unzip new_taxdump_2019-09-01.zip nodes.dmp names.dmp -d new_taxdump_2019-09-01; fi

thapbi_pict load-tax -d "sqlite:///:memory:" -t new_taxdump_2019-09-01

# Same taxonomy as database/build_ITS1_DB.sh via tests/test_build_db.sh (old style, ~50MB)
export `grep ^TAX= database/build_ITS1_DB.sh`
if [ ! -f "${TAX}.zip" ]; then curl -L -O "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/${TAX}.zip"; fi
if [ ! -d "database/${TAX}" ]; then unzip ${TAX}.zip nodes.dmp names.dmp -d database/${TAX}; fi

# 2025987 = Nothophytophthora
export DB=$TMP/Nothophytophthora.sqlite
thapbi_pict load-tax -d $DB -t "database/${TAX}" -a 2025987
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "0" ]; then echo "Wrong data_source count"; false; fi
# 1 genus plus 6 accepted species:
if [ `sqlite3 $DB "SELECT COUNT(*) FROM taxonomy;"` -ne 7 ]; then echo "Wrong taxonomy count"; false; fi
# 4 unclassified species unused, plus "unclassified Nothophytophthora" as synonyms of the genus:
if [ `sqlite3 $DB "SELECT COUNT(*) FROM synonym;"` -ne 1 ]; then echo "Wrong synonym count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(synonym.id) FROM synonym JOIN taxonomy ON synonym.taxonomy_id==taxonomy.id WHERE taxonomy.genus='Nothophytophthora' AND taxonomy.species='';"` -ne 1 ]; then echo "Wrong synonym target"; false; fi

# Defaults to all Oomycetes
thapbi_pict load-tax -d "sqlite:///:memory:" -t "database/${TAX}"

# Request Phytophthora and Peronospora
thapbi_pict load-tax -d sqlite:///:memory: -t "database/${TAX}" -a 4783,70742

# Request all Peronosporales
thapbi_pict load-tax -d sqlite:///:memory: -t "database/${TAX}" -a 4776

# 6231 = Nematoda
thapbi_pict load-tax -d sqlite:///:memory: -t "database/${TAX}" -a 6231

# Check this error condition
set +o pipefail
thapbi_pict load-tax -d sqlite:///:memory: -t "database/${TAX}" -a 12908 2>&1 | grep "Could not identify any genus names"
set -o pipefail

echo "$0 - test_load-tax.sh passed"
