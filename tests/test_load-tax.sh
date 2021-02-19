#!/bin/bash

# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

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
