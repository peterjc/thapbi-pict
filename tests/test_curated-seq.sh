#!/bin/bash

# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "===================="
echo "Checking curated-seq"
echo "===================="
set -x
thapbi_pict curated-seq 2>&1 | grep "the following arguments are required"

# No tsv file given:
thapbi_pict curated-seq -i tests/curated-import/dup_seqs.fasta 2>&1 | grep "Need \*.fasta files with matching \*.known.tsv classification"

# Not all sequences with a species:
thapbi_pict curated-seq -i tests/classify/P-infestans-T30-4.* -m identity 2>&1 | grep "No species for 4b639550662bfcf27f9face76af10d6b"

set -o pipefail

# Single FASTA and known.tsv pairing:
if [ `thapbi_pict curated-seq -i tests/classify/P-infestans-T30-4.* | grep "^>.* Phytophthora infestans" -c` -n 3 ]; then echo "Wrong FASTA record count"; false; fi

# Single FASTA and non-default method.tsv:
if [ `thapbi_pict curated-seq -i tests/classify/P-infestans-T30-4.* -m onebp | grep "^>.* Phytophthora andina;Phytophthora infestans;Phytophthora ipomoeae" -c` -n 3 ]; then echo "Wrong FASTA record count"; false; fi

# TODO - Test a pair with a wildcard in the TSV file

echo "$0 - test_curated-seq.sh passed"
