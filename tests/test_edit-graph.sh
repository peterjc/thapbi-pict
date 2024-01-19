#!/bin/bash

# Copyright 2019-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

echo "==================="
echo "Checking edit-graph"
echo "==================="
set -x
thapbi_pict edit-graph 2>&1 | grep "If not using -i / --input argument, require -k / --marker"
thapbi_pict edit-graph -d '' 2>&1 | grep "Require -d / --database and/or -i / --input argument"
set -o pipefail

# No database, small FASTA file, have to use explicit total abundance threshold
diff --strip-trailing-cr tests/edit-graph/DNAMIX_S95_L001.xgmml <(thapbi_pict edit-graph -d '' -i tests/sample-tally/DNAMIX_S95_L001.tally.tsv -t 200 --min-samples 1)
# Loaded 7 unique sequences from 1 FASTA files.
# Minimum total abundance threshold 200 left 7 sequences from FASTA files.
# Computed 42 Levenshtein edit distances between 7 sequences.
# Will draw 2 nodes with at least one edge (5 are isolated sequences).
# Including high abundance isolated sequences, will draw 7 nodes.
# 1

# Same example as above with default xgmml output, but here different output formats:
diff --strip-trailing-cr tests/edit-graph/DNAMIX_S95_L001.tsv <(thapbi_pict edit-graph -d '' -i tests/sample-tally/DNAMIX_S95_L001.tally.tsv -t 200 -f matrix)
if [ `thapbi_pict edit-graph -d '' -i tests/sample-tally/DNAMIX_S95_L001.tally.tsv -t 200 -f graphml | grep -c "<edge "` -ne 1 ]; then echo echo "Wrong edge count"; false; fi
if [ `thapbi_pict edit-graph -d '' -i tests/sample-tally/DNAMIX_S95_L001.tally.tsv -t 200 -f gexf | grep -c "<edge "` -ne 1 ]; then echo echo "Wrong edge count"; false; fi
if [ `thapbi_pict edit-graph -d '' -i tests/sample-tally/DNAMIX_S95_L001.tally.tsv -t 200 -f gml | grep -c "  edge \["` -ne 1 ]; then echo echo "Wrong edge count"; false; fi

# Same example, but PDF output (more dependencies):
thapbi_pict edit-graph -d '' -i tests/sample-tally/DNAMIX_S95_L001.tally.tsv -t 200 -f pdf | grep "%PDF-1.4"

echo "$0 - test_edit-graph.sh passed"
