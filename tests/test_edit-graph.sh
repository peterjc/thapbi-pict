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

echo "==================="
echo "Checking edit-graph"
echo "==================="
set -x
thapbi_pict edit-graph 2>&1 | grep "If not using -i / --input argument, require -s / --showdb"
thapbi_pict edit-graph -d '' 2>&1 | grep "Require -d / --database and/or -i / --input argument"
set -o pipefail

# No database, small FASTA file, have to use explicit total abundance threshold
if [ `thapbi_pict edit-graph -d '' -i tests/prepare-reads/DNAMIX_S95_L001.fasta -t 200 -m - | grep -c "<edge "` -ne 1 ]; then echo echo "Wrong edge count"; false; fi
# Loaded 7 unique sequences from 1 FASTA files.
# Minimum total abundance threshold 200 left 7 sequences from FASTA files.
# Computed 42 Levenshtein edit distances between 7 sequences.
# Will draw 2 nodes with at least one edge (5 are isolated sequences).
# Including high abundance isolated sequences, will draw 7 nodes.
# 1

# Same example as above with default xgmml output, but here different output formats:
if [ `thapbi_pict edit-graph -d '' -i tests/prepare-reads/DNAMIX_S95_L001.fasta -t 200 -m - -f graphml | grep -c "<edge "` -ne 1 ]; then echo echo "Wrong edge count"; false; fi
if [ `thapbi_pict edit-graph -d '' -i tests/prepare-reads/DNAMIX_S95_L001.fasta -t 200 -m - -f gexf | grep -c "<edge "` -ne 1 ]; then echo echo "Wrong edge count"; false; fi
# gml fails, https://github.com/networkx/networkx/issues/3471

# Same example, but PDF output (more dependencies):
thapbi_pict edit-graph -d '' -i tests/prepare-reads/DNAMIX_S95_L001.fasta -t 200 -m - -f pdf | grep "%PDF-1.4"


# No database, generic FASTA file, have to use explicit abundance thresholds of 1:
if [ `thapbi_pict edit-graph -d '' -i database/Phytophthora_ITS1_curated.fasta -a 1 -t 1 -m - | grep -c "<node "` -ne 205 ]; then echo "Wrong node count"; false; fi
# WARNING: Sequence(s) in database/Phytophthora_ITS1_curated.fasta not using MD5_abundance naming
# Loaded 167 unique sequences from 1 FASTA files.
# Minimum total abundance threshold 1 left 167 sequences from FASTA files.
# Computed 13861 Levenshtein edit distances between 167 sequences.
# Will draw 104 nodes with at least one edge (63 are isolated sequences).
# Including high abundance isolated sequences, will draw 167 nodes.
# 167

echo "$0 - test_edit-graph.sh passed"
