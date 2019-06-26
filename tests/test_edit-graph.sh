#!/bin/bash

# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eux
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking edit-graph"
thapbi_pict edit-graph 2>&1 | grep "If not using -i / --input argument, require -s / --showdb"
thapbi_pict edit-graph -d '' 2>&1 | grep "Require -d / --database and/or -i / --input argument"
set -o pipefail

# No database, small FASTA file, have to use explicit total abundance threshold
if [ `thapbi_pict edit-graph -d '' -i tests/prepare-reads/DNAMIX_S95_L001.fasta -t 200 | grep -c "<edge "` -ne 1 ]; then echo echo "Wrong edge count"; false; fi
# Loaded 7 unique sequences from 1 FASTA files.
# Minimum total abundance threshold 200 left 7 sequences from FASTA files.
# Computed 42 Levenshtein edit distances between 7 sequences.
# Will draw 2 nodes with at least one edge (5 are isolated sequences).
# Including high abundance isolated sequences, will draw 7 nodes.
# 1

# Same example as above with default xgmml output, but here different output formats:
if [ `thapbi_pict edit-graph -d '' -i tests/prepare-reads/DNAMIX_S95_L001.fasta -t 200 -f graphml | grep -c "<edge "` -ne 1 ]; then echo echo "Wrong edge count"; false; fi
if [ `thapbi_pict edit-graph -d '' -i tests/prepare-reads/DNAMIX_S95_L001.fasta -t 200 -f gexf | grep -c "<edge "` -ne 1 ]; then echo echo "Wrong edge count"; false; fi
# gml fails, https://github.com/networkx/networkx/issues/3471

# Same example, but PDF output (more dependencies):
thapbi_pict edit-graph -d '' -i tests/prepare-reads/DNAMIX_S95_L001.fasta -t 200 -f pdf | grep "%PDF-1.4"


# No database, generic FASTA file, have to use explicit abundance thresholds of 1:
if [ `thapbi_pict edit-graph -d '' -i database/legacy/Phytophthora_ITS_database_v0.005.fasta -a 1 -t 1 | grep -c "<node "` -ne 176 ]; then echo "Wrong node count"; false; fi
# WARNING: Sequence(s) in database/legacy/Phytophthora_ITS_database_v0.005.fasta not using MD5_abundance naming
# Loaded 176 unique sequences from 1 FASTA files.
# Minimum total abundance threshold 1 left 176 sequences from FASTA files.
# Computed 30800 Levenshtein edit distances between 176 sequences.
# Will draw 103 nodes with at least one edge (73 are isolated sequences).
# Including high abundance isolated sequences, will draw 176 nodes.
# 176

echo "$0 - test_edit-graph.sh passed"
