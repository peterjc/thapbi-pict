#!/bin/bash
IFS=$'\n\t'
set -eux
# Note not using "set -o pipefail" until after check error message with grep

# Note all tests here (initially) using default database:

echo "Checking edit-graph"
thapbi_pict edit-graph -d - 2>&1 | grep "Require -d / --database and/or -i / --input argument"
set -o pipefail

# No database, small FASTA file, have to use explciit total abundance threshold
if [ `thapbi_pict edit-graph -d - -i tests/prepare-reads/DNAMIX_S95_L001.fasta -t 200 | grep -c "<edge"` -ne 1 ]; then echo echo "Wrong edge count"; false; fi
# Loaded 7 unique sequences from 1 FASTA files.
# Minimum total abundance threshold 200 left 7 sequences from FASTA files.
# Computed 42 Levenshtein edit distances between 7 sequences.
# Max edit distance 3 gave 6 connected components (size 2 to 1).
# Dropped 0 sequences with no siblings within maximum edit distance 3.
# 1
#
# i.e. GraphXML is missing the single unconnected nodes, so just get
# two nodes with one edge in this example


# No database, generic FASTA file, have to use explicit abundance thresholds of 1:
if [ `thapbi_pict edit-graph -d - -i database/legacy/Phytophthora_ITS_database_v0.005.fasta -a 1 -t 1 | grep -c "<node"` -ne 103 ]; then echo "Wrong node count"; false; fi

echo "$0 - test_edit-graph.sh passed"
