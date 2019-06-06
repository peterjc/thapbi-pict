#!/bin/bash
IFS=$'\n\t'
set -eux
# Note not using "set -o pipefail" until after check error message with grep

echo "Checking edit-graph"
thapbi_pict edit-graph -d - 2>&1 | grep "Require -d / --database and/or -i / --input argument"
set -o pipefail

# No database, small FASTA file, have to use explciit total abundance threshold
if [ `thapbi_pict edit-graph -d - -i tests/prepare-reads/DNAMIX_S95_L001.fasta -t 200 | grep -c "<edge"` -ne 1 ]; then echo echo "Wrong edge count"; false; fi
# Loaded 7 unique sequences from 1 FASTA files.
# Minimum total abundance threshold 200 left 7 sequences from FASTA files.
# Computed 42 Levenshtein edit distances between 7 sequences.
# Will draw 2 nodes with at least one edge (5 are isolated sequences).
# Including high abundance isolated sequences, will draw 7 nodes.
# 1


# No database, generic FASTA file, have to use explicit abundance thresholds of 1:
if [ `thapbi_pict edit-graph -d - -i database/legacy/Phytophthora_ITS_database_v0.005.fasta -a 1 -t 1 | grep -c "<node"` -ne 176 ]; then echo "Wrong node count"; false; fi
# WARNING: Sequence(s) in database/legacy/Phytophthora_ITS_database_v0.005.fasta not using MD5_abundance naming
# Loaded 176 unique sequences from 1 FASTA files.
# Minimum total abundance threshold 1 left 176 sequences from FASTA files.
# Computed 30800 Levenshtein edit distances between 176 sequences.
# Will draw 103 nodes with at least one edge (73 are isolated sequences).
# Including high abundance isolated sequences, will draw 176 nodes.
# 176

echo "$0 - test_edit-graph.sh passed"
