#!/bin/bash

# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -euo pipefail

export TMP=${TMP:-/tmp}

echo "Preparing sample data for woody hosts example"

rm -rf $TMP/woody_hosts
mkdir $TMP/woody_hosts
mkdir $TMP/woody_hosts/intermediate
mkdir $TMP/woody_hosts/summary
mkdir $TMP/woody_hosts/positive_controls/
for f in examples/woody_hosts/expected/*.known.tsv; do ln -s $PWD/$f $TMP/woody_hosts/positive_controls/ ; done

# Idea here is to mimic what "thapbi_pict pipeline" would do if we had
# the FASTQ files here:
# thapbi_pict pipeline -i sample_data/raw_data/ \
#     -s $TMP/woody_hosts/intermediate \
#     -o $TMP/woody_hosts/summary -r woody-hosts \
#     -t examples/woody_hosts/metadata.tsv \
#     -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -x 16 -f 20


echo "=================================="
echo "Decompressing prepare-reads output"
echo "=================================="
time tar -jxvf examples/woody_hosts/intermediate.tar.bz2 -C $TMP/woody_hosts/ | wc -l

if [ -d tests/woody_hosts/raw_data/ ]; then
    echo "================================="
    echo "Running woody hosts prepare-reads"
    echo "================================="
    mkdir $TMP/woody_hosts/intermediate_new
    time thapbi_pict prepare-reads -i tests/woody_hosts/raw_data/ -o $TMP/woody_hosts/intermediate_new -n tests/woody_hosts/raw_data/NEGATIVE*.fastq.gz
    for f in $TMP/woody_hosts/intermediate/*.fasta; do
        echo diff $f $TMP/woody_hosts/intermediate_new/${f##*/}
        diff $f $TMP/woody_hosts/intermediate_new/${f##*/}
    done
fi

echo "======================================="
echo "Running woody hosts fasta-nr & classify"
echo "======================================="
thapbi_pict fasta-nr -i $TMP/woody_hosts/intermediate/*.fasta -o $TMP/woody_hosts/all.fasta
echo $TMP/woody_hosts/all.fasta tests/woody_hosts/all.fasta
diff $TMP/woody_hosts/all.fasta tests/woody_hosts/all.fasta
for M in onebp identity blast; do
    # Writing to stdout to set a single filename.
    # Discarding the comment column, and the header,
    # leaving the most stable core part of the output
    thapbi_pict classify -i $TMP/woody_hosts/all.fasta -o - -m $M | grep -v "^#" | cut -f 1-3 > $TMP/woody_hosts/all.$M.tsv
    echo diff $TMP/woody_hosts/all.$M.tsv tests/woody_hosts/all.$M.tsv
    diff $TMP/woody_hosts/all.$M.tsv tests/woody_hosts/all.$M.tsv
done

echo "============================"
echo "Running woody hosts classify"
echo "============================"
# Default for -o should be the same next to the inputs, which is fine
time thapbi_pict classify -i $TMP/woody_hosts/intermediate/

echo "==========================="
echo "Running woody hosts summary"
echo "==========================="
time thapbi_pict summary -i $TMP/woody_hosts/intermediate/ \
    -o $TMP/woody_hosts/summary/ -r no-metadata
ls $TMP/woody_hosts/summary/no-metadata.*
if [ `grep -c -v "^#" $TMP/woody_hosts/summary/no-metadata.reads.onebp.tsv` -ne 100 ]; then echo "Wrong unique sequence count"; false; fi
# Expect 99 + total line

time thapbi_pict summary -i $TMP/woody_hosts/intermediate/ \
    -o $TMP/woody_hosts/summary/ -r with-metadata \
    -t examples/woody_hosts/metadata.tsv \
    -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -x 16 -f 20
ls $TMP/woody_hosts/summary/with-metadata.*
if [ `grep -c "^Site: " "$TMP/woody_hosts/summary/with-metadata.samples.onebp.txt"` -ne 17 ]; then echo "Wrong site count"; false; fi
if [ `grep -c "^Sequencing sample: " "$TMP/woody_hosts/summary/with-metadata.samples.onebp.txt"` -ne 122 ]; then echo "Wrong sample count"; false; fi

# Should be identical apart from row order if discard extra leading columns
# Discarding the header row as only one will still have hash at start
diff <(grep -v "^#" $TMP/woody_hosts/summary/no-metadata.samples.onebp.tsv | sort) <(grep -v "^#" $TMP/woody_hosts/summary/with-metadata.samples.onebp.tsv | cut -f 16- | sort)

if [ `grep -c -v "^#" $TMP/woody_hosts/summary/with-metadata.reads.onebp.tsv` -ne 100 ]; then echo "Wrong unique sequence count"; false; fi
# Expect 99 + total line

echo "=============================="
echo "Running woody hosts edit-graph"
echo "=============================="
time thapbi_pict edit-graph -i $TMP/woody_hosts/intermediate/ -o $TMP/woody_hosts/summary/no-metadata.edit-graph.xgmml
if [ `grep -c "<node " $TMP/woody_hosts/summary/no-metadata.edit-graph.xgmml` -ne 99 ]; then echo "Wrong node count"; false; fi
if [ `grep -c "<edge " $TMP/woody_hosts/summary/no-metadata.edit-graph.xgmml` -ne 69 ]; then echo "Wrong edge count"; false; fi

echo "=========================="
echo "Running woody hosts assess"
echo "=========================="
time thapbi_pict assess -i $TMP/woody_hosts/positive_controls/ $TMP/woody_hosts/intermediate/ -o $TMP/woody_hosts/DNA_MIXES.assess.tsv
echo diff $TMP/woody_hosts/DNA_MIXES.assess.tsv tests/woody_hosts/DNA_MIXES.assess.tsv
diff $TMP/woody_hosts/DNA_MIXES.assess.tsv tests/woody_hosts/DNA_MIXES.assess.tsv

echo "$0 - test_woody_hosts.sh passed"
