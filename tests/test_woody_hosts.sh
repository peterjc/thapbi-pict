#!/bin/bash

# Copyright 2019-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -euo pipefail

export TMP=${TMP:-/tmp/thapbi_pict}/woody_hosts
rm -rf $TMP
mkdir -p $TMP

echo "Preparing sample data for woody hosts example"

mkdir $TMP/intermediate
mkdir $TMP/summary
mkdir $TMP/positive_controls/
for f in examples/woody_hosts/expected/*.known.tsv; do ln -s $PWD/$f $TMP/positive_controls/ ; done

# Idea here is to mimic what "thapbi_pict pipeline" would do if we had
# the FASTQ files here:
# thapbi_pict pipeline -i sample_data/raw_data/ \
#     -n tests/woody_hosts/raw_data/NEGATIVE*.fastq.gz \
#     -s $TMP/intermediate \
#     -o $TMP/summary/woody-hosts \
#     -t examples/woody_hosts/metadata.tsv \
#     -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -x 16


echo "=================================="
echo "Decompressing prepare-reads output"
echo "=================================="
time tar -jxvf tests/woody_hosts/intermediate.tar.bz2 -C $TMP/ | wc -l

if [ -d tests/woody_hosts/raw_data/ ]; then
    echo "================================="
    echo "Running woody hosts prepare-reads"
    echo "================================="
    mkdir $TMP/intermediate_new
    time thapbi_pict prepare-reads -i tests/woody_hosts/raw_data/ -o $TMP/intermediate_new -n tests/woody_hosts/raw_data/NEGATIVE*.fastq.gz
    for f in $TMP/intermediate/ITS1/*.fasta; do
        echo diff $f $TMP/intermediate_new/ITS1/${f##*/}
        diff $f $TMP/intermediate_new/ITS1/${f##*/}
    done
else
    echo "To check how much tests/woody_hosts/intermediate.tar.bz2 is out of date,"
    echo 'use: $ ln -s $PWD/examples/woody_hosts/raw_data tests/woody_hosts/'
    echo "and then re-run tests/test_woody_hosts.sh which will report first diff."
    echo "Replace the old FASTA files, then in the updated intermediate folder run:"
    echo '$ tar -cvjf intermediate.tar.bz2 intermediate/ITS1/*.fasta'
fi

echo "================================"
echo "Running woody hosts sample-tally"
echo "================================"
thapbi_pict sample-tally -i $TMP/intermediate/ITS1/*.fasta \
            -n $TMP/intermediate/ITS1/NEGATIVE_*.fasta \
            -o $TMP/woody_hosts.tally.tsv
echo diff $TMP/woody_hosts.tally.tsv tests/woody_hosts/all.tally.tsv
diff $TMP/woody_hosts.tally.tsv tests/woody_hosts/all.tally.tsv

echo "And with UNOISE..."
thapbi_pict sample-tally -i $TMP/intermediate/ITS1/*.fasta \
            -n $TMP/intermediate/ITS1/NEGATIVE_*.fasta \
            -o $TMP/woody_hosts.de-noise.tsv --denoise unoise-l
echo diff $TMP/woody_hosts.de-noise.tsv tests/woody_hosts/denoise.tally.tsv
diff $TMP/woody_hosts.de-noise.tsv tests/woody_hosts/denoise.tally.tsv

echo "============================"
echo "Running woody hosts fasta-nr"
echo "============================"
thapbi_pict fasta-nr -i $TMP/intermediate/ITS1/*.fasta -o $TMP/woody_hosts.all_reads.fasta
echo diff $TMP/woody_hosts.all_reads.fasta tests/woody_hosts/all.fasta
diff $TMP/woody_hosts.all_reads.fasta tests/woody_hosts/all.fasta

echo "============================"
echo "Running woody hosts classify"
echo "============================"
for M in onebp identity; do
    if [ -x "$(command -v biom)" ]; then
        thapbi_pict classify -i $TMP/woody_hosts.tally.tsv -m $M --biom
        biom validate-table -i $TMP/woody_hosts.$M.biom
    else
        thapbi_pict classify -i $TMP/woody_hosts.tally.tsv -m $M
    fi
    echo diff $TMP/woody_hosts.$M.tsv tests/woody_hosts/all.$M.tsv
    diff $TMP/woody_hosts.$M.tsv tests/woody_hosts/all.$M.tsv
done

echo "==========================="
echo "Running woody hosts summary"
echo "==========================="
time thapbi_pict summary -i $TMP/woody_hosts.onebp.tsv \
    -o $TMP/summary/no-metadata
ls $TMP/summary/no-metadata.*
if [ `grep -c -v "^#" $TMP/summary/no-metadata.reads.onebp.tsv` -ne 100 ]; then echo "Wrong unique sequence count"; false; fi
# Expect 99 + total line

time thapbi_pict summary -i $TMP/woody_hosts.onebp.tsv \
    -o $TMP/summary/with-metadata \
    -t examples/woody_hosts/metadata.tsv \
    -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -x 16
ls $TMP/summary/with-metadata.*

# Should be identical apart from row order if discard extra leading columns
# Discarding the header row as only one will still have hash at start
diff <(grep -v "^#" $TMP/summary/no-metadata.samples.onebp.tsv | sort) <(grep -v "^#" $TMP/summary/with-metadata.samples.onebp.tsv | cut -f 16- | sort)

if [ `grep -c -v "^#" $TMP/summary/with-metadata.reads.onebp.tsv` -ne 100 ]; then echo "Wrong unique sequence count"; false; fi
# Expect 99 + total line

echo "=============================="
echo "Running woody hosts edit-graph"
echo "=============================="
time thapbi_pict edit-graph -i $TMP/woody_hosts.onebp.tsv -o $TMP/summary/no-metadata.edit-graph.xgmml
if [ `grep -c "<node " $TMP/summary/no-metadata.edit-graph.xgmml` -ne 99 ]; then echo "Wrong node count"; false; fi
if [ `grep -c "<edge " $TMP/summary/no-metadata.edit-graph.xgmml` -ne 69 ]; then echo "Wrong edge count"; false; fi

echo "=========================="
echo "Running woody hosts assess"
echo "=========================="
time thapbi_pict assess -i $TMP/positive_controls/ $TMP/woody_hosts.onebp.tsv -o $TMP/DNA_MIXES.assess.tsv
echo diff $TMP/DNA_MIXES.assess.tsv tests/woody_hosts/DNA_MIXES.assess.tsv
diff $TMP/DNA_MIXES.assess.tsv tests/woody_hosts/DNA_MIXES.assess.tsv

echo "$0 - test_woody_hosts.sh passed"
