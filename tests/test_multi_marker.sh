#!/bin/bash

# Copyright 2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -euo pipefail

export TMP=${TMP:-/tmp/thapbi_pict}/multi_marker
rm -rf $TMP
mkdir -p $TMP

export DB=$TMP/pooled.sqlite
mkdir $TMP/intermediate

echo "====================================="
echo "Setting up empty DB for multi-markers"
echo "====================================="

# This is a cut down version of examples/endangered_species/run.sh
function import_marker {
    # Takes arguments via variable names
    echo "#${NAME}" > $TMP/empty.fasta
    thapbi_pict import -d $DB -i $TMP/empty.fasta -x -s ";" \
                -k $NAME --left $LEFT --right $RIGHT -x --minlen $MINLEN
}

NAME=16S
LEFT=CGCCTGTTTATCAAAAACAT
RIGHT=CCGGTCTGAACTCAGATCACGT
MINLEN=200
import_marker  # call function above

NAME=Mini-16S
LEFT=AYAAGACGAGAAGACCC
RIGHT=GATTGCGCTGTTATTCC
MINLEN=200
import_marker  # call function above

NAME=Mini-COI
LEFT=GGWACWGGWTGAACWGTWTAYCCYCC
RIGHT=TAIACYTCIGGRTGICCRAARAAYCA
MINLEN=200
import_marker  # call function above

NAME=Mini-cyt-b
LEFT=CCATCCAACATCTCAGCATGATGAAA
RIGHT=CCCTCAGAATGATATTTGTCCTCA
MINLEN=200
import_marker  # call function above

NAME=rbcL
LEFT=ATGTCACCACAAACAGAGACTAAAGC
RIGHT=GTAAAATCAAGTCCACCRCG
MINLEN=100  # Even shorter than author's 200 default
import_marker  # call function above

NAME=Mini-rbcL
LEFT=GTTGGATTCAAAGCTGGTGTTA
RIGHT=CVGTCCAMACAGTWGTCCATGT
MINLEN=140  # Shorter!
import_marker  # call function above

NAME=trnL-UAA
LEFT=CGAAATCGGTAGACGCTACG
RIGHT=GGGGATAGAGGGACTTGAAC
MINLEN=200
import_marker  # call function above

NAME=trnL-P6-loop
LEFT=GGGCAATCCTGAGCCAA
RIGHT=CCATTGAGTCTCTGCACCTATC
MINLEN=10  # Much shorter!
import_marker  # call function above

NAME=ITS2
LEFT=ATGCGATACTTGGTGTGAAT
RIGHT=GACGCTTCTCCAGACTACAAT
MINLEN=100  # Shorter!
import_marker  # call function above

echo "====================="
echo "Running prepare-reads"
echo "====================="

thapbi_pict prepare-reads \
    -i tests/multi_marker/raw_data/ \
    -o $TMP/intermediate \
    -d $DB -a 10

for f in tests/multi_marker/intermediate/*/EM_1_sample.fasta; do
    export m=${f%/*}
    export m=${m##*/}
    echo "Checking output for $m marker..."
    echo diff $TMP/intermediate/$m/EM_1_sample.fasta $f
    diff $TMP/intermediate/$m/EM_1_sample.fasta $f
done

echo "$0 - test_multi_marker.sh passed"
