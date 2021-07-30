#!/bin/bash
set -euo pipefail

echo "NOTE: Expected first time run time is a few hours"
echo "(or about 20 minutes from the merged reads cache),"
echo "repeat runs take seconds just to regenerate reports."
echo

mkdir -p tmp_merged/ intermediate/ summary/

# ILL_LEFT=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
# ILL_RIGHT=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG

# Quoting table 2:
#
# Primer Sequence (from 5ʹ end) Source
# NF1  GGTGGTGCATGGCCGTTCTTAGTT Porazinska et al. 2009
# 18Sr2b TACAAAGGGCAGGGACGTAAT
# SSUF04  GCTTGTCTCAAAGATTAAGCC Blaxter et al. 1998
# SSUR22 GCCTGCTGCCTTCCTTGGA
# D3FA  GACCCGTCTTGAAACACGGA Nunn 1992
# D3BR CGGAAGGAACCAGCTACTA
# JB3 TTTTTTGGGCATCCTGAGGTTTAT Bowles et al. 1992
# JB5GED  AGCACCTAAACTTAAAACATARTGRAARTG Derycke et al. 2010

if [ ! -f references/pooled.sqlite ]; then
    echo =================
    echo Creating database
    echo =================

    NAME=NF1-18Sr2b
    LEFT=GGTGGTGCATGGCCGTTCTTAGTT
    RIGHT=TACAAAGGGCAGGGACGTAAT

    thapbi_pict import -d references/pooled.sqlite \
                -i references/$NAME.fasta -s ";" -x \
                -k $NAME --left $LEFT --right $RIGHT

    NAME=SSUF04-SSUR22
    LEFT=GCTTGTCTCAAAGATTAAGCC
    RIGHT=GCCTGCTGCCTTCCTTGGA

    thapbi_pict import -d references/pooled.sqlite \
                -i references/$NAME.fasta -s ";" -x \
                -k $NAME --left $LEFT --right $RIGHT

    NAME=D3Af-D3Br
    LEFT=GACCCGTCTTGAAACACGGA
    RIGHT=CGGAAGGAACCAGCTACTA

    thapbi_pict import -d references/pooled.sqlite \
                -i references/$NAME.fasta -s ";" -x \
                -k $NAME --left $LEFT --right $RIGHT

    NAME=JB3-JB5GED
    LEFT=TTTTTTGGGCATCCTGAGGTTTAT
    RIGHT=AGCACCTAAACTTAAAACATARTGRAARTG

    thapbi_pict import -d references/pooled.sqlite \
                -i references/$NAME.fasta -s ";" -x \
                -k $NAME --left $LEFT --right $RIGHT

fi

echo ================
echo Running analysis
echo ================

# Fraction 0.001 means 0.1%
thapbi_pict pipeline -d references/pooled.sqlite --synthetic '' \
            -i raw_data/ expected/ --merged-cache tmp_merged/ \
            -s intermediate/ -o summary/ -a 25 -f 0.001 \
            -t metadata.tsv -x 1 -c 4,3

echo ====
echo Done
echo ====
