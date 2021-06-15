#!/bin/bash
set -eup pipeline

echo NOTE: Expected first time run time is about XX minutes,
echo repeat runs about 1 minute just to regenerate reports.
echo

mkdir -p tmp_merged/ intermediate/ summary/

ILL_LEFT=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
ILL_RIGHT=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG

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


# Takes arguments via variable names
function analyse {
    if [ ! -f references/$NAME.sqlite ]; then
        echo "Making DB $NAME.sqlite"
        thapbi_pict import -d references/$NAME.sqlite \
                    -i references/$NAME.fasta -s ";" -x
    fi

    echo "Running analysis"
    mkdir -p intermediate/$NAME/
    thapbi_pict pipeline -d references/$NAME.sqlite --merged-cache tmp_merged/ \
                 --left $ILL_LEFT$LEFT --right $ILL_RIGHT$RIGHT \
                -i raw_data/ expected/$NAME/ -s intermediate/$NAME/ \
                -o summary -r $NAME -t metadata.tsv -x 1 -c 4,3

    echo "$NAME done"
}

echo ========================
echo NF1/18Sr2b
echo ========================

NAME=NF1-18Sr2b
LEFT=GGTGGTGCATGGCCGTTCTTAGTT
RIGHT=TACAAAGGGCAGGGACGTAAT

analyse  # call function above

echo ========================
echo SSUF04/SSUR22
echo ========================

NAME=SSUF04-SSUR22
LEFT=GCTTGTCTCAAAGATTAAGCC
RIGHT=GCCTGCTGCCTTCCTTGGA

analyse  # call function above

echo ========================
echo D3Af/D3Br
echo ========================

NAME=D3Af-D3Br
LEFT=GACCCGTCTTGAAACACGGA
RIGHT=CGGAAGGAACCAGCTACTA

analyse  # call function above

echo ========================
echo JB3/JB5GED
echo ========================

NAME=JB3-JB5GED
LEFT=TTTTTTGGGCATCCTGAGGTTTAT
RIGHT=AGCACCTAAACTTAAAACATARTGRAARTG

analyse  # call function above

echo ====
echo Done
echo ====
