#!/bin/bash
set -eup pipeline

mkdir -p intermediate/ summary/

if [ ! -f raw_data/MD5SUM.txt ]; then
    echo "ERROR: Missing raw_data/MD5SUM.txt"
    false
fi
for ACC in `grep ^SRR PRJNA379165.txt| cut -f 1`; do
    # echo "Downloading $ACC"
    # Column 7 should have two URLs (R1 and R2), semi-colon separated:
    for URL in `grep ^$ACC PRJNA379165.txt | cut -f 7 | sed "s/;/ /g"`; do
        NAME=${URL##*/}
        FILE=raw_data/$NAME
        # Avoiding leaving partial FASTQ if wget is interupted
        rm -rf $FILE.tmp
        if [ -f $FILE ]; then
            echo "Already have $FILE"
        else
            echo "Downloading $FILE"
            wget -O "${FILE}.tmp" "$URL"
            mv "${FILE}.tmp" "${FILE}"
        fi
    done
done
