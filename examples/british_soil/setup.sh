#!/bin/bash
set -eup pipeline

mkdir -p intermediate/ summary/

# ITS1 = PRJNA690943, COI-5P = PRJNA691575
for ACC in `grep ITS PRJNA690943.tsv | cut -f 1`; do
    # echo "Downloading $ACC"
    # Column 5 should have two URLs (R1 and R2), semi-colon separated:
    for URL in `grep ^$ACC PRJNA690943.tsv | cut -f 5 | sed "s/;/ /g"`; do
        FILE=raw_data/${URL##*/}
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
