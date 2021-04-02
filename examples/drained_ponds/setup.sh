#!/bin/bash
set -eup pipeline

mkdir -p intermediate/ summary/ expected/

for ACC in `grep ^SRR PRJNA638011.tsv | cut -f 1`; do
    # echo "Downloading $ACC"
    # Column 4 should have two URLs (R1 and R2), semi-colon separated:
    for URL in `grep ^$ACC PRJNA638011.tsv | cut -f 3 | sed "s/;/ /g"`; do
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

# There are four positive controls:
for ACC in `grep "positive" PRJNA638011.tsv | cut -f 1`; do
    # echo "Setting up expected classification"
    FILE=expected/$ACC.known.tsv
    if [ -f $FILE ]; then
        echo "Already have $FILE"
    else
        echo "Linking $FILE to cichlid positive control"
        cd expected/
        ln -s ../cichlid_control.known.tsv $ACC.known.tsv
        cd ..
    fi
done

# There are six negative controls:
for ACC in `grep "negative" PRJNA638011.tsv | cut -f 1`; do
    FILE=expected/$ACC.known.tsv
    if [ -f $FILE ]; then
        echo "Already have $FILE"
    else
        echo "Linking $FILE to negative control"
        cd expected/
        ln -s ../negative_control.known.tsv $ACC.known.tsv
        cd ..
    fi
done
