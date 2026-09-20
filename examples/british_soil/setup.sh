#!/bin/bash
set -eup pipeline

mkdir -p intermediate/ tmp_merged/ summary/ expected/

# ITS1 = PRJNA690943, COI-5P = PRJNA691575
for ACC in `grep ITS PRJNA690943.tsv | cut -f 1`; do
    # echo "Downloading $ACC"
    # Column 5 should have two URLs (R1 and R2), semi-colon separated:
    for URL in `grep ^$ACC PRJNA690943.tsv | cut -f 5 | sed "s/;/ /g"`; do
        FILE=raw_data/${URL##*/}
        # Avoiding leaving partial FASTQ if wget is interrupted
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

# There are three positive controls:
for ACC in `grep "Control Phytophthora" PRJNA690943.tsv | cut -f 1`; do
    # echo "Setting up expected classification"
    FILE=expected/$ACC.known.tsv
    if [ -f $FILE ]; then
        echo "Already have $FILE"
    else
        echo "Linking $FILE to mock community control"
        cd expected/
        ln -s ../mock_community.known.tsv $ACC.known.tsv
        cd ..
    fi
done

# There are two negative controls:
for ACC in `grep "Negative" PRJNA690943.tsv | cut -f 1`; do
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
