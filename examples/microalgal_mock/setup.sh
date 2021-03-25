#!/bin/bash
set -eup pipeline

mkdir -p expected/ summary/
for LIBRARY in V4 V8V9 ; do
    echo "==========="
    echo $LIBRARY
    echo "==========="
    mkdir -p expected/$LIBRARY/
    if [ ! -f raw_data/$LIBRARY/MD5SUM.txt ]; then
        echo "ERROR: Missing raw_data/$LIBRARY/MD5SUM.txt"
        false
    fi
    for ACC in `grep ^SRR PRJNA314977.tsv | grep $LIBRARY | cut -f 1`; do
        # echo "Downloading $ACC"
        # Column 7 should have two URLs (R1 and R2), semi-colon separated:
        for URL in `grep ^$ACC PRJNA314977.tsv | cut -f 7 | sed "s/;/ /g"`; do
            NAME=${URL##*/}
            FILE=raw_data/$LIBRARY/$NAME
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
        # echo "Setting up expected classification"
        FILE=expected/$LIBRARY/$ACC.known.tsv
        if [ -f $FILE ]; then
            echo "Already have $FILE"
        elif grep $ACC PRJNA314977.tsv | grep --quiet "_MC1[1-3]"; then
            # i.e. Samples MC11, MC12, MC13
            echo "Linking $FILE to mock freshwater community control"
            cd expected/$LIBRARY/
            ln -s ../../mock_freshwater.known.tsv $ACC.known.tsv
            cd ../..
        elif grep $ACC PRJNA314977.tsv | grep --quiet "_MC7[1-3]"; then
            # i.e. Samples MC71, MC72, MC73
            echo "Linking $FILE to mock marine community control"
            cd expected/$LIBRARY/
            ln -s ../../mock_marine.known.tsv $ACC.known.tsv
            cd ../..
        elif grep $ACC PRJNA314977.tsv | grep --quiet "synthetic metagenome"; then
            echo "Linking $FILE to mock community control"
            cd expected/$LIBRARY/
            ln -s ../../mock_community.known.tsv $ACC.known.tsv
            cd ../..
        fi
    done
done
