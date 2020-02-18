#!/bin/bash
set -eup pipeline

for LIBRARY in V4 V8V9 ; do
    echo "==========="
    echo $LIBRARY
    echo "==========="
    mkdir -p $LIBRARY/expected
    if [ ! -f $LIBRARY/raw_data/MD5SUM.txt ]; then
        echo "ERROR: Missing $LIBRARY/raw_data/MD5SUM.txt"
        false
    fi
    for ACC in `grep ^SRR PRJNA314977.txt | grep $LIBRARY | cut -f 1`; do
        # echo "Downloading $ACC"
        # Column 7 should have two URLs (R1 and R2), semi-colon separated:
        for URL in `grep ^$ACC PRJNA314977.txt | cut -f 7 | sed "s/;/ /g"`; do
            NAME=${URL##*/}
            FILE=$LIBRARY/raw_data/$NAME
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
        FILE=$LIBRARY/expected/$ACC.known.tsv
    if [ -f $FILE ]; then
        echo "Already have $FILE"
    elif grep $ACC PRJNA314977.txt | grep --quiet "_MC1[1-3]"; then
        # i.e. Samples MC11, MC12, MC13
        echo "Linking $FILE to mock freshwater community control"
        cd $LIBRARY/expected/
        ln -s ../../mock_freshwater.known.tsv $ACC.known.tsv
        cd ../..
    elif grep $ACC PRJNA314977.txt | grep --quiet "_MC7[1-3]"; then
        # i.e. Samples MC71, MC72, MC73
        echo "Linking $FILE to mock marine community control"
        cd $LIBRARY/expected/
        ln -s ../../mock_marine.known.tsv $ACC.known.tsv
        cd ../..
    elif grep $ACC PRJNA314977.txt | grep --quiet "synthetic metagenome"; then
        echo "Linking $FILE to mock community control"
        cd $LIBRARY/expected/
        ln -s ../../mock_community.known.tsv $ACC.known.tsv
        cd ../..
    fi
    done
done
