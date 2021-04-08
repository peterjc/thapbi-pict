#!/bin/bash
set -eup pipeline

mkdir -p expected/ intermediate/ summary/
for LIBRARY in AL1 AL2; do
    echo "==="
    echo $LIBRARY
    echo "==="
    mkdir -p expected/$LIBRARY/
    if [ ! -f raw_data/$LIBRARY/MD5SUM.txt ]; then
        echo "ERROR: Missing $LIBRARY/raw_data/MD5SUM.txt"
        false
    fi
    for ACC in `grep ^SRR metadata_$LIBRARY.tsv | cut -f 1`; do
        # echo "Downloading $ACC"
        # Column 6 should have two URLs (R1 and R2), semi-colon separated:
        for URL in `grep ^$ACC PRJNA377530.tsv | cut -f 6 | sed "s/;/ /g"`; do
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
        elif grep $ACC metadata_$LIBRARY.tsv | grep --quiet negative; then
            echo "Linking $FILE to negative control"
            cd expected/$LIBRARY/
            ln -s ../../negative_control.known.tsv $ACC.known.tsv
            cd ../..
        else
            echo "Linking $FILE to mock community"
            cd expected/$LIBRARY/
            ln -s ../../mock_community.known.tsv $ACC.known.tsv
            cd ../..
        fi
    done
done
