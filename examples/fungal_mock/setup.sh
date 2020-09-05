#!/bin/bash
set -eup pipeline

for LIBRARY in amp_lib_one amp_lib_two; do
    echo "==========="
    echo $LIBRARY
    echo "==========="
    mkdir -p $LIBRARY/expected/ $LIBRARY/intermediate/ $LIBRARY/summary/
    if [ ! -f $LIBRARY/raw_data/MD5SUM.txt ]; then
        echo "ERROR: Missing $LIBRARY/raw_data/MD5SUM.txt"
        false
    fi
    for ACC in `grep ^SRR $LIBRARY/metadata.tsv | cut -f 1`; do
        # echo "Downloading $ACC"
        # Column 6 should have two URLs (R1 and R2), semi-colon separated:
        for URL in `grep ^$ACC PRJNA377530.txt | cut -f 6 | sed "s/;/ /g"`; do
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
        elif grep $ACC $LIBRARY/metadata.tsv | grep --quiet negative; then
            echo "Linking $FILE to negative control"
            cd $LIBRARY/expected/
            ln -s ../../negative_control.known.tsv $ACC.known.tsv
            cd ../..
        else
            echo "Linking $FILE to mock community"
            cd $LIBRARY/expected/
            ln -s ../../mock_community.known.tsv $ACC.known.tsv
            cd ../..
        fi
    done
done
