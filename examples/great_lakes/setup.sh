#!/bin/bash
set -eup pipeline

mkdir -p intermediate/ summary/

if [ ! -f raw_data/MD5SUM.txt ]; then
    echo "ERROR: Missing raw_data/MD5SUM.txt"
    false
fi
for NAME in MOL16S SPH16S; do
    mkdir -p expected/$NAME/
    for ACC in `grep Blank PRJNA379165.tsv | grep $NAME | cut -f 1`; do
        FILE=expected/$NAME/$ACC.known.tsv
        if [ -f $FILE ]; then
            echo "Already have $FILE"
        else
            echo "Linking $FILE to negative control"
            cd expected/$NAME/
            ln -s ../negative_control.known.tsv $ACC.known.tsv
            cd ../..
        fi
    done
    for ACC in `grep Mock PRJNA379165.tsv | grep -v $NAME | cut -f 1`; do
        FILE=expected/$NAME/$ACC.known.tsv
        if [ -f $FILE ]; then
            echo "Already have $FILE"
        else
            echo "Linking $FILE to negative control"
            cd expected/$NAME/
            ln -s ../negative_control.known.tsv $ACC.known.tsv
            cd ../..
        fi
    done
    for ACC in `grep Mock PRJNA379165.tsv | grep $NAME | cut -f 1`; do
        FILE=expected/$NAME/$ACC.known.tsv
        if [ -f $FILE ]; then
            echo "Already have $FILE"
        else
            echo "Linking $FILE to mock community"
            cd expected/$NAME/
            ln -s ../mock_community.$NAME.known.tsv $ACC.known.tsv
            cd ../..
        fi
    done
done

for ACC in `grep ^SRR PRJNA379165.tsv| cut -f 1`; do
    # echo "Downloading $ACC"
    # Column 7 should have two URLs (R1 and R2), semi-colon separated:
    for URL in `grep ^$ACC PRJNA379165.tsv | cut -f 7 | sed "s/;/ /g"`; do
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
