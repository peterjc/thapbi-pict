#!/bin/bash
set -euo pipefail

mkdir -p intermediate/ summary/ expected/

for ACC in `grep ^SRR PRJNA638011.tsv | cut -f 1`; do
    # echo "Downloading $ACC"
    # Column 4 should have two URLs (R1 and R2), semi-colon separated:
    for URL in `grep ^$ACC PRJNA638011.tsv | cut -f 3 | sed "s/;/ /g"`; do
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

# There are four positive controls:
echo "Setting up positive controls"
for ACC in `grep "positive" PRJNA638011.tsv | cut -f 1`; do
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

# There are six "negative" and eight "blank" controls:
echo "Setting up negative and blank controls"
for ACC in `grep -E "(negative|blank)" PRJNA638011.tsv | cut -f 1`; do
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

# There are nine "New Lake" samples (plus 5 controls done above):
echo "Setting up expected New Lake classification"
for ACC in `grep "New Lake" PRJNA638011.tsv | cut -f 1`; do
    FILE=expected/$ACC.known.tsv
    if [ -f $FILE ]; then
        echo "Already have $FILE"
    else
        echo "Linking $FILE to cichlid positive control"
        cd expected/
        ln -s ../new_lake.known.tsv $ACC.known.tsv
        cd ..
    fi
done

echo "Setting up expected Middle Lake classification"
for ACC in `grep -E "(STX|MCE)" PRJNA638011.tsv | cut -f 1`; do
    FILE=expected/$ACC.known.tsv
    if [ -f $FILE ]; then
    echo "Already have $FILE"
    else
        echo "Linking $FILE to cichlid positive control"
    cd expected/
    ln -s ../middle_lake.known.tsv $ACC.known.tsv
    cd ..
    fi
done

echo "Setup done for drained ponds example"
