#!/bin/bash
set -eup pipeline

if [ ! -f COI_430_bats.fasta ]; then
    echo "Fetching reference FASTA file of COI sequence for 430 bats,"
    echo "Supplmentary S2 from Walker et al. (2016)"
    wget -O COI_430_bats.fasta.tmp "https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0162342.s002&type=supplementary"
    # Replace underscores with spaces
    cat COI_430_bats.fasta.tmp | sed "s#_# #g" > COI_430_bats.fasta
    rm COI_430_bats.fasta.tmp
fi

mkdir -p expected/ intermediate/ summary/

if [ ! -f raw_data/MD5SUM.txt ]; then
    echo "ERROR: Missing raw_data/MD5SUM.txt"
    false
fi
for ACC in `grep ^SRR PRJNA574765.tsv | cut -f 1`; do
    # echo "Downloading $ACC"
    # Column 5 should have two URLs (R1 and R2), semi-colon separated:
    for URL in `grep ^$ACC PRJNA574765.tsv | cut -f 5 | sed "s/;/ /g"`; do
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
    # echo "Setting up expected classification"
    FILE=expected/$ACC.known.tsv
    if [ -f $FILE ]; then
        echo "Already have $FILE"
    else
        echo "Linking $FILE to mock community"
        cd expected/
        ln -s ../mock_community.known.tsv $ACC.known.tsv
        cd ..
    fi
done
