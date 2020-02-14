#!/bin/bash
set -eup pipeline

if [ ! -f raw_data/MD5SUM.txt ]; then
    echo "ERROR: Missing raw_data/MD5SUM.txt"
    false
fi
for ACC in `grep ^SRR metadata.tsv | cut -f 1`; do
    # echo "Downloading $ACC"
    # Column 6 should have two URLs (R1 and R2), semi-colon separated:
    for URL in `grep ^$ACC PRJNA417859.txt | cut -f 6 | sed "s/;/ /g"`; do
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
