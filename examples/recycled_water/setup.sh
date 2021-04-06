#!/bin/bash
set -eup pipeline

mkdir -p intermediate/ summary/

if [ -f taxdmp_2019-12-01/names.dmp ]; then
    echo "Already have taxdmp_2019-12-01/"
else
    echo "Downloading NCBI taxonomy"
    wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2019-12-01.zip
    echo "Decompressing NCBI taxonomy"
    # Only need names.dmp and nodes.dmp
    unzip -n -d taxdmp_2019-12-01 taxdmp_2019-12-01.zip names.dmp nodes.dmp
fi


if [ ! -f raw_data/MD5SUM.txt ]; then
    echo "ERROR: Missing raw_data/MD5SUM.txt"
    false
fi
for ACC in `grep ^SRR PRJNA417859.tsv| cut -f 1`; do
    # echo "Downloading $ACC"
    # Column 6 should have two URLs (R1 and R2), semi-colon separated:
    for URL in `grep ^$ACC PRJNA417859.tsv | cut -f 6 | sed "s/;/ /g"`; do
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
