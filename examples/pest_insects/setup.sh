#!/bin/bash
set -euo pipefail

mkdir -p raw_data/ expected/ tmp_merged/ intermediate/ summary/
if [ ! -f raw_data/MD5SUM.txt ]; then
    echo "ERROR: Missing raw_data/MD5SUM.txt"
    #false
fi

# This FASTA file is also available via https://zenodo.org/record/3557020
wget -nc "https://github.com/alexpiper/HemipteraMetabarcodingMS/raw/440e08ef8c0ec357fe130a2116144bc1453dc00a/reference/merged_arthropoda_rdp_species.fa.gz"

echo "=========================="
echo "Downloading reads from ENA"
echo "=========================="
# The original FASTQ files are also available in data.rar from
# https://zenodo.org/record/5171623
for ACC in $(grep ^SRR PRJNA716058.tsv | cut -f 1); do
    # echo "Downloading $ACC"
    # Column 6 should have two URLs (R1 and R2), semi-colon separated:
    for URL in $(grep ^$ACC PRJNA716058.tsv | cut -f 6 | sed "s/;/ /g"); do
        NAME=${URL##*/}
        FILE=raw_data/$NAME
        # Avoiding leaving partial FASTQ if wget is interrupted
        rm -rf "$FILE.tmp"
        if [ -f "$FILE" ]; then
            echo "Already have $FILE"
        else
            echo "Downloading $FILE"
            wget -O "${FILE}.tmp" "ftp://$URL"
            mv "${FILE}.tmp" "${FILE}"
        fi
    done
done

echo "=================================="
echo "Setting up expected classification"
echo "=================================="
for POOL in {1..5}; do
    for ACC in $(grep Pool-${POOL}_ PRJNA716058.tsv | cut -f 1); do
        FILE=expected/$ACC.known.tsv
        if [ -f "$FILE" ]; then
            echo "Already have $FILE"
        else
            echo "Linking $FILE to Pool $POOL mock community"
            cd expected/
            ln -s ../mock_community_${POOL}.known.tsv $ACC.known.tsv
            cd ..
        fi
    done
done

echo "=========="
echo "Setup done"
echo "=========="
