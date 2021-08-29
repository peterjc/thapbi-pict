#!/bin/bash
set -euo pipefail

mkdir -p raw_data/ expected/ tmp_merged/ intermediate/ intermediate_pool/ summary/
if [ ! -f raw_data/MD5SUM.txt ]; then
    echo "ERROR: Missing raw_data/MD5SUM.txt"
    #false
fi

echo "=========================="
echo "Downloading reads from ENA"
echo "=========================="
for ACC in `grep ^ERR PRJEB27581.tsv | cut -f 1`; do
    # echo "Downloading $ACC"
    # Column 4 should have two URLs (R1 and R2), semi-colon separated:
    for URL in `grep ^$ACC PRJEB27581.tsv | cut -f 4 | sed "s/;/ /g"` ; do
        NAME=${URL##*/}
        FILE=raw_data/$NAME
        # Avoiding leaving partial FASTQ if wget is interrupted
        rm -rf $FILE.tmp
        if [ -f $FILE ]; then
            echo "Already have $FILE"
        else
            echo "Downloading $FILE"
            wget -O "${FILE}.tmp" "ftp://$URL"
            mv "${FILE}.tmp" "${FILE}"
        fi
    done
    # Beware 1 vs I typo in ENA metadata for NF1-18Sr2b in PRJEB27581.tsv
    # The paper and our metadata.tsv uses the digit 1.
    for MARKER in NF1-18Sr2b SSUF04-SSUR22 D3Af-D3Br JB3-JB5GED; do
        mkdir -p expected/$MARKER/
        if [ ! -f expected/$MARKER/$ACC.known.tsv ]; then
            echo "Creating link for expected/$MARKER/$ACC.known.tsv"
            cd expected/$MARKER/
            if [[ `grep $ACC ../../metadata.tsv | cut -f 3` == "Blank" ]]; then
                ln -s ../../negative_control.known.tsv $ACC.known.tsv
            elif [[ `grep $ACC ../../metadata.tsv | cut -f 4` == "$MARKER" ]]; then
                ln -s ../../mock_community.known.tsv $ACC.known.tsv
            else
                # Should not contain this marker
                ln -s ../../negative_control.known.tsv $ACC.known.tsv
            fi
            cd ../..
        fi
    done
done

echo "=========="
echo "Setup done"
echo "=========="
