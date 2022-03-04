#!/bin/bash
set -euo pipefail

echo "=========================="
echo "Downloading reads from ENA"
echo "=========================="
# We are ignoring most of the Illumina sequencing runs
# (m4A, m4b, and m4C), and all the Ion Torrent data too.
if [ ! -f raw_data/MD5SUM.txt ]; then
    echo "ERROR: Missing raw_data/MD5SUM.txt"
    false
fi
# Only want the specific Illumina plates (using library prefix)
for ACC in `grep "Illumina MiSeq\tm6-" PRJNA305924.tsv | cut -f 1`; do
    # echo "Downloading $ACC"
    # Column 6 should have two URLs (R1 and R2), semi-colon separated:
    for URL in `grep ^$ACC PRJNA305924.tsv | cut -f 6 | sed "s/;/ /g"` ; do
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
done

# BioMock mixtures
mkdir -p expected
echo "Setting up positive controls"
for EXPT in BioMockStds BioMock SynMock; do
    for ACC in `grep "Illumina MiSeq\tm6-" PRJNA305924.tsv | grep "\t$EXPT$" | cut -f 1`; do
        FILE=expected/$ACC.known.tsv
        if [ -f $FILE ]; then
            echo "Already have $FILE"
        else
            echo "Linking $FILE to $EXPT mixture"
            cd expected/
            ln -s ../$EXPT.known.tsv $ACC.known.tsv
            cd ..
        fi
    done
done

echo "=========="
echo "Setup done"
echo "=========="
