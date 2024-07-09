#!/bin/bash
set -euo pipefail

mkdir -p expected/ intermediate/ summary/

if [ ! -f raw_data/MD5SUM.txt ]; then
    echo "ERROR: Missing raw_data/MD5SUM.txt"
    false
fi

if [ ! -f 2020-10-22_release_1_rps10.fasta ]; then
    echo "Downloading rps10 references (need primer trimming)"
    # TODO - This was a pre-release with only 809 sequences,
    # actual release 1 was dated 2021-03-01, with 885 sequences
    wget "https://github.com/grunwaldlab/OomyceteDB/raw/master/website/2020-10-22_release_1_rps10.fasta"
fi

for ACC in $(grep ^SRR PRJNA699663.tsv | cut -f 1); do
    if grep ^$ACC metadata.tsv | grep -q "rps10_Felipe"; then
        echo "Ignoring $ACC as using alternative rps10 primers"
        # The "Felipe" rps10 primers differ in the left primer, giving a
        # smaller product. This makes it impossible to separate the amplicons
        # with an unanchored primer search - so since the paper ignores these
        # samples, we will too.
    else
        # Column 7 should have two URLs (R1 and R2), semi-colon separated:
        for URL in $(grep ^$ACC PRJNA699663.tsv | cut -f 7 | sed "s/;/ /g"); do
            NAME=${URL##*/}
            FILE=raw_data/$NAME
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
        # Quoting https://github.com/grunwaldlab/rps10_barcode/blob/main/06--mock_community.Rmd
        # Only `mock2` will be used in the paper. `mock1` was another, older mock community.
        if grep ^$ACC PRJNA699663.tsv | grep -q "mock community 2"; then
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
        fi
    fi
done
