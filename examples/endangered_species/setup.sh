#!/bin/bash
set -eup pipeline

mkdir -p raw_data/ expected/ tmp_merged/ intermediate/ intermediate_pool/ summary/
if [ ! -f raw_download/MD5SUM.txt ]; then
    echo "ERROR: Missing raw_download/MD5SUM.txt"
    false
fi

echo "=========================="
echo "Downloading reads from ENA"
echo "=========================="
for ACC in `grep ERR PRJEB18620.tsv | cut -f 1`; do
    # echo "Downloading $ACC"
    SAMPLE=`grep ^$ACC PRJEB18620.tsv | cut -f 4`
    # Column 3 should have two URLs (R1 and R2), semi-colon separated:
    for URL in `grep ^$ACC PRJEB18620.tsv | cut -f 3 | sed "s/;/ /g"`; do
        NAME=${URL##*/}
        FILE=raw_download/$NAME
        if [ "EM_" == ${NAME:0:3} ]; then
            # Rename EM_* files as really *zipped* FASTQ, not *gzipped*
            FILE=raw_download/${NAME%.gz}.zip
        fi
        # Avoiding leaving partial FASTQ if wget is interupted
        rm -rf $FILE.tmp
        if [ -f $FILE ]; then
            echo "Already have $FILE"
        else
            echo "Downloading $FILE"
            wget -O "${FILE}.tmp" "$URL"
            mv "${FILE}.tmp" "${FILE}"
        fi
        if [ "EM_" == ${NAME:0:3} ]; then
            # The EM_* files are really *zipped* FASTQ, not *gzipped*
            # The filename does however already match the short sample
            # (and we're ignoring the long name used within the zip file)
            NEW=raw_data/$NAME
            if [ ! -f $NEW ]; then
                echo "Unzipping $FILE and instead gzipping to make $NEW"
                unzip -p "$FILE" | gzip > "${NEW}.tmp"
                mv "${NEW}.tmp" "${NEW}"
            fi
        else
            # The lab files are already gzipped, but have long names
            if [[ "$NAME" =~ "_R1_" ]]; then
                NEW=raw_data/${SAMPLE}_R1.fastq.gz
            elif [[ "$NAME" =~ "_R2_" ]]; then
                NEW=raw_data/${SAMPLE}_R2.fastq.gz
            else
                echo "ERROR: $NAME lacks R1 or R2"
                false
            fi
            if [ ! -f $NEW ]; then
                ln -s ../$FILE $NEW
            fi
        fi
    done
done

echo "=================================="
echo "Setting up expected classification"
echo "=================================="
cd expected/
# Samples 3 and 8 are traditional medicines, not mock communities
for SAMPLE in 1 2 4 5 6 7 9 10; do
    echo "Sample $SAMPLE"
    for LAB in {1..16}; do
        FILE="S${SAMPLE}_Lab_${LAB}.known.tsv"
        if [ ! -f "$FILE" ]; then
            echo "Linking $FILE to mock community"
            ln -s "S${SAMPLE}.template.tsv" "$FILE"
        fi
    done
done
cd ../

echo "=========="
echo "Setup done"
echo "=========="
