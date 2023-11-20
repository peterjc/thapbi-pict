#!/bin/bash
VERSION=`thapbi_pict -v | sed "s/THAPBI PICT //g"`
echo "Using THAPBI PICT $VERSION"
set -euo pipefail
if ! command -v efetch &> /dev/null; then
    echo "ERROR: NCBI efetch could not be found, try 'conda install entrez-direct'"
    exit 1
fi
if ! command -v esearch &> /dev/null; then
    echo "ERROR: NCBI esearch could not be found, try 'conda install entrez-direct'"
    exit 1
fi

export TMP=${TMP:-/tmp/thapbi_pict}/fetch
rm -rf $TMP
mkdir -p $TMP

echo "Starting NCBI search..."
esearch -db nuccore -sort accession \
     -query "(Peronosporales[organism] OR Pythiales[organism])\
     AND ((internal AND transcribed AND spacer) OR its1)\
     AND 150:10000[sequence length]" > $TMP/search.xml

echo "Fetching NCBI matches (likely to take over an hour)..."
efetch -format fasta < $TMP/search.xml > $TMP/search.fasta

COUNT=`grep -c "^>" $TMP/search.fasta`
if [ `grep -oh "<Count>\d*</Count>" /tmp/search.xml` != "<Count>$COUNT</Count>" ]; then
    echo "ERROR: FASTA file had $COUNT entries, but search said:"
    grep -oh "<Count>\d*</Count>" /tmp/search.xml
    exit 1
fi

echo "Downloaded $COUNT ITS1 sequences."
mv $TMP/search.fasta Oomycota_ITS1_search.fasta
echo

echo "Selecting and trimming those with expected 32bp leader..."
cutadapt -a GYRGGGACGAAAGTCYYTGC Oomycota_ITS1_search.fasta \
  --discard-untrimmed -e 0.2 --quiet \
  | sed "s/^TTCCGTAGGTGAAC/tTTCCGTAGGTGAAC/" \
  | sed  "s/^TCCGTAGGTGAAC/ttTCCGTAGGTGAAC/"  \
  | cutadapt -g GAAGGTGAAGTCGTAACAAGG --quiet /dev/stdin \
  | cutadapt -g TTTCCGTAGGTGAACCTGCGGAAGGATCATTA -O 32 --action retain \
  --discard-untrimmed -M 450 --quiet /dev/stdin \
  -o Oomycota_ITS1_w32.fasta

echo "Checking for observed NCBI entries without the expected 32bp leader..."
if [ ! -f unknowns.fasta ]; then
    echo "ERROR: Missing unknowns.fasta, use:"
    echo "../scripts/unknowns.py -i thapbi-pict.ITS1.reads.identity.tsv \\"
    echo "-a 1000 -s 5 -o unknowns.fasta"
    exit 1
fi
../scripts/missed_refs.py -i unknowns.fasta \
  -f Oomycota_ITS1_search.fasta \
  -x Oomycota_ITS1_w32.fasta \
  -o Oomycota_ITS1_obs.fasta

echo "Done - review and commit Oomycota_ITS1_w32.fasta & Oomycota_ITS1_obs.fasta"
