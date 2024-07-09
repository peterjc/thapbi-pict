#!/bin/bash
set -euo pipefail

echo "NOTE: Expected first time run time is about XX minutes,"
echo "repeat runs are a few seconds (regenerating reports etc)."

mkdir -p tmp_merged/ intermediate/ summary/

if [ ! -f 2020-10-22_release_1_rps10.fasta ]; then
    echo "Missing rps10 references, run setup.sh"
    false
fi

if [ ! -f pooled.sqlite ]; then
    echo =================
    echo Creating database
    echo =================

    # Same ITS1 primers, take a copy of the default DB
    cp ../../thapbi_pict/ITS1_DB.sqlite pooled.sqlite

    # The rps10 assay is a multiplex PCR reaction comprising two rps10 forward primers that differ
    # slightly in sequence but anneal to the same position in the tRNA-Phe gene (rps10_F1 and
    # rps10_F2) and seven rps10 reverse primers that differ slightly in sequence but anneal to the
    # same position in the tRNA-Arg gene (rps10_R1 through rps10_R7) (Table 1).

    # However, follow the authors in using IUPAC codes to approximate the primers:
    # https://github.com/grunwaldlab/rps10_barcode/blob/main/raw_data/primer_data.csv

    # Doing the left and right primer trimming separately:
    cutadapt --quiet -g GTTGGTTAGAGYARAAGACT 2020-10-22_release_1_rps10.fasta |
        cutadapt --quiet -a AGTTCRARTCTTTCTARRYAT /dev/stdin |
        ./rename_fasta.py > rps10.fasta
    thapbi_pict import -d pooled.sqlite \
        -i rps10.fasta -c obitools -x \
        -k rps10 --left GTTGGTTAGAGYARAAGACT --right ATRYYTAGAAAGAYTYGAACT

    # Note we are ignoring the samples for the alternative Felipe rps10 primers,
    # where the left primer is ~25bp later but right is same (overlapping product)
fi

echo ================
echo Running pipeline
echo ================

# Drop from -a 100 -f 0.001
mkdir -p intermediate/ summary/
thapbi_pict pipeline -d pooled.sqlite --synthetic '' -m 1s3g \
    -a 50 -f 0.0001 -i raw_data/ expected/ \
    --merged-cache tmp_merged/ \
    -s intermediate/ -o summary/ \
    -t metadata.tsv -x 1 -c 7,6,3,4
#-t PRJNA699663.tsv -x 1 -c 5,3,4

echo ====
echo Done
echo ====
