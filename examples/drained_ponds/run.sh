#!/bin/bash
set -eup pipeline

echo NOTE: Expected first time run time is under 5 minutes,
echo repeat runs take seconds just to regenerate reports
echo

if [ ! -f NCBI_12S.sqlite ]; then
    echo "Building 12S database from NCBI sequences"
    thapbi_pict import -d NCBI_12S.sqlite -i NCBI_12S.fasta -x -s ";"
fi

# Primers, quoting Muri et al:
#
#    In the first round of PCR, indexed primers targeting a 106 bp
#    region within the mitochondrial 12S gene were used (Riaz et al.
#    2011; Kelly et al. 2014).
#
# Quoting Kelly et al. 2014:
#
#    We amplified target samples using PCR primers designed by Riaz and
#    coauthors to amplify vertebrate-specific fragments from the
#    mitochondrial 12S rRNA gene (Riaz et al. 2011). A 106 bp fragment
#    from a variable region of the 12S rRNA gene was amplified with the
#    primers F1 (5′-ACTGGGATTAGATACCCC-3′) and R1 (5′- TAGAACAGGCTCCTCTAG-3′).

echo "Running the pipeline"
mkdir -p intermediate/ summary/
thapbi_pict pipeline -i raw_data/ expected/ -s intermediate/ \
            -o summary/ -r drained_ponds \
            --minlen 80 --maxlen 130 -a 50 -d NCBI_12S.sqlite \
            --left ACTGGGATTAGATACCCC --right TAGAACAGGCTCCTCTAG \
            -t metadata.tsv -x 1 -c 5,6,7,8,9,10,4,2,3

echo "Done"
