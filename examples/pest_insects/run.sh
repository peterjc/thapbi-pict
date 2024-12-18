#!/bin/bash

#The following will be used on SLURM via sbatch:
#======================================================
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=2:30:00
#SBATCH --job-name=pest_insects
#======================================================

set -euo pipefail

echo "NOTE: Expected first time run time is about two hours,"
echo "about 15 minutes from the merged reads cache."
echo "Repeat runs take a few minutes to regenerate reports."
echo

mkdir -p references/ intermediate/ summary/

function import_marker {
    # Takes arguments via variable names
    if [ ! -f merged_arthropoda_rdp_species.fasta ]; then
        echo "Unzipping and de-duplicating merged_arthropoda_rdp_species.fa.gz"
        cat merged_arthropoda_rdp_species.fa.gz | gunzip | ../../scripts/make_nr.py /dev/stdin > merged_arthropoda_rdp_species.fasta
    fi
    if [ ! -f references/$NAME.fasta ]; then
        echo "Trimming reference FASTA for $NAME"
        RIGHT_RC=$(python -c "from Bio.Seq import reverse_complement as rc; print(rc('$RIGHT'))")
        cutadapt --discard-untrimmed --quiet \
            -g $LEFT...$RIGHT_RC -o references/$NAME.fasta \
            merged_arthropoda_rdp_species.fasta
    fi
    echo "Importing $NAME into DB"
    thapbi_pict import -d references/merged_arthropoda.sqlite \
        -i references/${NAME}.fasta -x -s ";" \
        -k $NAME --left $LEFT --right $RIGHT
}

if [ ! -f references/merged_arthropoda.sqlite ]; then
    echo =================
    echo Creating database
    echo =================

    # Quoting hemiptera_metabarcoding.md from their GitHub,
    #
    # Name                    Illumina overhang adapter           Primer sequences
    # Sterno18S_F2_tail  ACACTCTTTCCCTACACGACGCTCTTCCGATCT  ATGCATGTCTCAGTGCAAG
    # Sterno18S_R1_tail  GACTGGAGTTCAGACGTGTGCTCTTCCGATC    TCGACAGTTGATAAGGCAGAC
    # Sterno12S_F2_tail  ACACTCTTTCCCTACACGACGCTCTTCCGATCT  CAYCTTGACYTAACAT
    # Sterno12S_R2_tail  GACTGGAGTTCAGACGTGTGCTCTTCCGATC    TAAAYYAGGATTAGATACCC
    # SternoCOI_F1_tail  ACACTCTTTCCCTACACGACGCTCTTCCGATCT  ATTGGWGGWTTYGGAAAYTG
    # SternoCOI_R1_tail  GACTGGAGTTCAGACGTGTGCTCTTCCGATC    TATRAARTTRATWGCTCCTA
    #

    #ILL_LEFT=ACACTCTTTCCCTACACGACGCTCTTCCGATCT
    #ILL_RIGHT=GACTGGAGTTCAGACGTGTGCTCTTCCGATC

    # 18S - Sterno18S_F/Sterno18S_R
    NAME=18S
    LEFT=ATGCATGTCTCAGTGCAAG
    RIGHT=TCGACAGTTGATAAGGCAGAC
    import_marker # calls function defined above

    # 12S - Sterno12S_F/Sterno12S_R
    NAME=12S
    LEFT=CAYCTTGACYTAACAT
    RIGHT=TAAAYYAGGATTAGATACCC
    import_marker # calls function defined above

    # COI - SternoCOI_F/SternoCOI_R
    NAME=COI
    LEFT=ATTGGWGGWTTYGGAAAYTG
    RIGHT=TATRAARTTRATWGCTCCTA
    import_marker # calls function defined above

fi

echo ================
echo Running pipeline
echo ================

thapbi_pict pipeline --cpu 4 \
    -d references/merged_arthropoda.sqlite \
    --synthetic "" -a 50 -f 0.0001 \
    -i raw_data/ expected/ \
    --merged-cache tmp_merged/ -s intermediate/ \
    -t metadata.tsv -c 3,4,2 -x 1 \
    -o summary/

echo ==============
echo Making figures
echo ==============

# First redraw with original data
if [ -f figure3original.tsv ]; then
    echo "Recreating original Figure 3 with Batovska et al. (2021) data"
    ./recreate_figure3.py -i figure3original.tsv -o figure3original.png
    ./recreate_figure3.py -i figure3original.tsv -o figure3original.pdf
else
    # We provide the TSV file, so shouldn't need to do this...
    echo "Missing figure3.tsv data file. Using RStudio run Batovska et al."
    echo "(2021)'s R code to generate Figure 3, then run our figure3.R script"
    echo "to export the figure's data."
fi

# Now draw with our data
./figure3reproduction.py
./recreate_figure3.py -i figure3reproduction.tsv -o figure3reproduction.png
./recreate_figure3.py -i figure3reproduction.tsv -o figure3reproduction.pdf

echo ====
echo Done
echo ====
