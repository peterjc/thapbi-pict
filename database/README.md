Database of ITS1 sequences for use as a molecular barcode.

Currently this folder only contains ITS1_DB_v0.1.0.sqlite.bz2
(compressed from 7MB file ITS1_DB_v0.1.0.sqlite) created using
thapbi_pict v0.0.15 using three sets of sequences.
    
 - Curated Phytophthora ITS1 sequences (at species level) from
   file ``legacy/Phytophthora_ITS_database_v0.005.fasta`` (in
   a subdirectory within the source code repository).
    
 - NCBI Peronosporaceae (including Phytophthora) at genus level,
   file ``2019-04-03-ITS_Peronosporaceae_8336.fasta`` with 8336
   entries created with and NCBI Entrez search run on 2019-04-03:

```
((internal AND transcribed AND spacer) OR its1) AND
150:800[sequence length] AND peronosporaceae[organism]
```
    
 - Observed ITS1 sequences from single isolate positive controls
   run on a MiSeq plate via ``thapbi_pict prepare-reads`` with
   default settings (plate level minimum abundance was 535,
   but in anycase import minimum default was 1000 copies).

 - This used the NCBI taxonomy as of 2019-01-01, which means
   we rejected some of the curated FASTA file entries.

The DB was created via the following script:
    

```bash
#!/bin/bash
VERSION=`thapbi_pict -v | sed "s/THAPBI PICT //g"`
echo "Using THAPBI PICT $VERSION"
set -euo pipefail
TAX=new_taxdump_2019-01-01
DB=L5-and-8336-Peronosporaceae-and-PosCtrl-$VERSION
rm -rf "$DB.sqlite" "$DB.fasta"
    
thapbi_pict load-tax -d "$DB.sqlite" -t "$TAX"
thapbi_pict legacy-import -d "$DB.sqlite" $HOME/repositories/thapbi-pict/database/legacy/Phytophthora_ITS_database_v0.005.fasta
thapbi_pict ncbi-import -d "$DB.sqlite" 2019-04-03-ITS_Peronosporaceae_8336.fasta -g
thapbi_pict seq-import -d "$DB.sqlite" thapbi20180709p1_MetaControls/prepared_reads_${VERSION}/*.fasta thapbi20180709p1_MetaControls/positive_controls/*.known.tsv
    
thapbi_pict dump -d "$DB.sqlite" -o "$DB.txt"
thapbi_pict dump -f fasta -d "$DB.sqlite" -o "$DB.fasta"
```
