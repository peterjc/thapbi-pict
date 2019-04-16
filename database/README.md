Database of ITS1 sequences for use as a molecular barcode.

Currently this folder only contains ITS1_DB.sql (plain text SQL
dump from sqlite3), from which we generate the binary DB using:

```bash
$ sqlite3 ITS1_DB.sqlite < ITS1_DB.sql
```

This database was created using the following three sets of
sequences:

 - Curated Phytophthora ITS1 sequences (at species level) from
   file ``legacy/Phytophthora_ITS_database_v0.005.fasta`` (in
   a subdirectory within the source code repository).

 - NCBI Peronosporales (including Phytophthora) at genus level,
   file ``2019-04-03-ITS_Peronosporales_16394.fasta`` with 16394
   entries created with an NCBI Entrez search run on 2019-04-16:

```
((internal AND transcribed AND spacer) OR its1) AND
150:10000[sequence length] AND Peronosporales[organism]
```

 - Observed ITS1 sequences from single isolate positive controls
   run on a MiSeq plate via ``thapbi_pict prepare-reads`` with
   default settings (plate level minimum abundance was 545,
   but in anycase import minimum default was 1000 copies).

 - This used the NCBI taxonomy as of 2019-01-01, which means
   we rejected some of the curated FASTA file entries.

The DB was created via the ``build_ITS1_DB.sh`` script.
