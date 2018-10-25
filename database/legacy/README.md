# Legacy databases

These are FASTA format files of *Phytophthora* ITS1 sequences from earlier
in the THAPBI Phyto-Threats project.

# database.fasta circa 2016

File ``database.fasta`` was originally checked into the GitHub repository
https://github.com/widdowquinn/THAPBI-pycits/ ago on 24 October 2016 as the
following commit:

https://github.com/widdowquinn/THAPBI-pycits/commit/05502ddb15b7172c58f495d7df81f244e540925c

This is a FASTA file of 153 sequences, ranging in length from 163bp to 838bp
(i.e. some are much longer than the core ITS1 region of our interest).

Peter Cock understands this dataset was originally compiled by Santiago Garcia
(Forest Research).

The entry names like ``1Phytophthora_infestans_CBS36651`` take the format of
the clade name (here ``1``), species (here *Phytophthora infestans*) and an
identifier (here voucher number ``CBS36651`` whose 795bp nucleotide sequence
is available on the NCBI as accession ``HQ643247.1``).

# Phytophthora_ITS_database_v0.004.fasta from Feb 2018

File ``Phytophthora_ITS_database_v0.004.fasta`` taken from here, with the
spelling of *Phytophthora* corrected), created 28 Feb 2018:

https://github.com/peterthorpe5/public_scripts/blob/master/metapy/data/Phytophora_ITS_database_v0.004.fasta

This is a FASTA file of 170 sequences, ranging in length from 136bp to 203bp
(i.e. trimmed to the core ITS1 region of out interest).

# Phytophthora_ITS_database_v0.005.fasta from May 2018

File ``Phytophthora_ITS_database_v0.005.fasta`` taken from here, with the
spelling of *Phytophthora* corrected), created 30 May 2018:

https://github.com/peterthorpe5/public_scripts/commits/master/metapy/data/Phytophora_ITS_database_v0.005.fasta

This is a FASTA file of 176 sequences, ranging in length from 136bp to 203bp
(i.e. trimmed to the core ITS1 region of out interest).

This change added six sequences (but also changed the order so this cannot
be shown easily using ``diff``), which included four control sequences:

* 7a_Phytophthora_alni_subsp._alni_DQ335646
* 7b_Phytophthora_cinnamomi_KU726740
* Control_05 C1
* Control_14 C2
* Control_22 C3
* Control_42 C4
