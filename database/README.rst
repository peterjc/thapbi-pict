Default ITS1 database
=====================

We include a default ITS1 database with THAPBI PICT using the following four
sets of biological sequences (subject to taxonomy filtering):

- Curated *Phytophthora* ITS1 sequences (mostly at species level) from file
  ``database/Phytophthora_ITS1_curated.fasta``. This is sorted by sequence.
  Examples sharing the same trimmed sequence are recorded as a single record
  with multiple semi-colon separated *accession genus species* entries.
  Sadly many published sequences are missing the first 32bp which our default
  primers would amplify, namely ``TTTCCGTAGGTGAACCTGCGGAAGGATCATTA``. This
  has been assumed and added to the FASTA in lower case as needed - at the
  time of writing no exceptions are known.

- Curated *Nothophytophthora* ITS1 sequences (mostly at species level) from
  file ``Nothophytophthora_ITS1_curated.fasta``. Again, most have been
  extended with the assumed 32bp leader.
- NCBI Peronosporales (including *Phytophthora*) or Pythiales at genus level,
  in two files extracted from 34111 sequences downloaded from an NCBI Entrez
  search on 2022-05-07::

      (Peronosporales[organism] OR Pythiales[organism])
      AND ((internal AND transcribed AND spacer) OR its1)
      AND 150:10000[sequence length]

  Require and remove the right primer, remove any left primer if present, and
  require and trim the start to the expected 32bp leader as follows::

      $ cutadapt -a GYRGGGACGAAAGTCYYTGC 2022-07-05_ITS1_Oomycota_34111.fasta \
        --discard-untrimmed -e 0.2 --quiet \
        | cutadapt -g GAAGGTGAAGTCGTAACAAGG --quiet /dev/stdin \
        | cutadapt -g TTTCCGTAGGTGAACCTGCGGAAGGATCATTA -O 30 --action retain \
        --discard-untrimmed --quiet /dev/stdin -o 2022-07-05_ITS1_Oomycota_w32.fasta

  Then to catch useful references without the typical 32bp leader in full,
  where observed in at least five of our samples::

      $ ../scripts/unknowns.py -i thapbi-pict.ITS1.reads.identity.tsv \
        -a 1000 -s 5 -o unknowns.fasta
      $ ../scripts/missed_refs.py -i unknowns.fasta \
        -f 2022-07-05_ITS1_Oomycota_34111.fasta \
        -x 2022-07-05_ITS1_Oomycota_w32.fasta \
        -o 2022-07-05_ITS1_Oomycota_obs.fasta

  Note the import command discards uncultured entries.

- Observed ITS1 sequence files ``single_isolates/*.fasta`` from single isolate
  positive controls run on MiSeq plates. Created using
  ``thapbi_pict prepare-reads`` and ``thapbi_pict curated-seq`` with default
  settings (default minimum abundance 1000 was well above the plate level
  minimum abundance set via negative controls).

Additionally:

- Four G-BLOCKS synthetic controls in file ``database/controls.fasta``

- This used the NCBI taxonomy as of 2022-07-01, which means a handful of
  unrecognised curated entries are imported without an NCBI taxid.

The database is created with the ``database/build_ITS1_DB.sh`` script:

.. code:: console

    $ cd database/
    $ ./build_ITS1_DB.sh
    $ chmod a-w ITS1_DB.sqlite
    $ cp ITS1_DB.sqlite ../thapbi_pict/ITS1_DB.sqlite

We include the binary SQLite3 file as ``thapbi_pict/ITS1_DB.sqlite`` in the
software releases, but that binary SQLite3 database file itself is not under
version control. Instead, we track a plain text FASTA dump of the database as
``database/ITS1_DB.fasta`` giving a meaningful change history (we initially
tracked this as a plain text SQL dump, ``database/ITS1_DB.sql``).
