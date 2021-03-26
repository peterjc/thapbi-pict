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
- NCBI Oomycota (including *Phytophthora*) at genus level, 4907 entries in
  file ``database/2021-01-28-ITS_Oomycota_w32.fasta`` trimmed to start at the
  expected 32bp leader, and any obvious right primer removed (which the import
  would do anyway - but this gives a smaller more useful intermediate file),
  using::

      $ cutadapt -g TTTCCGTAGGTGAACCTGCGGAAGGATCATTA -O 30 --action retain \
      --discard-untrimmed 2021-02-23-ITS_Oomycota_36025.fasta \
      | cutadapt -a GYRGGGACGAAAGTCYYTGC --fasta /dev/stdin \
      -o 2021-02-23-ITS_Oomycota_w32.fasta

  Started from 36025 downloaded from an NCBI Entrez search on 2021-02-23::

      ((internal AND transcribed AND spacer) OR its1) AND
      150:10000[sequence length] AND Oomycota[organism]

  Note the NCBI import command discards uncultured entries.

- Observed ITS1 sequence files ``single_isolates/*.fasta`` from single isolate
  positive controls run on MiSeq plates. Created using
  ``thapbi_pict prepare-reads`` and ``thapbi_pict curated-seq`` with default
  settings (default minimum abundance 1000 was well above the plate level
  minimum abundance set via negative controls).

Additionally:

- Four G-BLOCKS synthetic controls in file ``database/controls.fasta``

- This used the NCBI taxonomy as of 2021-01-01, which means we rejected some
  of the curated FASTA file entries, or just used them at genus level.

The database is created with the ``database/build_ITS1_DB.sh`` script:

.. code:: console

    $ cd database/
    $ ./build_ITS1_DB.sh
    $ chmod a-w ITS1_DB.sqlite
    $ cp ITS1_DB.sqlite ../thapbi_pict/ITS1_DB.sqlite

We include the binary SQLite3 file as ``thapbi_pict/ITS1_DB.sqlite`` in the
software releases, but that binary SQLite3 database file itself is not under
version control. Instead, we track a plain text dump of the database as
``database/ITS1_DB.sql`` giving a meaningful change history.

As an alternative to regenerating the database, the text dump can be converted
into a binary SQLite3 database file as follows:

.. code:: console

   $ sqlite3 thapbi_pict/ITS1_DB.sqlite < database/ITS1_DB.sql
   $ chmod a-w thapbi_pict/ITS1_DB.sqlite
