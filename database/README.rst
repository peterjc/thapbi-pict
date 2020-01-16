Default ITS1 database
=====================

We include a default ITS1 database with THAPBI PICT using the following three
sets of sequences (subject to taxonomy filtering):

- Curated *Phytophthora* ITS1 sequences (at species level) from file
  ``database/Phytophthora_ITS1_curated.fasta``. This is sorted by sequence.
  Examples sharing the same trimmed sequence are recorded as a single record
  with multiple semi-colon separated *accession genus species* entries.
  Sadly many published sequences are missing the first 32bp which our default
  primers would amplify, namely ``TTTCCGTAGGTGAACCTGCGGAAGGATCATTA``. This
  has been assumed and added to the FASTA in lower case as needed - at the
  time of writing no exceptions are known.

- NCBI *Peronosporales* (including *Phytophthora*) at genus level,
  file ``database/2019-04-03-ITS_Peronosporales_16394.fasta`` with 16394
  entries created with an NCBI Entrez search run on 2019-04-16::

      ((internal AND transcribed AND spacer) OR its1) AND
      150:10000[sequence length] AND Peronosporales[organism]

- Observed ITS1 sequences from single isolate positive controls run on a MiSeq
  plate via ``thapbi_pict prepare-reads`` with default settings (plate level
  minimum abundance was 545, but in anycase import minimum default is 1000).

- Four G-BLOCKS synthetic controls in file ``database/controls.fasta``

- This used the NCBI taxonomy as of 2019-11-01, which means we rejected some
  of the curated FASTA file entries or just used them at genus level..

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
