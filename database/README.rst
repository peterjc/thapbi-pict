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
  4466 entries in file ``database/2021-11-22-ITS_Oomycota_w32.fasta``
  extracted from 33064 downloaded from an NCBI Entrez search on 2021-11-22::

      (Peronosporales[organism] OR Pythiales[organism])
      AND ((internal AND transcribed AND spacer) OR its1)
      AND 150:10000[sequence length]

  Any obvious left and right primers were removed, and the sequence trimmed to
  the expected 32bp leader as follows::

      $ cutadapt -a GYRGGGACGAAAGTCYYTGC 2021-11-22-ITS_Oomycota_33064.fasta \
        --quiet | cutadapt -g GAAGGTGAAGTCGTAACAAGG --quiet /dev/stdin \
        | cutadapt -g TTTCCGTAGGTGAACCTGCGGAAGGATCATTA -O 30 --action retain \
        --discard-untrimmed --quiet /dev/stdin -o 2021-11-22-ITS_Oomycota_w32.fasta

  Note the import command discards uncultured entries.

- Observed ITS1 sequence files ``single_isolates/*.fasta`` from single isolate
  positive controls run on MiSeq plates. Created using
  ``thapbi_pict prepare-reads`` and ``thapbi_pict curated-seq`` with default
  settings (default minimum abundance 1000 was well above the plate level
  minimum abundance set via negative controls).

Additionally:

- Four G-BLOCKS synthetic controls in file ``database/controls.fasta``

- This used the NCBI taxonomy as of 2021-10-01, which means a handful of
  unrecognised species are imported at just genus level.

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
