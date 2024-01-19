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
  time of writing exceptions are rare.

- Curated *Nothophytophthora* ITS1 sequences (mostly at species level) from
  file ``Nothophytophthora_ITS1_curated.fasta``. Again, most have been
  extended with the assumed 32bp leader.

- Curated other Peronosporales ITS1 sequences including *Peronosclerospora*
  from file ``Peronosporales_ITS1_curated.fasta``, currently imported only to
  genus. Most have been extended with the assumed 32bp leader, or in some
  cases the final dozen bases.

- NCBI Peronosporales (including *Phytophthora*) or Pythiales at genus level,
  in two files extracted from tens of thousands of candidate sequences in file
  ``Oomycota_ITS1_search.fasta`` (not under version control) downloaded from
  an NCBI Entrez search via script ``Oomycota_ITS1_search.sh``::

      (Peronosporales[organism] OR Pythiales[organism])
      AND ((internal AND transcribed AND spacer) OR its1)
      AND 150:10000[sequence length]

  The search results are filtered to extract the amplicons with the expected
  32bp leader (allowing for the leading T or TT to be missing) with cutadapt
  as file ``Oomycota_ITS1_w32.fasta`` (under version control).

  Then script ``../scripts/missed_refs.py`` is used to catch useful references
  without the typical 32bp leader in full as file ``Oomycota_ITS1_obs.fasta``
  (under version control), where observed in at least five of our samples
  (using file ``unknowns.fasta`` which is not under version controls but is
  generated from our sequencing data using script ``../scripts/unknowns.py``).

  Note the import command discards uncultured entries.

- Observed ITS1 sequence files ``single_isolates/*.fasta`` with at least
  abundance 1000 (well above plate negative controls) from single isolate
  positive controls run on MiSeq plates.

Additionally:

- Four G-BLOCKS synthetic controls in file ``database/controls.fasta``

- This used the NCBI taxonomy recorded in the script, which typically means
  a handful of unrecognised curated entries may be imported without an NCBI
  taxid.

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
