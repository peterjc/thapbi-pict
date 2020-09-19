Reference database
==================

Introduction
------------

THAPBI PICT has been designed as a framework which can be applied to multiple
biological contexts, demonstrated in the :ref:`worked examples
<worked_examples>`. Each new marker (i.e. PCR primer target) will require a
new reference database be compiled, most likely starting from published
sequences, but we also sequenced culture collections.

Applied to environmental samples, some primer pairs will amplify a much wider
sequence space than others, either reflecting a more diverse genome region, or
simply a longer sequence. Related to this, the fraction of observed sequences
with a published reference will also vary - and thus the density of the
references in sequence space. This in turn will can change which classifier
algorithm is most appropriate. Inspecting the edit-graph produced for all your
samples and including your initial database entries can help interpret this.

The default classifier allows perfect matches, or a single base pair
difference (substitution, insertion or deletion). This requires good database
coverage with unambiguous sequences, which we have been able to achieve for
the *Phytophthora* ITS1 region targeted by default.

Provided database
-----------------

THAPBI provides a default database which is used when the command line ``-d``
or ``--database`` setting is omitted. This is intended for use with the
*Phytophthora* ITS1 region targeted by the default left and right primers
assumed during preparation of the paired FASTQ files. This is used in the
first :ref:`worked example <woody_hosts>`.

The default database is compiled from the following sets of sequences:

- Curated *Phytophthora* ITS1 sequences (at species level) in a FASTA file,
  imported using the ``thapbi_pict curated-import`` command.
- NCBI *Peronosporales* (including *Phytophthora*) at genus level, using
  the ``thapbi_pict ncbi-import`` command with ``-g`` or ``--genus``.
  i.e. We discard the author provided species information as we found too many
  were misclassified, keeping only the stated genus.
- Observed ITS1 sequences from single isolate positive controls run on a MiSeq
  plate, using the ``thapbi_pict seq-import`` command.
- Four G-BLOCK synthetic controls in a FASTA file, imported using the
  ``thapbi_pict curated-import`` command.

These are vetted against the NCBI taxonomy, which rejects some entries (e.g.
unknown species, or currently listed under unclassified *Phytophthora*).

For further details see the ``database/README.rst`` file in the source code,
and script ``database/build_ITS1_DB.sh`` which automates this.

Ambiguous bases in database
---------------------------

Ideally all the reference sequences in your database will have unambiguous
sequences only (A, C, G and T). However, some published species sequences will
contain IUPAC ambiguity codes, especially if capillary sequenced. How this is
handled will depend on the classifier algorithm used.

For example *Phytophthora condilina* accession ``KJ372262`` has a single ``W``
meaning ``A`` or ``T``. In this case for *P. condilina* in our curated set, we
could select the unambiguous accession ``MG707826`` instead.

With the strictest ``identity`` classifier, the ``W`` will never be matched
(since the Illumina platform does not produce any ambiguous bases other than
``N``). With the default ``onebp`` classifier, this can match but the ``W``
would be the single allowed mismatch (and any database entry with more than
one ambiguity would never be matched). The ``blast`` classifier uses NCBI
BLAST+ internally, and would handle the base as expected.

Conflicting taxonomic assignments
---------------------------------

With any amplicon marker, it is possible that distinct species will share the
exact same sequence. For example, this happens with our ITS1 marker for model
organism *Phytophthora infestans* and sister species *P. andina* and
*P. ipomoeae*. In cases like this where the classifier finds multiple equally
valid taxonomic assignments in the database, they are **all** reported. Should
the user wish however, their database could record a single assignment like
*Phytophthora infestans*-complex.

Our default primers for *Phytophthora* can amplify related genera, not just
*Peronosporales*, but also other *Oomycetes*. Expanding the database coverage
runs into two main problems. First, with less published sequences available,
the default strict classifier fails to match many sequences to a published
sequence. Second, the taxonomic annotation becomes less consistent, with the
same sequence sometimes assigned to multiple genera. This could be due to
misannotation, historical taxonomy not drawing on molecular data, or simply
that this region is more highly conserved in some branches of the *Oomycetes*.

The ``thapbi_pict conflicts`` subcommand can be used to report any conflicts
at species or genus level.
