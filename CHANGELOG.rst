Release History
===============

======= ========== ============================================================================
Version Date       Notes
======= ========== ============================================================================
v0.9.0  *Pending*  Dropped use of Trimmomatic, faster and slightly higher read counts.
v0.8.4  2021-04-13 Speed up re-running by delaying method setup until and if required.
v0.8.3  2021-04-13 Include abundance threshold in summary reports (if varied by sample).
v0.8.2  2021-04-13 Sample report pooling script. Fixed ``-p`` in ``prepare-reads``.
v0.8.1  2021-04-09 Dropped species list embedded in intermediate TSV, ``assess`` needs DB now.
v0.8.0  2021-04-06 Revised genus/species columns in sample report. Added ``scripts/`` folder.
v0.7.11 2021-03-30 ``assess`` now only at sample level. Abundance threshold in ``classify``.
v0.7.10 2021-03-24 Pipeline includes ``fasta-nr`` command making non-redundant FASTA file.
v0.7.9  2021-03-15 Option to show unsequenced entries in summary sample report (``-u``).
v0.7.8  2021-03-11 Only import IUPAC DNA characters to DB. Fixed *N. valdiviana* in default DB.
v0.7.7  2021-02-24 Revised default ITS1 DB: NCBI Oomycetes, more curation & single isolates.
v0.7.6  2021-02-17 Replaced ``seq-import`` with ``curated-seq``, used when building default DB.
v0.7.5  2021-02-16 Refined default DB by adjusting how genus-level NCBI import was trimmed.
v0.7.4  2021-02-15 Edit-graph genus-only labels. New ``1s2g``, ``1s4g`` & ``1s5g`` classifiers.
v0.7.3  2021-01-29 Updated NCBI import & taxonomy. New ``1s3g`` classifier. Use cutadapt v3.0+.
v0.7.2  2020-10-06 Added ``ena-submit`` command for use with interactive ENA read submission.
v0.7.1  2020-09-29 Curated *Phytophthora* DB minor updates. Classifier output in edit-graph.
v0.7.0  2020-04-02 Read counts etc as a header in intermediate FASTA files; shown in reports.
v0.6.15 2020-03-12 Fix regression in read report column sorting.
v0.6.14 2020-03-12 Merge ``read-summary`` & ``sample-summary`` into new ``summary`` command.
v0.6.13 2020-03-09 New classifier method ``substr`` for testing with poorly trimmed DB content.
v0.6.12 2020-03-09 New advanced setting ``--merged-cache`` intended for multiple marker use.
v0.6.11 2020-03-02 Updated genus-level only NCBI import, restrict to those with 32bp leader.
v0.6.10 2020-02-24 Treat I (for inosine as in tRNA) in primers as N (IUPAC code for any base).
v0.6.9  2020-02-20 Allow pre-primer-trimmed FASTQ. Fixed row coloring when missing samples.
v0.6.8  2020-02-17 Metadata ``-x`` default now column 1. Fixed read report metadata captions.
v0.6.7  2020-02-13 Method in ``pipeline`` filenames; max sample abundance in read reports.
v0.6.6  2020-02-05 Coloring groups in ``sample-report``. Can call assessment from ``pipeline``.
v0.6.5  2020-01-27 Do ``--flip`` in ``prepare-reads`` using cutadapt v2.8 or later.
v0.6.4  2020-01-23 ``curated-import`` accepts primers. Reduced memory usage for ``onebp``.
v0.6.3  2020-01-20 Treat NCBI taxonomy "includes" as synonyms, adds 396 new species aliases.
v0.6.2  2020-01-14 Memory optimisation to the default ``onebp`` classifier.
v0.6.1  2020-01-08 Requires at least Python 3.6 as now using f-strings (internal change only).
v0.6.0  2020-01-08 Stop discarding normally conserved 32bp start of *Phytophthora* ITS1 marker.
v0.5.8  2019-12-11 Correction to start of a *P. parsiana* curated sequence in our DB.
v0.5.7  2019-12-09 Replace min bit score with min percentage coverage in ``blast`` classifier.
v0.5.6  2019-12-04 Import species under "unclassified *Phytophthora*" as genus *Phytophthora*.
v0.5.5  2019-12-03 Updated NCBI taxonomy, adds *Phytophthora caryae* and *P. pseudopolonica*.
v0.5.4  2019-12-02 Only use HMM to detect synthetic read negative controls.
v0.5.3  2019-11-25 Replace HMM filter on importing to the database with length check only.
v0.5.2  2019-11-25 Removed redundant use of HMM filter in ``seq-import`` command.
v0.5.1  2019-11-22 Updated NCBI taxonomy, adds *Phytophthora oreophila* and *P. cacuminis*.
v0.5.0  2019-11-21 Only use HMM as a filter, not for trimming in DB import or classify steps.
v0.4.19 2019-11-19 Additional curated entries in default ITS1 database.
v0.4.18 2019-11-19 Reworked ``sample-summary`` table output, now samples vs species with Excel.
v0.4.17 2019-11-15 Control based minimum abundance threshold applied at folders level.
v0.4.16 2019-11-15 Bug fix in ``fasta-nr`` when using input records with descriptions.
v0.4.15 2019-11-04 Harmonised ``dump`` FASTA & ``curated-import`` with semi-colon separator.
v0.4.14 2019-10-23 Configurable FASTA entry separator for ``curated-import`` & ``ncbi-import``.
v0.4.13 2019-10-22 Fix 5 cases missing ``A`` near end, ``...CTGAAAACT`` to ``...CTGAAAAACT``.
v0.4.12 2019-10-22 Removed now unused ``legacy-import`` and ``database/legacy/`` files.
v0.4.11 2019-10-21 Updated the curated DB entries, focused on truncated sequences.
v0.4.10 2019-10-21 New ``curated-import`` command, reworked handling of curated DB entries.
v0.4.9  2019-10-17 New ``sample-summary`` switch ``-q`` / ``--requiremeta``. NetworkX v2.4 fix.
v0.4.8  2019-10-11 New ``fasta-nr`` command for use in alternatives to ``prepare-reads``.
v0.4.7  2019-10-10 New ``--minlen`` & ``--maxlen`` args for ``prepare-reads`` and ``pipeline``.
v0.4.6  2019-10-02 Forgot to include updated DB with the PyPI release.
v0.4.5  2019-10-02 Apply primer trimming to ``ncbi-import`` (crop if primers found).
v0.4.4  2019-10-02 New ``--hmm`` & ``--flip`` arguments for ``prepare-reads`` and ``pipeline``.
v0.4.3  2019-09-26 New ``conflicts`` command reporting genus/species level conflicts in DB.
v0.4.2  2019-09-26 Drop clade from taxonomy table, require unique species entries.
v0.4.1  2019-09-16 Include NCBI strains/variants/etc and their synonyms as species synonyms.
v0.4.0  2019-09-12 NCBI taxonomy synonym support; taxonomy import defaults to all *Oomycetes*.
v0.3.12 2019-09-12 New ``thapbi_pict dump`` option ``-m`` /  ``--minimal`` for DB comparison.
v0.3.11 2019-09-09 Updated default DB and tests to use September 2019 NCBI taxonomy.
v0.3.10 2019-09-05 Handle missing or empty input FASTQ files more gracefully.
v0.3.9  2019-08-14 Log BLAST bit score. Merge ``thapbi assess`` warnings, 3dp for ad-hoc loss.
v0.3.8  2019-08-09 The ``blast`` classifier now applies a minimum BLAST bit score of 100.
v0.3.7  2019-08-05 Added Python API to the main documentation.
v0.3.6  2019-07-19 Added Zenodo FASTQ link to worked example, now includes ``assess`` command.
v0.3.5  2019-07-12 Added missing ``T`` or ``CT`` to 11 of the legacy ITS1 sequences in the DB.
v0.3.4  2019-07-08 Worked example using woody hosts dataset from Riddell *et al.* (2019).
v0.3.3  2019-07-04 Fixed regression in group coloring for ``read-summary`` Excel output.
v0.3.2  2019-07-04 Read The Docs; use ``-i`` / ``--input`` consistently - no positional args.
v0.3.1  2019-06-27 Reformatted documentation to use reStructuredText rather than Markdown.
v0.3.0  2019-06-26 Include four gBlocks synthetic negative controls in DB and pipeline.
v0.2.6  2019-06-25 *Phytophthora* ITS1 HMM threshold now set within model file, not the code.
v0.2.5  2019-06-21 Include XGMML edit-graph (for Cytoscape use) in ``pipeline`` output.
v0.2.4  2019-06-21 Fixed 3 *Hyaloperonospora* also in *Peronospora* in default DB.
v0.2.3  2019-06-18 Sample count rather than total read abundance for node size in edit-graph.
v0.2.2  2019-06-12 New ``edit-graph`` command for use with Cytoscape etc, or PDF via GraphViz.
v0.2.1  2019-05-27 Cope better with multiple (short) ITS1 fragments during classification.
v0.2.0  2019-05-14 Limit ITS1 length, 100 to 250bp. Exclude uncultured NCBI entries from DB.
v0.1.12 2019-05-09 Sort ``read-summary`` output by species. Set coloring group at command line.
v0.1.11 2019-05-06 Excel output from ``read-summary`` with formatting applied.
v0.1.10 2019-05-03 Tweaking command line API, renamed ``plate-summary`` to ``read-summary``.
v0.1.9  2019-05-02 Implemented ``pipeline`` subcommand (prepare reads, classify, and report).
v0.1.8  2019-05-01 Standard errors for missing external tools; versions logged in verbose mode.
v0.1.7  2019-05-01 Changed default classifier method from ``identity`` to more fuzzy ``onebp``.
v0.1.6  2019-04-30 Include ready to use binary ITS1 database in source tar-ball & wheel files.
v0.1.5  2019-04-29 Reworked optional metadata integration and its display in summary reports.
v0.1.4  2019-04-25 Sort samples using the optional metadata fields requested in reports.
v0.1.3  2019-04-24 Can optionally display sample metadata from TSV file in summary reports.
v0.1.2  2019-04-17 Keep searching if ``onebp`` classifier perfect match is at genus-level only.
v0.1.1  2019-04-16 Expand default taxonomy and database from Peronosporaceae to Peronosporales.
v0.1.0  2019-04-04 Include a bundled ITS1 database.
v0.0.15 2019-04-03 Support for genus-level only entries in the database.
v0.0.14 2019-04-01 MD5 in dump output. Fixed importing sequences failing taxonomic validation.
v0.0.13 2019-03-22 Remove conserved 32bp when primer trim. Assess at sample level by default.
v0.0.12 2019-03-11 Fixed bug in ``swarmid`` classifier.
v0.0.11 2019-03-08 Speed up FASTQ preparation by using ``flash`` instead of ``pear`` v0.9.6.
v0.0.10 2019-03-06 Replace primer code allowing only 1bp differences with ``cutadapt``.
v0.0.9  2019-03-05 Looks for expected primers, discards mismatches. Caches HMM files locally.
v0.0.8  2019-02-21 Fix multi-class TN under-counting. New loss metric, ``swarmid`` classifier.
v0.0.7  2019-02-12 Added ``plate-summary`` command, ``onebp`` classifier.
v0.0.6  2019-02-07 Misc. cleanup and import fixes.
v0.0.5  2019-02-06 Hamming Loss in assessment output.
v0.0.4  2019-01-24 Added ``seq-import`` command, ``blast`` classifier, multi-taxon predictions.
v0.0.3  2019-01-22 Simplified generated filenames.
v0.0.2  2019-01-21 Added ``assess`` command.
v0.0.1  2019-01-17 Initial framework with ``identity`` and ``swarm`` classifiers.
======= ========== ============================================================================
