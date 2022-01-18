Release History
===============

======= ========== ============================================================================
Version Date       Notes
======= ========== ============================================================================
v0.11.1 2022-01-18 Using ``rapidfuzz`` rather than ``python-Levenshtein``.
v0.11.0 2022-01-13 Multi-marker reports, pooling predictions from each marker.
v0.10.6 2022-01-12 Fixed slow-down in v0.10.0 on large datasets with small DB.
v0.10.5 2021-12-23 Default for ``-f`` / ``--abundance-fraction`` is now 0.001, meaning 0.1%.
v0.10.4 2021-11-24 Updates to default curated DB, including newer NCBI taxonomy.
v0.10.3 2021-11-19 New ``-f`` / ``--abundance-fraction`` setting, off by default.
v0.10.2 2021-11-05 Updates to default curated DB, and small changes to NCBI taxonomy loading.
v0.10.1 2021-07-28 Fix for using SQLAlchemy v1.3 (previous release needed v1.4).
v0.10.0 2021-07-28 Rework to handle larger DB and multiple markers. Modifies DB schema.
v0.9.9  2021-07-08 Drop SWARM based classifiers. Single intermediate TSV file in pipeline.
v0.9.8  2021-06-17 Drop edit-graph in pipeline. Require full length primers in merged reads.
v0.9.7  2021-06-04 USEARCH SINTAX & OBITools FASTA conventions in ``import`` command.
v0.9.6  2021-05-21 Update default DB taxonomy, Peronosporales & Pythiales up to 450bp only.
v0.9.5  2021-05-10 Simplify to just one ``import`` command for pre-trimmed FASTA input.
v0.9.4  2021-05-05 Drop unused metadata fields from DB schema. Fix GML format edit graphs.
v0.9.3  2021-05-04 Drop HMM for spike-in control detection, now via DB & *k*-mer counting.
v0.9.2  2021-04-28 Fix obscure problem using relative versions of absolute paths.
v0.9.1  2021-04-20 Set metadata encoding. Warn if low abundance threshold. HMM default off.
v0.9.0  2021-04-19 Drop use of Trimmomatic, faster and slightly higher read counts.
v0.8.4  2021-04-13 Sped up re-running by delaying method setup until and if required.
v0.8.3  2021-04-13 Include abundance threshold in summary reports (if varied by sample).
v0.8.2  2021-04-13 Sample report pooling script. Fix ``-p`` in ``prepare-reads``.
v0.8.1  2021-04-09 Drop species list embedded in intermediate TSV, ``assess`` needs DB now.
v0.8.0  2021-04-06 Revise genus/species columns in sample report. Add ``scripts/`` folder.
v0.7.11 2021-03-30 ``assess`` now only at sample level. Abundance threshold in ``classify``.
v0.7.10 2021-03-24 Pipeline includes ``fasta-nr`` command making non-redundant FASTA file.
v0.7.9  2021-03-15 Option to show unsequenced entries in summary sample report (``-u``).
v0.7.8  2021-03-11 Only import IUPAC DNA characters to DB. Fix *N. valdiviana* in default DB.
v0.7.7  2021-02-24 Revise default ITS1 DB: NCBI Oomycetes, more curation & single isolates.
v0.7.6  2021-02-17 ``curated-seq`` replaces ``seq-import``, used when building default DB.
v0.7.5  2021-02-16 Refine default DB by adjusting how genus-level NCBI import trimmed.
v0.7.4  2021-02-15 Edit-graph genus-only labels. New ``1s2g``, ``1s4g`` & ``1s5g`` classifiers.
v0.7.3  2021-01-29 Update NCBI import & taxonomy. New ``1s3g`` classifier. Use cutadapt v3.0+.
v0.7.2  2020-10-06 New ``ena-submit`` command for use with interactive ENA read submission.
v0.7.1  2020-09-29 Curated *Phytophthora* DB minor updates. Classifier output in edit-graph.
v0.7.0  2020-04-02 Read counts etc as a header in intermediate FASTA files; shown in reports.
v0.6.15 2020-03-12 Fix regression in read report column sorting.
v0.6.14 2020-03-12 Merge ``read-summary`` & ``sample-summary`` into new ``summary`` command.
v0.6.13 2020-03-09 New classifier method ``substr`` for testing with poorly trimmed DB content.
v0.6.12 2020-03-09 New advanced setting ``--merged-cache`` intended for multiple marker use.
v0.6.11 2020-03-02 Update genus-level only NCBI import, restrict to those with 32bp leader.
v0.6.10 2020-02-24 Treat I (for inosine as in tRNA) in primers as N (IUPAC code for any base).
v0.6.9  2020-02-20 Allow pre-primer-trimmed FASTQ. Fix row coloring when missing samples.
v0.6.8  2020-02-17 Metadata ``-x`` default now column 1. Fix read report metadata captions.
v0.6.7  2020-02-13 Method in ``pipeline`` filenames; max sample abundance in read reports.
v0.6.6  2020-02-05 Coloring groups in ``sample-report``. Can call assessment from ``pipeline``.
v0.6.5  2020-01-27 Do ``--flip`` in ``prepare-reads`` using cutadapt v2.8 or later.
v0.6.4  2020-01-23 ``curated-import`` accepts primers. Reduce memory usage for ``onebp``.
v0.6.3  2020-01-20 Treat NCBI taxonomy "includes" as synonyms, adds 396 new species aliases.
v0.6.2  2020-01-14 Memory optimisation to the default ``onebp`` classifier.
v0.6.1  2020-01-08 Requires at least Python 3.6 as now using f-strings (internal change only).
v0.6.0  2020-01-08 Stop discarding normally conserved *Phytophthora* ITS1 marker 32bp start.
v0.5.8  2019-12-11 Correction to start of a *P. parsiana* curated sequence in our DB.
v0.5.7  2019-12-09 Replace min bit score with min percentage coverage in ``blast`` classifier.
v0.5.6  2019-12-04 Import species under "unclassified *Phytophthora*" as genus *Phytophthora*.
v0.5.5  2019-12-03 Update NCBI taxonomy, adds *Phytophthora caryae* and *P. pseudopolonica*.
v0.5.4  2019-12-02 Only use HMM to detect synthetic read negative controls.
v0.5.3  2019-11-25 Replace HMM filter on importing to the database with length check only.
v0.5.2  2019-11-25 Remove redundant use of HMM filter in ``seq-import`` command.
v0.5.1  2019-11-22 Update NCBI taxonomy, adds *Phytophthora oreophila* and *P. cacuminis*.
v0.5.0  2019-11-21 Only use HMM as a filter, not for trimming in DB import or classify steps.
v0.4.19 2019-11-19 Additional curated entries in default ITS1 database.
v0.4.18 2019-11-19 Rework ``sample-summary`` table output, now samples vs species with Excel.
v0.4.17 2019-11-15 Control based minimum abundance threshold applied at folders level.
v0.4.16 2019-11-15 Bug fix in ``fasta-nr`` when using input records with descriptions.
v0.4.15 2019-11-04 Harmonise ``dump`` FASTA & ``curated-import`` with semi-colon separator.
v0.4.14 2019-10-23 Configurable FASTA entry separator for ``curated-import`` & ``ncbi-import``.
v0.4.13 2019-10-22 Fix 5 cases missing ``A`` near end, ``...CTGAAAACT`` to ``...CTGAAAAACT``.
v0.4.12 2019-10-22 Remove now unused ``legacy-import`` and ``database/legacy/`` files.
v0.4.11 2019-10-21 Update curated DB entries, focused on truncated sequences.
v0.4.10 2019-10-21 New ``curated-import`` command, rework handling of curated DB entries.
v0.4.9  2019-10-17 New ``sample-summary`` switch ``-q`` / ``--requiremeta``. NetworkX v2.4 fix.
v0.4.8  2019-10-11 New ``fasta-nr`` command for use in alternatives to ``prepare-reads``.
v0.4.7  2019-10-10 New ``--minlen`` & ``--maxlen`` args for ``prepare-reads`` and ``pipeline``.
v0.4.6  2019-10-02 Forgot to include updated DB with the PyPI release.
v0.4.5  2019-10-02 Apply primer trimming to ``ncbi-import`` (crop if primers found).
v0.4.4  2019-10-02 New ``--hmm`` & ``--flip`` arguments for ``prepare-reads`` and ``pipeline``.
v0.4.3  2019-09-26 New ``conflicts`` command reports genus/species level conflicts in DB.
v0.4.2  2019-09-26 Drop clade from taxonomy table, require unique species entries.
v0.4.1  2019-09-16 Include NCBI strains/variants/etc & their synonyms as species synonyms.
v0.4.0  2019-09-12 NCBI taxonomy synonym support; *Oomycetes* default taxonomy import.
v0.3.12 2019-09-12 New ``dump`` option ``-m`` /  ``--minimal`` for DB comparison.
v0.3.11 2019-09-09 Update default DB and tests to use September 2019 NCBI taxonomy.
v0.3.10 2019-09-05 Handle missing or empty input FASTQ files more gracefully.
v0.3.9  2019-08-14 Log BLAST bit score, merge ``assess`` warnings, 3dp for ad-hoc loss.
v0.3.8  2019-08-09 The ``blast`` classifier now applies a minimum BLAST bit score of 100.
v0.3.7  2019-08-05 Add Python API to the main documentation.
v0.3.6  2019-07-19 Add Zenodo FASTQ link to worked example and use ``assess`` command.
v0.3.5  2019-07-12 Add missing ``T`` or ``CT`` to 11 of the legacy ITS1 sequences in the DB.
v0.3.4  2019-07-08 Worked example using woody hosts dataset from Riddell *et al.* (2019).
v0.3.3  2019-07-04 Fix regression in group coloring for ``read-summary`` Excel output.
v0.3.2  2019-07-04 Read The Docs; use ``-i`` / ``--input`` consistently - no positional args.
v0.3.1  2019-06-27 Reformat documentation to use reStructuredText rather than Markdown.
v0.3.0  2019-06-26 Include four gBlocks synthetic negative controls in DB and pipeline.
v0.2.6  2019-06-25 *Phytophthora* ITS1 HMM threshold set within model file, not in code.
v0.2.5  2019-06-21 Include XGMML edit-graph (for Cytoscape use) in ``pipeline`` output.
v0.2.4  2019-06-21 Fix 3 *Hyaloperonospora* also in *Peronospora* in default DB.
v0.2.3  2019-06-18 Sample count rather than total read abundance for edit-graph node size.
v0.2.2  2019-06-12 New ``edit-graph`` command. Use Cytoscape etc, or PDF via GraphViz.
v0.2.1  2019-05-27 Cope better with multiple (short) ITS1 fragments during classification.
v0.2.0  2019-05-14 Limit ITS1 length, 100-250bp. Exclude uncultured NCBI entries from DB.
v0.1.12 2019-05-09 Sort ``read-summary`` by species. Set coloring group at command line.
v0.1.11 2019-05-06 Excel output from ``read-summary`` with formatting applied.
v0.1.10 2019-05-03 Tweak command line API, renamed ``plate-summary`` to ``read-summary``.
v0.1.9  2019-05-02 New ``pipeline`` subcommand (prepare reads, classify, and report).
v0.1.8  2019-05-01 Standard errors for missing external tools. Log versions in verbose mode.
v0.1.7  2019-05-01 Chang default classifier method from ``identity`` to more fuzzy ``onebp``.
v0.1.6  2019-04-30 Include ready to use binary ITS1 DB in source tar-ball & wheel files.
v0.1.5  2019-04-29 Rework optional metadata integration and its display in summary reports.
v0.1.4  2019-04-25 Sort samples using the optional metadata fields requested in reports.
v0.1.3  2019-04-24 Can optionally display sample metadata from TSV file in summary reports.
v0.1.2  2019-04-17 Keep searching if ``onebp`` classifier perfect match is at genus-level only.
v0.1.1  2019-04-16 Expand default taxonomy & DB from Peronosporaceae to Peronosporales.
v0.1.0  2019-04-04 Include a bundled ITS1 DB.
v0.0.15 2019-04-03 Support for genus-level only entries in the DB.
v0.0.14 2019-04-01 MD5 in dump output. Fix importing sequences failing taxonomic validation.
v0.0.13 2019-03-22 Drop conserved 32bp when primer trim. Assess at sample level by default.
v0.0.12 2019-03-11 Fix bug in ``swarmid`` classifier.
v0.0.11 2019-03-08 Sped up FASTQ preparation by using ``flash`` instead of ``pear`` v0.9.6.
v0.0.10 2019-03-06 Replace primer code allowing only 1bp differences with ``cutadapt``.
v0.0.9  2019-03-05 Look for expected primers, discards mismatches. Cache HMM files locally.
v0.0.8  2019-02-21 Fix multi-class TN under-counting. New loss metric, ``swarmid`` classifier.
v0.0.7  2019-02-12 New ``plate-summary`` command, ``onebp`` classifier.
v0.0.6  2019-02-07 Misc. cleanup and import fixes.
v0.0.5  2019-02-06 Hamming Loss in assessment output.
v0.0.4  2019-01-24 New ``seq-import`` command, ``blast`` classifier, multi-taxon predictions.
v0.0.3  2019-01-22 Simplify generated filenames.
v0.0.2  2019-01-21 New ``assess`` command.
v0.0.1  2019-01-17 Initial framework with ``identity`` and ``swarm`` classifiers.
======= ========== ============================================================================
