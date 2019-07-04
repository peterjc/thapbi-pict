Release History
===============

======= ========== ============================================================================
Version Date       Notes
======= ========== ============================================================================
v0.3.2  2019-07-04 Read The Docs; using ``-i`` / ``--input`` consistently.
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
