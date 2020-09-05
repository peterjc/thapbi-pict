.. _endangered_species_sample_data:

Introduction
============

Data source
-----------

This example is based on the multiple amplicon sequencing libraries from this
paper:

    Arulandhu *et al.* (2017) Development and validation of a multi-locus DNA
    metabarcoding method to identify endangered species in complex samples.
    https://doi.org/10.1093/gigascience/gix080
    https://www.ebi.ac.uk/ena/data/view/PRJEB18620
    https://www.ncbi.nlm.nih.gov/bioproject/PRJEB18620/

There are 177 sequenced samples (6.5GB), 17 experimental mixtures (only two
with replicates, 1.1GB, unusually given as zip files) and 160 inter-laboratory
trials (16 samples with 10 laboratories, 5.4GB).

All the samples were all amplified with a dozen primers (see Table 1), meaning
we need to run THAPBI PICT many times - which is not ideal. Also, to run this
properly you would need a well curated database for each marker.

Provided files
--------------

Either clone the THAPBI PICT source code repository, or decompress the
latest source code release (``.tar.gz`` file). You should find it contains
a directory ``examples/endangered_species/`` which is for this example.

File ``PRJEB18620.txt`` was download from the ENA and includes the raw data
checksums, URLs, but lacks any sample metadata.

Files ``references/*.fasta`` were compiled by hand on an *ad hoc* basis to
use for pre-trimmed reference databases. They *should not* be used as is in
any serious analysis. In many cases ambiguous matches have been omitted in
preference of just species expected in the control mixtures. For example, only
recording *Brassica napus* and *Brassica oleracea* despite some markers being
shared by *Brassica juncea* or *Brassica nigra* etc. In more extreme cases,
markers are clearly not even genus specific, but again only the control
mixture representative appears - e.g. *Carica papaya*, *Glycine max*,
*Gossypium hirsutum*, *Lactuca sativa*, *Solanum lycopersicum*. Deliberately
reducing the false positives from these ambiguous marker sequences was done
for illustrative purposes only.

Note that false positives remain, for example an ITS2 sequence most likely
from *Lactuca sativa* in the control mixture is just one base pair away from
a published sequence from that species (KM210323.1), but perfectly matches
published sequences from *Lactuca altaica*, *L. serriola* and *L. virosa*.

Files ``expected/*.known.tsv`` were compiled by hand from the species content
of the experimental samples (using the PRJEB18620 sample descriptions on the
NCBI and Table 7).

Shell scripts ``setup.py`` and ``run.sh`` should reproduce the analysis
discussed in the THAPBI PICT documentation.


Setup
-----

We assume you have aquired the THAPBI PICT source code, and have your command
line terminal open in the ``examples/endangered_species/`` folder. First we run
the ``setup.py`` script:

.. code:: console

   $ ./setup.py

This first downloads files from the ENA under ``raw_downloads/`` (a mix of
``*.zip`` and ``*.fastq.gz`` files), and then sets up consistently named and
compressed entries under ``raw_data/*.fastq.gz`` instead.

If you have the ``md5sum`` tool installed (standard on Linux), verify the files
downloaded correctly:

.. code:: console

    $ cd raw_download/
    $ md5sum -c MD5SUM.txt
    $ cd ..

There is no need to decompress the files.

Running the pipeline
--------------------

Next, you can run the ``run.py`` script which will call THAPBI PICT multiple
times. Under the ``intermediate/`` folder will be a subdirectory for each of
the primer settings, and the primer name is used as a prefix for the reports
in ``summary/``.

Compared to the other examples, there is an additional ``tmp_merged/``
subfolder which contains gzipped FASTA files after quality trimming and merging
overlapping paired reads into single sequences - but prior to applying the
various primers and abundance thresholds.

Metadata
--------

The sample metadata on the ENA is minimal, although the NCBI SRA has longer
descriptions. For example run ``ERS1545972`` from sample ``SAMEA80893168`` aka
``EM_1`` has title "Experimental mixture 1" but only the NCBI has description
"Experimental mixture containing 99% Bos taurus and 1% Lactuca sativa". Or,
run ``ERS1546502`` from sample ``SAMEA81290668``  aka ``S1_Lab_1`` has title
"Interlaboratory trial" while the NCBI also has the description "Experimental
mixture containing 1% Zea mays, 1% Glycine max, 1% Aloe variegata, 1%
Dendrobium sp., 1% Huso Dauricus, 1% Crocodylus niloticus, 47% Brassica
oleracea and 47% Bos taurus, in dry weight percentages".

File ``PRJEB18620.txt`` with the descriptions on the NCBI SRA, supplemented by
Table 7, was used to write ``metadata.tsv``, which has the following columns:

1. run_accessions, e.g. "ERR1824060;ERR1824061;...;ERR1824075"
2. run_names, e.g. "EM_1" or "S1_Lab_1;S1_Lab_2;...;S1_Lab_16"
3. group, "Experimental mixture" or "Interlaboratory trial"
4. sample, e.g. "EM_1" or "S1"
5. description, e.g. "Experimental mixture containing 1% Zea mays, 1% Glycine
   max, 1% Aloe variegata, 1% Dendrobium sp., 1% Huso Dauricus, 1% Crocodylus
   niloticus, 47% Brassica oleracea and 47% Bos taurus, in dry weight
   percentages."

Note we have a single row for the replicates (two cases in the initial
"Experimental mixture" set, and 16 laboratories for each of the 10
"Interlaboratory trial" samples), cross referenced to the individual runs
with semi-colon separated lists in columns 1 (accession) and 2 (filename).

When calling THAPBI PICT, the meta data commands are given as follows:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -c 3,4,5 -x 2 -g 4

Argument ``-c 3,4,5`` says which columns to display and sort by. This means
group, sample, description. Given the sample prefix naming, putting the group
first is not essential for sorting, but is logical.

Argument ``-x 2`` indicates the filename stem can be found in column 2. Unlike
most of the worked examples, we are not using the accession filenames here.

Argument ``-g 4`` means assign colour bands using sample. This gives 15 thin
bands for the "Experimental mixture" set, and then 10 wide bands for the
"Interlaboratory trial" samples. By chance the two traditional medicine
samples both get wide green bands in the Excel reports.
