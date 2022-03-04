Minimum Abundance Threshold
===========================

With less samples multiplexed per sample than our own work (which guided the
default settings), these samples were sequenced at much higher depth.
The m4A plate has from 237,351 to 334,852 raw reads:

.. code:: console

    $ grep -E "(name|m4A)" metadata.tsv | cut -f 1,2,4
    <SEE TABLE BELOW>

As a table:

============= =================== ==========
run_accession library_name        read_count
============= =================== ==========
SRR7109334    m4A-301-1           270478
SRR7109346    m4A-mock3-32000b    290925
SRR7109347    m4A-766-2           324127
SRR7109348    m4A-mock3-32000am4A 296832
SRR7109349    m4A-mock3-16000m4A  284316
SRR7109350    m4A-766-1           302345
SRR7109352    m4A-744-2           258138
SRR7109353    m4A-301-2           296226
SRR7109355    m4A-757-2           277106
SRR7109358    m4A-stds            237351
SRR7109377    m4A-736-1           279896
SRR7109379    m4A-500-1           334852
SRR7109380    m4A-757-1           282080
SRR7109382    m4A-755-2           296389
SRR7109383    m4A-500-2           275990
SRR7109384    m4A-744-1           284153
SRR7109385    m4A-736-2           264055
SRR7109386    m4A-712-2           242161
SRR7109387    m4A-712-1           265985
SRR7109388    m4A-755-1           305617
============= =================== ==========

The second plate we are looking at, m6, has even higher depth:

.. code:: console

    $ grep -E "(name|m6)" metadata.tsv | cut -f 1,2,4
    <SEE TABLE BELOW>

As a table:

============= =============== ==========
run_accession library_name    read_count
============= =============== ==========
SRR7109326    m6-stds         817764
SRR7109327    m6-301-1        890561
SRR7109328    m6-mock3-32000b 943839
SRR7109329    m6-766-1        840068
SRR7109330    m6-744-2        704173
SRR7109331    m6-500-1        911793
SRR7109341    m6-712-2        872265
SRR7109342    m6-500-2        879762
SRR7109343    m6-757-1        903886
SRR7109344    m6-757-2        1210627
SRR7109345    m6-mock3-16000  922440
SRR7109406    m6-712-1        897159
SRR7109408    m6-744-1        778090
SRR7109409    m6-mock3-32000a 1125275
SRR7109411    m6-766-2        785776
SRR7109412    m6-755-1        957067
SRR7109414    m6-736-1        998817
SRR7109415    m6-301-2        1181567
SRR7109417    m6-755-2        1071829
SRR7109418    m6-736-2        919363
SRR7109420    m6-SynMock      1299238
============= =============== ==========

The defaults are an absolute abundance threshold of 100, and a fractional
threshold of 0.1% (i.e. ``-a 100 -f 0.001``).

After merging overlapping reads and primer matching we could expect over 200,000
reads per m4A sample which with the default 0.1% becomes a threshold over 200
reads, and over 600,000 reads per m6 samples giving a read threshold over 600.

So, with this coverage the default fractional abundance threshold of 0.1% makes
the default absolute abundance threshold of 100 redundant. However, this is
still quite cautious, and control samples can help set thresholds objectively.

In this dataset there is a single synthetic control for m6 sequencing run,
library ``SynMock`` aka ``SRR7109420``. We can tell THAPBI PICT at the command
line to use this to set the fractional abundance threshold via ``-y`` or
``--synctrls``, or set the absolute abundance threshold via ``-n`` or
``--negctrls`` (with a list of control file names).
