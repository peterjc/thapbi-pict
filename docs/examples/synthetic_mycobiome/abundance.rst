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
``--negctrls`` (with a list of control file names). It turns out however that
with the default thresholds the control is clean (no unwanted non-synthetic ITS2
reads).

So, there is scope to lower the default thresholds - but how low? We will start
by reproducing the Illumina part of Figure 6, which was based on the m6 MiSeq
sequencing run. This figure explores tag-switching in the demultiplexing, and
in the authors analysis goes as low as 5 reads.

So, the ``run.sh`` example starts by running the pipeline on the m6 dataset with
``-f 0 -a 5`` meaning no fractional threshold and the very low absolute
threshold of 5 reads. This first analysis does *not* use the synthetic control
to raise the threshold on the rest of the samples - we want to see any low level
mixing. We then can compare our sample report against Figure 6.

Look at ``summary/a5.ITS2.samples.onebp.xlsx`` or working at the command line
with the TSV file:


.. code:: console

    $ cut -f 1,3-8 summary/a5.ITS2.samples.onebp.tsv
    <SEE TABLE BELOW>

As a table:

================== ========= ======= ======== ============= ============ ==========
#Sequencing sample Raw FASTQ Flash   Cutadapt Max non-spike Max spike-in Read count
================== ========= ======= ======== ============= ============ ==========
SRR7109326         817764    740627  736334   35300         0            602782
SRR7109327         890561    812674  807956   348111        0            667575
SRR7109328         943839    872263  866253   56120         0            721272
SRR7109329         840068    794475  792260   526536        0            699349
SRR7109330         704173    654661  651528   136471        0            552216
SRR7109331         911793    823392  819469   289230        0            666624
SRR7109341         872265    800475  796363   299243        0            660047
SRR7109342         879762    817277  813470   214155        0            679563
SRR7109343         903886    839105  835431   281533        0            705159
SRR7109344         1210627   1105530 1099959  224635        0            927042
SRR7109345         922440    859262  846519   65686         0            711511
SRR7109406         897159    823034  820146   131937        0            675482
SRR7109408         778090    710762  706659   358092        0            602283
SRR7109409         1125275   1047383 1023234  84748         0            861570
SRR7109411         785776    714894  711176   251097        0            589576
SRR7109412         957067    891942  887650   462496        15           760868
SRR7109414         998817    948348  943426   349965        15           813439
SRR7109415         1181567   1113606 1108129  457441        0            955977
SRR7109417         1071829   987280  982087   589121        0            830400
SRR7109418         919363    858915  854919   282133        0            740813
SRR7109420         1299238   1204532 1199806  187           103014       1028751
================== ========= ======= ======== ============= ============ ==========

Here ``SRR7109420`` is the synthetic control, and it has some non-spike-in
reads present, the most abundant at 187 copies. Conversely, samples
``SRR7109412`` aka ``m6-755-1`` and ``SRR7109414`` aka ``m6-736-1`` have trace
levels of unwanted synthetic spike-in reads, the most abundant at 15 copies.
The counts differ, but these are the same samples highlighted in Figure 6
(sharing the same Illumina i7 index for multiplexing).

As percentages, 187/1199806 gives 0.015% which is nearly ten times lower than
our default of 0.1%. The numbers the other way round are even lower,
15/462496 gives 0.003% and 15/349965 gives 0.004%.
