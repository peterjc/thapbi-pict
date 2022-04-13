Minimum Abundance Threshold
===========================

With less samples multiplexed per sample than our own work (which guided the
default settings), these samples were sequenced at much higher depth:

.. code:: console

    $ cut -f 1,2,5 metadata.tsv
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
threshold of 0.1% (i.e. ``-a 100 -f 0.001``). After merging overlapping reads
and primer matching we could expect over 650,000 reads per m6 sample, giving a
threshold over 650 reads.

So, with this coverage the default fractional abundance threshold of 0.1% makes
the default absolute abundance threshold of 100 redundant. However, on this
dataset our defaults are quite cautious, and control samples can help set
thresholds objectively.

In this dataset there is a single synthetic control for m6 sequencing run,
library ``SynMock`` aka ``SRR7109420``. We can tell THAPBI PICT at the command
line to use this to set the fractional abundance threshold via ``-y`` or
``--synctrls``, or set the absolute abundance threshold via ``-n`` or
``--negctrls`` (with a list of control file names). It turns out however that
with the default thresholds the control is clean (no unwanted non-synthetic
ITS2 reads).

So, there is scope to lower the default thresholds - but how low? We will start
by reproducing the Illumina part of Figure 6, which was based on the m6 MiSeq
sequencing run. This figure explores tag-switching in the demultiplexing, and
in the authors' analysis goes as low as 5 reads.

The ``run.sh`` example starts by running the pipeline on the m6 dataset with
``-f 0 -a 2`` to accept everything except singletons (sequences which are only
seen once in a sample; including them gives about ten times as many unique
sequences which slows everything down). This first analysis does *not* use the
synthetic control to raise the threshold on the rest of the samples - we want
to see any low level mixing. We then can compare our sample report against
Figure 6.

Look at ``summary/a2.ITS2.samples.onebp.xlsx`` or working at the command line
with the TSV file:

.. code:: console

    $ cut -f 1,3,5-9,11 summary/a2.ITS2.samples.onebp.tsv
    <SEE TABLE BELOW>

As a table:

============= ================= ========= ======= ======== ============= ============ ==========
#sample_alias Sequencing sample Raw FASTQ Flash   Cutadapt Max non-spike Max spike-in Read count
============= ================= ========= ======= ======== ============= ============ ==========
301-1         SRR7109327        890561    812674  807956   348111        0            687951
301-2         SRR7109415        1181567   1113606 1108129  457441        0            977004
500-1         SRR7109331        911793    823392  819469   289230        0            689176
500-2         SRR7109342        879762    817277  813470   214155        0            699634
712-1         SRR7109406        897159    823034  820146   131937        0            703189
712-2         SRR7109341        872265    800475  796363   299243        0            683058
736-1         SRR7109414        998817    948348  943426   349965        15           834461
736-2         SRR7109418        919363    858915  854919   282133        0            757098
744-1         SRR7109408        778090    710762  706659   358092        0            614989
744-2         SRR7109330        704173    654661  651528   136471        0            564238
755-1         SRR7109412        957067    891942  887650   462496        15           782053
755-2         SRR7109417        1071829   987280  982087   589121        0            848794
757-1         SRR7109343        903886    839105  835431   281533        0            725058
757-2         SRR7109344        1210627   1105530 1099959  224635        0            950458
766-1         SRR7109329        840068    794475  792260   526536        0            712127
766-2         SRR7109411        785776    714894  711176   251097        0            606887
BioMock       SRR7109328        943839    872263  866253   56120         0            744007
BioMock       SRR7109345        922440    859262  846519   65686         0            733784
BioMock       SRR7109409        1125275   1047383 1023234  84748         3            884515
BioMockStds   SRR7109326        817764    740627  736334   35300         0            628576
SynMock       SRR7109420        1299238   1204532 1199806  187           103014       1043533
============= ================= ========= ======= ======== ============= ============ ==========

Here ``SynMock`` (``SRR7109420``) is the synthetic control, and it has some
non-spike-in reads present, the most abundant at 187 copies. Conversely,
samples ``755-1`` (``SRR7109412``), ``736-1`` (``SRR7109414``), and one of the
BioMock samples (``SRR7109409``) have trace levels of unwanted synthetic
spike-in reads, the most abundant at 15, 15 and 3 copies respectively. The
counts differ, but these are all samples highlighted in Figure 6 (sharing the
same Illumina i7 or i5 index for multiplexing). We don't see this in a second
BioMock samples, but our pipeline appears slightly more stringent.

As percentages, 187/1199806 gives 0.015% which is nearly ten times lower than
our default of 0.1%. The numbers the other way round are all even lower,
15/462496 gives 0.003%, 15/349965 gives 0.004%, and 3/1023234 gives 0.003%.

Finally the ``run.sh`` example uses the ``SynMock`` synthetic control to
automatically raise the fractional abundance threshold to 0.015% by including
``-y raw_data/SRR7109420_*.fastq.gz`` in the command line. This brings down
the unique sequence count enough to allow use of a slower but more lenient
classifier as well.

Look at ``summary/ctrl.ITS2.samples.1s5g.xlsx`` or working at the command line
with the TSV file:

.. code:: console

    $ cut -f 1,3,7-10,12 summary/ctrl.ITS2.samples.1s5g.tsv
    <SEE TABLE BELOW>

Note we now get a threshold column showing the absolute threshold applied to
each sample (using the inferred percentage), all above the absolute default of
100, and you can see the total read count has dropped:

============= ================= ======== ========= ============= ============ ==========
#sample_alias Sequencing sample Cutadapt Threshold Max non-spike Max spike-in Read count
============= ================= ======== ========= ============= ============ ==========
301-1         SRR7109327        807956   126       348111        0            579503
301-2         SRR7109415        1108129  173       457441        0            829871
500-1         SRR7109331        819469   128       289230        0            568338
500-2         SRR7109342        813470   127       214155        0            578432
712-1         SRR7109406        820146   128       131937        0            569100
712-2         SRR7109341        796363   125       299243        0            570492
736-1         SRR7109414        943426   148       349965        0            708900
736-2         SRR7109418        854919   134       282133        0            653754
744-1         SRR7109408        706659   111       358092        0            540600
744-2         SRR7109330        651528   102       136471        0            472785
755-1         SRR7109412        887650   139       462496        0            694277
755-2         SRR7109417        982087   154       589121        0            754929
757-1         SRR7109343        835431   131       281533        0            610580
757-2         SRR7109344        1099959  172       224635        0            781213
766-1         SRR7109329        792260   124       526536        0            648525
766-2         SRR7109411        711176   111       251097        0            508838
BioMock       SRR7109328        866253   136       56120         0            607401
BioMock       SRR7109345        846519   132       65686         0            603188
BioMock       SRR7109409        1023234  160       84748         0            718661
BioMockStds   SRR7109326        736334   115       35300         0            526317
SynMock       SRR7109420        1199806  100       187           103014       885058
============= ================= ======== ========= ============= ============ ==========
Note that Palmer *et al.* (2018) apply a threshold to unique sequences, but
the thresholding strategy in THAPBI PICT applies the fractional threshold to
all the samples (given in the same sub-folder as input, so you can separate
your MiSeq runs, or your PCR plates, or just apply a global threshold).

In fact, looking at the read report ``summary/ctrl.ITS2.reads.1s5g.tsv`` it is
clear that while this threshold may have excluded Illumina tag-switching, it
has *not* excluded PCR noise - there are hundreds of low abundance sequences
unique to a single sample. To address that we have to use a considerably
higher threshold, and the default 0.1% is a reasonable choice here.
