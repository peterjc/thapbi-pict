Minimum Abundance Threshold
===========================

Having looked at the positive and negative controls to help set the minimum
abundance threshold at 50, how do the THAPBI PICT results compare to the
authors' own results for the ITS1 marker? See Landa *et al.* (2021) Table 1,
column 5, reproduced in file ``metadata.tsv`` provided with this example.

Samples without any Phytophthora
--------------------------------

There are multiple "-" entries Landa *et al.* (2021) Table 1, covering both
samples which were not Illumina sequenced, and samples which were sequenced
but not found to contain any Phytophthora.

The example (script ``run.sh``) runs the THAPBI PICT pipeline without using
the ``-u`` or ``--unsequenced`` option which would include the unsequenced
samples in the reports. We will initially focus on those without any
Phytophthora.

Open ``summary/british_soil_its1.samples.onebp.xlsx`` in Excel, or looking at
the plain text TSV equivalent, you can sort or filter for a "-" in column 5,
the author's ITS result. Here all the THAPBI PICT output *Phytophthora*
columns are zero.

.. code:: console

    $ cut -f 2,5,8,9,13 summary/british_soil_its1.samples.onebp.tsv | grep -E "(ITS|\t-\t)"
    <SEE TABLE BELOW>

As a table:

======== === ================= ====================================== ==========
Location ITS Sequencing sample Classification summary                 Read count
======== === ================= ====================================== ==========
CN1      -   SRR13393836       -                                      0
CN2      -   SRR13393824       -                                      0
DS1-9    -   SRR13393831       Unknown                                469
DS2-6    -   SRR13393827       Peronospora (unknown species), Unknown 4415
DS3-9    -   SRR13393818       Unknown                                3447
DS4-1    -   SRR13393815       Unknown                                4266
DS4-8    -   SRR13393810       Unknown                                1011
DS4-9    -   SRR13393809       Unknown                                321
DS5-9    -   SRR13393804       Unknown                                137
DS6-1    -   SRR13393803       Unknown                                3134
DS6-10   -   SRR13393798       Bremia (unknown species), Unknown      5076
DS6-2    -   SRR13393801       Unknown                                2209
DS6-3    -   SRR13393800       Unknown                                1319
DS9-3    -   SRR13393786       Unknown                                2820
US2-10   -   SRR13393757       Unknown                                3024
US2-8    -   SRR13393761       Unknown                                2624
US2-9    -   SRR13393759       Unknown                                1439
US5-2    -   SRR13393738       Peronospora (unknown species), Unknown 3676
======== === ================= ====================================== ==========

Column "Read count" gives the number of reads from the ``prepare-reads`` step,
after applying the abundance threshold. This is zero for the two negative
controls and low for several of the disturbed soil sites, but in the thousands
for the most of the locations.

Where Landa *et al.* (2021) reported no detectable *Phytophthora*, THAPBI PICT
agrees. In one case it has identified *Bremia*, and twice *Peronospora*, but
in general we have lots of unknown (unclassified) sequences.

However, in four files we unexpectedly found only unknowns (i.e. no detectable
*Phytophthora*), while Landa *et al.* (2021) reported *Phytophthora*:

.. code:: console

    $ cut -f 2,5,8,9,13 summary/british_soil_its1.samples.onebp.tsv | grep -E "(ITS|\tUnknown\t)" | grep -v "\t-\t"
    ...

With minor abbreviations to fit the table nicely:

======== ================ ================= ====================== ==========
Location ITS              Sequencing sample Classification summary Read count
======== ================ ================= ====================== ==========
DS1-5    P. austrocedri   SRR13393833       Unknown                2917
DS5-3    P. siskiyouensis SRR13393807       Unknown                2771
DS9-5    P. fallax, ...   SRR13393784       Unknown                2376
DS9-7    P. austrocedri   SRR13393783       Unknown                946
======== ================ ================= ====================== ==========

Looking at ``SRR13393833`` from ``DS1-5``, THAPBI PICT does not report any
*P. austrocedri* at this abundance threshold, but can find traces with just a
few reads supporting matching reads. An NCBI BLAST search of the dominant
sequences in this sample suggest *Plasmopara*, but that is inconclusive with
the current NCBI NT database. Having checked with the authors, they too found
only very low levels of *P. austrocedri* in this sample.

Looking at ``SRR13393807`` from ``DS5-3``, THAPBI PICT does not report any
*P. siskiyouensis* at this abundance threshold - but finds a trace with a
sequence found at 20 copies. However, the dominant sequence in this sample
gives a good NCBI BLAST match to *P. siskiyouensis* at 99% identical, 170/172
matches. That sequence is in our database, but the match is just outside the
1bp threshold of our default ``onebp`` classifier. A more relaxed classifier
method could be used instead, perhaps ``1s3g`` which would assign
*Phytophthora* (unknown species) here.

The last two are easily explained, DS9-5 and DS9-7 were Illumina sequenced
twice, and here we have DS-9-5a (``SRR13393784``) and DS-9-7a
(``SRR13393783``) which were apparently failures and resequenced as DS-9-5b
(``SRR13393777``) and DS-9-7b (``SRR13393775``).

Phytophthora Results
--------------------

Direct comparison of results is complicated by several factors. Landa *et al.*
(2021) assigned labels to candidate species as follows:

    The other 21 ASVs were named as *Phytophthora* sp. uncultured 1a to 21a.
    Their sequences showed high similarity with sequences of already
    known/undescribed *Phytophthora* species such as *P. alni/uniformis*,
    *P. cambivora*, *P. capsici/glovera*, *P. europaea/megasperma*,
    *P. iranica/clandestina*, *P. melonis/sinensis*, *P. quercina/P. sp.
    ohioensis*, *P. sojae* and *P. uliginosa*, but their ITS sequence homology
    was below 99%, ...

In contrast, in most of these cases we have found a match only one base pair
away in our curated database - and thus many more named species are reported.

There are also differences in granularity. For example, they have clumped
*P. plurivora/citricola* which we have not, while conversely THAPBI PICT has
for example clumped *P. capsici* with *P. glovera* as indistinguishable from
the ITS1 marker.

However, there is one interesting difference to highlight. Quoting from their
results text:

    Eight *Phytophthora* species were detected in the soil samples when using
    both ITS and COI regions including *P. cactorum*, *P. cinnamomi*,
    *P. megasperma*, *P. plurivora/citricola*, *P. primulae*,
    *P. pseudosyringae*, *P. ramorum* and *P. syringae*. On the other hand,
    12 species were only detected when using the ITS region including
    *P. austrocedri*, *P. capsici*, *P. castaneae*, *P. fallax*,
    *P. foliorum*, *P. idaei*, *P. kernoviae*, *P. lacustris*, *P. obscura*,
    *P. rubi/fragariae*, *P. siskiyouensis* and *P. uniformis*, whereas four
    species were detected only when using the COI region including
    *P. europaea*, *P. gonapodyides*, *P. quercina*, and *P. uliginosa*

That last point is interesting as we do find sequences matching those four
species in the ITS1 data. You might prefer to look at the read or sample
reports in Excel to find the rows and columns highlighted here at the command
line:

.. code:: console

    $ grep -E "(^#Marker|Phytophthora europaea|\tLocation)" \
      summary/british_soil_its1.reads.onebp.tsv | cut -f 1,4,22,30
    #                                               DS2-2        DS3-4
    #Marker-MD5                       Sample-count  SRR13393830  SRR13393820
    3d3fa2fd6fe0f183cad80771f5950b27  2             3501         630


Sample ``DS2-2`` (where the authors found *P. europaea* in the COI data) and
``DS3-4`` have the same sequence ``3d3fa2fd6fe0f183cad80771f5950b27`` which
is an equally good match to three species: *P. europaea*, *P. flexuosa* and
*P. tyrrhenica*.

For *P. gonapodyides* there were several different sequence variants, but they
all matched this species only. They came from four samples including ``DS9-6``
and ``DS1-8`` (but not ``US4-6``) where the authors could isolate it:

.. code:: console

    $ grep -E "(^#Marker|Phytophthora gonapodyides|\tLocation)" \
      summary/british_soil_its1.reads.onebp.tsv | cut -f 1,4,50,68,78,85
    #                                               DS7-1        DS9-6        US1-8        US2-5
    #Marker-MD5                       Sample-count  SRR13393797  SRR13393776  SRR13393770  SRR13393762
    ed15fefb7a3655147115fc28a8d6d671  3             116          998          0            1054
    c1a720b2005f101a9858107545726123  1             0            0            325          0
    96e0e2f0475bd1617a4b05e778bb04c9  1             0            0            178          0
    b7ca9a8e6388b39fa2d886e19b8f67ba  1             0            149          0            0
    69ecbd0dba57c8bf258f21109bd81917  1             0            0            0            66
    2cc48f88c9174b8d3c38f4546c2d402e  1             0            62           0            0


Only a single unique sequence matched *P. quercina* but it also matches
*P. castanetorum*, found in three samples, ``DS3-3``, ``DS7-1`` and ``DS9-9``
(including ``DS7-1`` where the authors found *P. quercina* in the COI data):

.. code:: console

    $ grep -E "(^#Marker|Phytophthora quercina|\tLocation)" \
      summary/british_soil_its1.reads.onebp.tsv | cut -f 1,4,29,50,71
    #                                               DS3-3        DS7-1        DS9-9
    #Marker-MD5                       Sample-count  SRR13393821  SRR13393797  SRR13393774
    ec642d4b8148085bb3f426829665755d  3             4230         2178         1080


And then two unique sequences matching just *P. uliginosa* from three samples,
``DS2-4`` (where the authors detected it with COI), ``US2-1`` and ``US3-2``:

.. code:: console

    $ grep -E "(^#Marker|Phytophthora uliginosa|\tLocation)" \
      summary/british_soil_its1.reads.onebp.tsv | cut -f 1,4,24,80,90
    #                                               DS2-4        US2-1        US3-2
    #Marker-MD5                       Sample-count  SRR13393828  SRR13393765  SRR13393755
    cab2875a481e358a4993910872ec53a6  2             953          0            195
    b403844130bc3d20fcc5356a31c19bca  1             0            5476         0

So with our pipeline we could find those species in the ITS1 marker data,
although not always unambiguously as sadly this marker is not always species
specific for *Phytophthora*.
