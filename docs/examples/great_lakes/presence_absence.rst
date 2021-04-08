Presence and absence
====================

This example includes mock communities which are a controlled setup where we
know what the classifier ought ideally to report for every sample - and all
their expected marker sequences are in the classification database.

The ``thapbi_pict assess`` command run via example script ``run.sh`` uses a
configuration file with all the mock community species for MOL16S, and the
three sphaeriid mussel species for SPH16S - regardless of the target copy
number in the mixture (see Klymus *et al.* (2017) Table 2), or
presence/absence of the fish block.

Of course, just as in the original author's analysis, not everything we think
was present is detected. And *vice versa*, we see some things which are not
classified.

SPH16S
------

This was the more specific primer pair, expected to only amplify sphaeriid
mussel species, so in general we expect less unique sequences than with the
more general MOL16S primers.

Only three members of the mock community should match. Looking at the
``summary/SPH16S.assess.onebp.tsv`` output file in Excel or at the command
line, when run at a minimum abundance threshold of 10, these are the key
numbers:

.. code:: console

    $ cut -f 1-5 summary/SPH16S.assess.onebp.tsv
    <SEE TABLE BELOW>

Or open this in Excel. You should find:

====================== == == == ===
#Species               TP FP FN TN
====================== == == == ===
OVERALL                9  5  0  656
Pisidium compressum    3  0  0  7
Sphaerium corneum      3  0  0  7
Sphaerium nucleus      0  3  0  7
Sphaerium simile       3  1  0  6
Sphaerium striatinum   0  1  0  9
OTHER 62 SPECIES IN DB 0  0  0  620
====================== == == == ===

No false negatives (but we have set the threshold very low), but 5 false
positives: Three cases of *Sphaerium nucleus*, and one each of *S. simile*
and *S. striatinum*.

The *S. nucleus* matches are simply down to an ambiguous sequence in the
database from both this and expected species *S. corneum*. See also the output
from ``thapbi_pict conflicts -d SPH16S.sqlite`` which can report this.

The *S. striatinum* prediction came from ``SPSC3PRO1`` aka ``SRR5534978``, and
is down to several sequences one base pair away the expected *S. simile*
reference, but also one base pair away from an *S. striatinum* database entry.

We already discussed the trace level of 10 reads for *Sphaerium simile* in
mock community sample ``NFSC3PRO3`` using the SOL16S primers. As suggested,
raising the minimum abundance threshold to at least 20 reads would solve this,
but the other false positives here are limitations of the reference set.

MOL16S
------

Looking at the ``summary/MOL16S.assess.onebp.tsv`` output file in Excel or
at the command line, when run at a minimum abundance threshold of 10, these
are the key numbers:

.. code:: console

    $ cut -f 1-5 summary/MOL16S.assess.onebp.tsv
    <SEE TABLE BELOW>

Or open this in Excel. You should find:

========================= == == == ====
#Species                  TP FP FN TN
========================= == == == ====
OVERALL                   74 25 3  1218
Cipangopaludina chinensis 7  0  0  4
Corbicula fluminea        0  1  0  10
Dreissena bugensis        0  8  0  3
Dreissena polymorpha      7  1  0  3
Dreissena rostriformis    7  1  0  3
Gillia altilis            7  0  0  4
Melanoides tuberculata    7  0  0  4
Mytilopsis leucophaeata   7  0  0  4
Pisidium compressum       7  1  0  3
Potamopyrgus antipodarum  7  0  0  4
Sander vitreus            4  0  3  4
Sphaerium corneum         7  1  0  3
Sphaerium nucleus         0  8  0  3
Sphaerium simile          7  3  0  1
Sphaerium striatinum      0  1  0  10
OTHER 105 SPECIES IN DB   0  0  0  1155
========================= == == == ====

This time we do have false negatives - in five of the seven samples the
lowest abundance *Melanoides tuberculata* was not found. Given this was
intended to be present at only 18 copies per sample, this is easily lost with
a minimum abundance threshold designed to screen out noise.

Also, three of the seven samples are missing *Sander vitreus*. Two of these
are from Community 3 where this is intended to be at only 14 copies, the third
was ``SC3PRO2`` aka ``SRR5534972`` for Mock Community 2 MOL16S with Fish Block
Primer, with a target abundance of 72 copies. Here the fish block worked.

Again we have lots of false positives, mostly sister species which reflects
limitations of the reference set.

The exceptions are *Corbicula fluminea* and *Pisidium compressum*, both can be
found in the sample summary ``MOL16S.samples.onebp.xlsx``. The only sample
showing any *Corbicula fluminea* is from ``SC3PRO1`` aka ``SRR5534973``, and
at low abundance. This species was present in the aquaria sample sediment, but
did not amplify there - so cross-contamination seem less likely.

The *Pisidium compressum* came from ``SPSC3PRO2`` aka ``SRR5534981``, which
was meant to be only SPH16S amplicons. As discussed earlier, this seems to be
from primer mixing.

Unknowns
--------

Looking at ``SPH16S.samples.onebp.xlsx`` and ``MOL16S.samples.onebp.xlsx``
even our controls have unknown reads. To study these, next I'd look at the
edit-graphs.
