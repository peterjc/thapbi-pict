Unexpected sequences
====================

In the previous section, we highlighted several unexpected contaminants in the
negative controls which could not be explained as cross-contamination from the
mock community. Likewise the read reports show plenty of unassigned sequences,
things which did not match the very narrow databases built from ``ITS1.fasta``
or ``ITS2.fasta`` containing markers expected from the mock community *only*.

Some unexpected sequences might reflect additional alternative copies of ITS1
or ITS2 in the genomes. Others are likely external contamination - after all
there are fungi practically everywhere. This seems to have happened on
amplicon library one in the high PCR cycle negative control at least.
Meanwhile, amplicon library two does not have any obvious external
contamination.

Amplicon library one - ITS1 (BITS/B58S3)
----------------------------------------

From the first amplicon library for ITS1 we saw the following sequences in the
negative controls (and by chance, not in any mock community samples) - shown
here with their highest single sample abundance, which supports using a
minimum abundance threshold higher than 10:

==================================== === ==================================
MD5 checksum                         Max Species
------------------------------------ --- ----------------------------------
``daadc4126b5747c43511bd3be0ea2438``  32 *Wallemia muriae*
``e5b7a8b5dc0da33108cc8a881eb409f5``  10 *Wallemia muriae*; *Wallemia sebi*
``5194a4ae3a27d987892a8fee7b1669b9``  17 *Trichosporon asahii*
``702929cef71042156acb3a28270d8831``  14 *Candida tropicalis*
==================================== === ==================================

Here are the reads from entries with a maximum sample abundance over 75
which the ``onebp`` and in some cases ``blast`` based classifier failed to
match, along with the most likely match from reviewing an online NCBI BLAST
search. You can easily extract these entries (and their sequences) from the
bottom of the ``summary/AL1_BITS_B58S3.reads.*.tsv`` files:

==================================== ==== ========================================
MD5 checksum                         Max  Species
------------------------------------ ---- ----------------------------------------
``5ca0acd7dd9d76fdd32c61c13ca5c881`` 4562 *Epicoccum nigrum*; *Epicoccum layuense*
``ee5382b80607f0f052a3ad3c4e87d0ce``  575 *glomeromycetes*, perhaps *Rhizophagus*
``880007c5a18be69c3f444efd144fc450``  236 *Ascochyta* or *Neoascochyta*?
``8e74f38b058222c58943fc6211d277fe``  149 *Fusarium*
``cae29429b90fc6539c440a140494aa25``  114 *glomeromycetes*, perhaps *Rhizophagus*
``85775735614d45d056ce5f1b67f8d2b2``  109 *Fusarium*
==================================== ==== ========================================

The sequence with the top abundance, ``5ca0acd7dd9d76fdd32c61c13ca5c881``,
perfectly matches fungus *Epicoccum nigrum* and *Epicoccum layuense*. Present
at low levels in multiple samples, this was the dominant sequence in
``SRR5314339`` aka ``FMockE.HC1_S178``, which was a *high PCR cycle number*
replicate of the even mixture. Perhaps this was a stray fragment of
*Epicoccum* which by chance was amplified early in the PCR? This example was
not highlighed in the original paper, but is exactly the kind of thing you
should worry about with a high PCR cycle number.

Next ``ee5382b80607f0f052a3ad3c4e87d0ce`` and the less abundant sequence
``cae29429b90fc6539c440a140494aa25`` looks like *glomeromycetes*, perhaps
*Rhizophagus* (from the mock community), but could be from a *Glomus* species.
Using the ``blast`` classifier and the minimal curated reference set matches
this to *Rhizophagus irregularis*, but the situation would be ambiguous in a
more complete database.

Sequence ``880007c5a18be69c3f444efd144fc450`` has perfect matches to lots of
unclassified fungi, and conflicting perfect matches including *Ascochyta* or
*Neoascochyta*. This was seen only in the high PCR cycle number sample
``SRR5314339`` as above.

Next ``8e74f38b058222c58943fc6211d277fe`` and
``85775735614d45d056ce5f1b67f8d2b2`` have good BLAST matches to several
different *Fusarium* species, so could also be from the mock community.

You can find all six of these sequence on the edit-graph, most as isolated grey
nodes along the bottom except ``cae29429b90fc6539c440a140494aa25`` which is 3bp
away from *Rhizophagus irregularis* and linked to it with a dashed line.

So some of the ITS1 sequences in amplicon library one are likely external
contamination - particularly with the high PCR cycle negative control (which
was likely included exactly because of this risk).

Amplicon library two - ITS1 (ITS1f/ITS2)
----------------------------------------

Using our ``blast`` classifier with the 19 species database, everything was
assigned a match. The default ``onebp`` classifier was stricter. For example
while the very common ``f1b689ef7d0db7b0d303e9c9206ee5ad`` (which with the
BITS/B58S3 primers gave ``bb28f2b57f8fddefe6e7b5d01eca8aea``) was matched to
*Fusarium oxysporum*, all the variations of this were too far away from the
database entries for a match.

These primers amplified a larger fragment to that in amplicon library one.
Focusing on those with a sample-abundance over 75 (as in the edit-graphs)
which the ``onebp`` classifier did not match to the curated reference set:

==================================== === =======================================
Long sequence MD5 (ITS1f/ITS2).      Max Species
------------------------------------ --- ---------------------------------------
``57b06dff740b38bd6a0375abd9db3972`` 640 *glomeromycetes*, perhaps *Rhizophagus*
``eed6e5c3881a233cca219f7ffd886bbe`` 315 *glomeromycetes*, perhaps *Rhizophagus*
``05007e829ab71427b49743994a14105f`` 154 *glomeromycetes*, perhaps *Rhizophagus*
``93b2d56429637947243e1b5d54a065cf`` 132 *Fusarium*
``610caedb1a5699836310fce9dbb9c5fa``  96 *Fusarium*
``54aecb27334809f56b7f940b9ca060a3``  93 *Fusarium*
``bd30cf52b7031ddd96e3d7588c1f0e1c``  90 *Fusarium*
``c40cad2530d633430c3805be3740c9a4``  88 *Fusarium*
``d44cd471b11f15e2e42070806737e5d1``  86 *Fusarium*
``831acf596cca4ef840c5543d82e23d16``  82 *Fusarium*
``d4145ba9e3ed6c8c2138ed15b147152d``  81 *Fusarium*
==================================== === =======================================

You can find all of these sequence on the edit-graph, most of those labelled as
likely *Fusarium* are a 1bp edit away from large grey node ``f1b689`` top left
(except ``610caedb1a5699836310fce9dbb9c5fa`` which is an isolated node placed
bottom middle). Those labelled *glomeromycetes* are in the middle near, and in
once case connected to, a dark red *Rhizophagus irregularis* node.

i.e. None of the ITS1 sequences in amplicon library two are clear cut external
contamination.

Amplicon library two - ITS2
---------------------------

Finally, amplicon library two using the ITS3-KYO and ITS4-KYO3 primers for
ITS2. Again, the ``blast`` based classifier matched everything to an entry in
the mock community database. The stricter ``onebp`` classifier assigned most
reads. Here are those few it failed to match with a maximum read abundance
over 75:

==================================== === =============
MD5 checksum                         Max Species
------------------------------------ --- -------------
``d1bb95fff4a7e9958fa3c7f13cc51343`` 211 *Fusarium*
``2ef33e6acd8079d729b81d24b91fcf88`` 133 *Fusarium*
``8edbf2c168b11f910458b0e567ae5fc6``  78 *Aspergillus*
==================================== === =============

These three all appears on the edit-graph separated from a red node (database
entry) by a dashed or dotted line indicating a 2bp or 3bp edit away.

Using an online NCBI BLAST search didn't pin any of these down to species
level, but they do all seem to be fungi. Again, quite a few *Fusarium* matches
which could be alternative ITS2 sequences in the genomes but not in the
curated reference set. Likewise the *Aspergillus* like sequence could be from
the *Aspergillus flavus* in the mock community.

i.e. None of the ITS2 sequences in amplicon library two are clear cut external
contamination.
