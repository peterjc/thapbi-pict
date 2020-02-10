Unexpected sequences
====================

In the previous section, we highlighted several unexpected contaminants in the
negative controls which could not be explained as cross-contamination from the
mock community. Likewise the read reports show plenty of unassigned sequences,
things which did not match the very narrow databases built from ``ITS1.fasta``
or ``ITS2.fasta`` containing markers expected from the mock community *only*.

Some unexpected sequences might reflect additional alternative copies of ITS1
in the genomes. Others are likely external contamination - after all there are
fungi practically everywhere.

Amplicon library one - ITS1
---------------------------

Specifically from the first amplicon library for ITS1 we saw the following -
which appear to be unique to the negative controls (not in any mock community
samples) - shown here with their highest single sample abundance, which
supports using a higher minimum abundance threshold:

==================================== ================================== ===
MD5 check sum                        Species                            Max
------------------------------------ ---------------------------------- ---
``daadc4126b5747c43511bd3be0ea2438`` *Wallemia muriae*                   32
``e5b7a8b5dc0da33108cc8a881eb409f5`` *Wallemia muriae*; *Wallemia sebi*  10
``5194a4ae3a27d987892a8fee7b1669b9`` *Trichosporon asahii*               17
``702929cef71042156acb3a28270d8831`` *Candida tropicalis*                14
==================================== ================================== ===

Here are a selection of unclassified reads from the mock communities:

==================================== ======================================== ====
MD5 check sum                        Species                                   Max
------------------------------------ ---------------------------------------- ----
``ee5382b80607f0f052a3ad3c4e87d0ce`` *glomeromycetes*, perhaps *Rhizophagus*   575
``e055cb2efa2e1e0eb32d201e018b8609`` *glomeromycetes*, perhaps *Rhizophagus*    63
``85775735614d45d056ce5f1b67f8d2b2`` *Fusarium*                                109
``5ca0acd7dd9d76fdd32c61c13ca5c881`` *Epicoccum nigrum*; *Epicoccum layuense* 4562
==================================== ======================================== ====

Listed first, ``ee5382b80607f0f052a3ad3c4e87d0ce`` and the less abundant
sequence ``e055cb2efa2e1e0eb32d201e018b8609`` (a 3bp edit away) look like
*glomeromycetes*, perhaps a *Rhizophagus* (from the mock community), but
could be from a *Glomus* species.

Or, ``85775735614d45d056ce5f1b67f8d2b2`` which has good BLAST matches to
several different *Fusarium* species (could be from the mock community).

On the other hand,  ``5ca0acd7dd9d76fdd32c61c13ca5c881`` which perfectly matches
fungus *Epicoccum nigrum* and *Epicoccum layuense*. Present at low levels in
multiple samples, this was the dominant sequence in ``SRR5314339`` aka
``FMockE.HC1_S178``, which was a *high PCR cycle number* replicate of the even
mixture. Perhaps this was a stray fragment of *Epicoccum* which by chance was
amplified early in the PCR? This example was not highlighed in the original
paper but is exactly the kind of thing you should worry about with a high PCR
cycle number.

Amplicon library two - ITS1
---------------------------

These primers amplified a larger fragment to that in the first amplicon library.
Many of the unexpected sequences are therefore the same. For instance, here
``57b06dff740b38bd6a0375abd9db3972`` and 3bp away ``05007e829ab71427b49743994a14105f``
are *glomeromycetes* matches - perhaps *Rhizophagus* (from the mock community).
While ``7e31840276e0b10db32a26cef0abda33`` and 1bp different
``e6f59a48cd979a95dc672dd54c9fcf02`` are probably *Fusarium*. As is
``610caedb1a5699836310fce9dbb9c5fa``.
