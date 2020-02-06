Unexpected sequences
====================

In the previous section, we highlighted several unexpected contaminants in the
negative controls which could not be explained as cross-contamination from the
mock community.

Specifically from the first amplicon library for ITS1 we saw the following -
which appear to be unique to the negative controls (not in any mock community
samples):

==================================== ==================================
MD5 check sum                        Species
------------------------------------ ----------------------------------
``daadc4126b5747c43511bd3be0ea2438`` *Wallemia muriae*
``e5b7a8b5dc0da33108cc8a881eb409f5`` *Wallemia muriae*; *Wallemia sebi*
``5194a4ae3a27d987892a8fee7b1669b9`` *Trichosporon asahii*
``702929cef71042156acb3a28270d8831`` *Candida tropicalis*
==================================== ==================================

Conversly the read reports show plenty of unassigned sequences - things which
did not match the very narrow database built from ``ITS1.fasta`` or
``ITS2.fasta`` containing the markers expected from the mock community *only*.

For example, ``5ca0acd7dd9d76fdd32c61c13ca5c881`` which perfectly matches
fungus *Epicoccum nigrum* and *Epicoccum layuense*. Present at low levels in
multiple samples, This was the dominant sequence in ``SRR5314339`` aka
``FMockE.HC1_S178``, which was a high PCR cycle number replicate of the even
mixture. Perhaps this was a stray fragment of *Epicoccum* which by chance was
amplified early in the PCR? This example was not highlighed in the original
paper.
