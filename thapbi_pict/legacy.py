# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

"""Code for THAPBI PICT to deal with some legacy FASTA databases.

This is a single place to collect code specific to the legacy ITS1 sequence
FASTA file databases under ``database/legacy/*.fasta`` which hold collections
of ITS1 sequences, with their own naming convention to capture the clade,
species, accession.

This code is needed initially for the specific task of importing those legacy
files into the new database schema.

File ``database/legacy/database.fasta`` mainly contains entries named like
``1Phytophthora_sp_andina_EC3163`` which use the convention of *Clade*,
no space, *Species* (using underscores for spaces), underscore, *Accession*.
There is sometimes a space and then free text. The clade can be ommitted.
Sometimes ``Phytophthora`` is abbreviated to ``P.`` for the genus names.
And there are some further free format entries. Automated parsing of these
is not currently planned.

File ``database/legacy/Phytophthora_ITS_database_v0.004.fasta`` mainly has
entries like ``8a_Phytophthora_sansomeana_EU925375`` which use the convention
of *Clade*, underscore, *Species* (using underscores for spaces), underscore,
*Accession*. Thus in general we can split on the underscores, take the first
entry as the clade, the last as the accession, and those in the middle as the
species name.

It also features 16 composite entries where an ITS1 sequence matches more than
one species, e.g. ``4_Phytophthora_alticola_HQ013214_4_P._arenaria_HQ013219``,
which should be interpretted as ``4_Phytophthora_alticola_HQ013214`` and
``4_Phytophthora_arenaria_HQ013219``. While a different separator would have
been cleaner, looking for underscore, clade, underscore, ``P.``, underscore
will work instead.

File ``database/legacy/Phytophthora_ITS_database_v0.005.fasta`` follows the
same naming pattern but adds four synthetic control sequences named
``Control_05 C1``, ``Control_14 C2``, ``Control_22 C3``, and finally
``Control_42 C4``. Note the use of a space here.
"""

import re

from .db_import import import_fasta_file

# Define a regular expression for the clade naming
clade_re = re.compile(r"\d+[a-z]*")  # note - RE is not start anchored
# Define a regular expression for the composite entry split (using clade)
composite_re = re.compile(r"_\d+[a-z]*_P\._")


# Quick in-line regular expression sanity tests
for _ in ("1", "1a", "1b", "999", "999a", "999z"):
    assert clade_re.fullmatch(_), _  # Needs Python 3.4 or later!
for _ in ("1P", "A"):
    assert not clade_re.fullmatch(_), _
for _ in ("_4_P._", "_7a_P._"):
    assert composite_re.fullmatch(_), _
for _ in ("_P._", "_7a_Phytophthora_"):
    assert not composite_re.fullmatch(_), _
for _ in (
    "8b_Phytophthora_brassicae_HQ643158_" "8b_P._brassicae_voucher_HQ643156",
    "4_Phytophthora_alticola_HQ013214_4_P._arenaria_HQ013219",
):
    assert composite_re.search(_), _


def split_composite_entry(text):
    """Split possibly composite FASTA description into list of strings.

    >>> split_composite_entry('4_Phytophthora_alticola_HQ013214')
    ['4_Phytophthora_alticola_HQ013214']
    >>> split_composite_entry('4_P._alticola_HQ013214_4_P._arenaria_HQ013219')
    ['4_P._alticola_HQ013214', '4_P._arenaria_HQ013219']

    """
    if text.startswith("Control_"):
        # Single entry, want the text after the space preserved.
        return [text]
    rest = text.split(None, 1)[0]  # Only look at the first word
    answer = []
    while True:
        split = composite_re.search(rest)
        if not split:
            return answer + [rest]
        index = split.start()
        answer.append(rest[:index])
        # plus one to omit the underscore:
        rest = rest[index + 1 :]


assert split_composite_entry("4_Phytophthora_alticola_HQ013214") == [
    "4_Phytophthora_alticola_HQ013214"
]
assert split_composite_entry(
    "4_Phytophthora_alticola_HQ013214_" "4_P._arenaria_HQ013219"
) == ["4_Phytophthora_alticola_HQ013214", "4_P._arenaria_HQ013219"]
assert split_composite_entry(
    "4_Phytophthora_alticola_HQ013214_" "4_P._arenaria_HQ013219_4a_P._other_XYZ"
) == ["4_Phytophthora_alticola_HQ013214", "4_P._arenaria_HQ013219", "4a_P._other_XYZ"]


def parse_fasta_entry(text):
    """Split an entry of Clade_SpeciesName_Accession into fields.

    Returns a tuple of taxid (always zero), clade, Species, and an
    empty string (since any extra text in the species is already
    included, unlike the NCBI FASTA parser where the end of the
    species is not well defined). Discards the accession.

    Will also convert the underscores in the species name into spaces:

    >>> parse_fasta_entry('4_P._arenaria_HQ013219')
    ('4', 'P. arenaria', '')

    Note - this assumes function ``split_composite_entry`` has already
    been used to break up any multiple species composite entries.

    Will also cope with the even older legacy format missing the
    underscore between the clade and species:

    >>> parse_fasta_entry('1Phytophthora_aff_infestans_P13660')
    ('1', 'Phytophthora aff infestans', '')

    And the old variant without any clade at the start, e.g.

    >>> parse_fasta_entry('P._amnicola_CBS131652')
    ('', 'P. amnicola', '')

    And the old variant with just an accession, e.g.

    >>> parse_fasta_entry('VHS17779')
    ('', '', '')

    Dividing the species name into genus, species, strain etc
    is not handled here.

    Handles the special case of the synthetic controls as follows:

    >>> parse_fasta_entry("Control_05 C1")
    (0, '', 'synthetic construct C1', '')
    """
    taxid = 0

    # The rstrip is to handle this legacy entry:
    # >PPhytophthora_fluvialis_CBS129424_JF701436_ P. fluvialis
    # See also fixing PPhytophthora later on...
    parts = text.rstrip("_").split("_")

    if parts[0] == "Control":
        # txid32630 is "synthetic construct", seems best taxonomy match
        if " " in text:
            return (
                0,
                "",
                ("synthetic construct %s" % text.split(" ", 1)[1]).rstrip(),
                "",
            )
        else:
            return (
                0,
                "",
                ("synthetic construct %s" % " ".join(parts[1:])).rstrip(),
                "",
            )

    if len(parts) == 1:
        # Legacy variant with just an accession
        return (taxid, "", "", "")

    clade = parts[0]

    if clade in ("PPhytophthora", "P."):
        # Legacy error no clade and either extra P, or just P.
        clade = "Phytophthora"

    if "Phytophthora" in clade:
        # Legacy variant missing the first underscore
        index = clade.index("Phytophthora")
        # Split parts[0] into two list entries
        parts[0] = clade[index:]
        clade = clade[:index]
        parts = [clade] + parts

    name = parts[1:-1]
    # acc = parts[-1]

    if name[0] == "P.":
        name[0] = "Phytophthora"

    # Crop species name at "voucher", e.g. 1_Phytophthora_idaei_voucher_HQ643246
    for crop in ["Voucher", "voucher"]:
        if crop in name:
            name = name[: name.index(crop)]

    # NCBI uses lower case "x" for hybrid species names,
    # while the legacy FASTA files used upper case "X":
    for i in range(len(name)):
        if name[i] == "X":
            name[i] = "x"

    # Handle a couple of synonyms here, where our legacy file and the NCBI
    # do not agree. This is a temporary stop-gap without full synonym support.
    synonyms = {
        "Phytophthora_austrocedri": "Phytophthora austrocedrae",  # taxid 631361
        "Phytophthora_taxon_cyperaceae": "Phytophthora balyanboodja",  # taxid 2054050
        "Phytophthora_alni_subsp._multiformis": "Phytophthora x multiformis",  # 299394
        "Phytophthora_alni_subsp._uniformis": "Phytophthora uniformis",  # taxid 299393
        "Phytophthora_alni_subsp._alni": "Phytophthora x alni",  # taxid 299392
    }
    if "_".join(name) in synonyms:
        name = synonyms["_".join(name)].split()

    if clade and not clade_re.fullmatch(clade):
        raise ValueError("Clade %s not recognised from %r" % (clade, text))
    return (taxid, clade, " ".join(name), "")


assert parse_fasta_entry("Control_05 C1") == (0, "", "synthetic construct C1", "")
assert parse_fasta_entry("7a_Phytophthora_alni_subsp._uniformis_AF139367") == (
    0,
    "7a",
    "Phytophthora uniformis",
    "",
)
assert parse_fasta_entry("8d_Phytophthora_austrocedri_DQ995184") == (
    0,
    "8d",
    "Phytophthora austrocedrae",
    "",
)
assert parse_fasta_entry("6_Phytophthora_taxon_cyperaceae_KJ372258") == (
    0,
    "6",
    "Phytophthora balyanboodja",
    "",
)
assert parse_fasta_entry("4_P._arenaria_HQ013219") == (
    0,
    "4",
    "Phytophthora arenaria",
    "",
)
assert parse_fasta_entry("1Phytophthora_aff_infestans_P13660") == (
    0,
    "1",
    "Phytophthora aff infestans",
    "",
)
assert parse_fasta_entry("P._amnicola_CBS131652") == (
    0,
    "",
    "Phytophthora amnicola",
    "",
)
assert parse_fasta_entry("10_Phytophthora_boehmeriae_Voucher_HQ643149") == (
    0,
    "10",
    "Phytophthora boehmeriae",
    "",
)
assert parse_fasta_entry("ACC-ONLY") == (0, "", "", "")


def main(
    fasta_file, db_url, name=None, validate_species=False, genus_only=False, debug=True
):
    """Implement the thapbi_pict legacy-import command."""
    return import_fasta_file(
        fasta_file,
        db_url,
        name=name,
        debug=debug,
        fasta_entry_fn=split_composite_entry,
        entry_taxonomy_fn=parse_fasta_entry,
        validate_species=validate_species,
        genus_only=genus_only,
    )
