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
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser


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
for _ in ("8b_Phytophthora_brassicae_HQ643158_"
          "8b_P._brassicae_voucher_HQ643156",
          "4_Phytophthora_alticola_HQ013214_4_P._arenaria_HQ013219"):
    assert composite_re.search(_), _


def split_composite_entry(text):
    """Split possibly composite FASTA description into list of strings.

    >>> split_composite_entry('4_Phytophthora_alticola_HQ013214')
    ['4_Phytophthora_alticola_HQ013214']
    >>> split_composite_entry('4_P._alticola_HQ013214_4_P._arenaria_HQ013219')
    ['4_P._alticola_HQ013214', '4_P._arenaria_HQ013219']

    """
    rest = text
    answer = []
    while True:
        split = composite_re.search(rest)
        if not split:
            return answer + [rest]
        index = split.start()
        answer.append(rest[:index])
        # plus one to omit the underscore:
        rest = rest[index + 1:]


assert (split_composite_entry('4_Phytophthora_alticola_HQ013214')
        == ['4_Phytophthora_alticola_HQ013214'])
assert (split_composite_entry('4_Phytophthora_alticola_HQ013214_'
                              '4_P._arenaria_HQ013219')
        == ['4_Phytophthora_alticola_HQ013214', '4_P._arenaria_HQ013219'])
assert (split_composite_entry('4_Phytophthora_alticola_HQ013214_'
                              '4_P._arenaria_HQ013219_4a_P._other_XYZ')
        == ['4_Phytophthora_alticola_HQ013214',
            '4_P._arenaria_HQ013219',
            '4a_P._other_XYZ'])


def main(fasta_files):
    """Run the script with command line arguments."""
    seq_count = 0
    entry_count = 0
    for filename in fasta_files:
        with open(filename) as handle:
            for title, seq in SimpleFastaParser(handle):
                seq_count += 1
                entries = split_composite_entry(title)
                for idn in entries:
                    entry_count += 1
                    print(idn)
                assert len(entries) == (title.count("Phytophthora_") +
                                        title.count("_P._")), title
    print("%i sequences, %i entries" % (seq_count, entry_count))


if __name__ == "__main__":
    main(sys.argv[1:])
