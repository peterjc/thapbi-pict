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


def main(fasta_files):
    """Run the script with command line arguments."""
    print("TODO")
    count = 0
    composite = 0
    for filename in fasta_files:
        with open(filename) as handle:
            for title, seq in SimpleFastaParser(handle):
                count += 1
                idn = title.split(None, 1)[0]  # First word only
                if composite_re.search(idn):
                    composite += 1
                    assert (idn.count("Phytophthora_") +
                            idn.count("_P._") > 1), idn
                else:
                    assert (idn.count("Phytophthora_") +
                            idn.count("_P._") <= 1), idn
                assert clade_re.match(idn)
    print("%i composite" % composite)
    print("%i sequences" % count)


if __name__ == "__main__":
    main(sys.argv[1:])
