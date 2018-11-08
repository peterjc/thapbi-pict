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

import hashlib
import os
import re
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser

from . import __version__
from .db_orm import DataSource, ITS1, SequenceSource, connect_to_db


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


def parse_fasta_entry(text):
    """Split an entry of Clade_SpeciesName_Accession into fields.

    Will also convert the underscores in the species name into spaces:

    >>> parse_fasta_entry('4_P._arenaria_HQ013219')
    ('4', 'P. arenaria', 'HQ013219')

    Note - this assumes function ``split_composite_entry`` has already
    been used to break up any multiple species composite entries.

    Will also cope with the even older legacy format missing the
    underscore between the clade and species:

    >>> parse_fasta_entry('1Phytophthora_aff_infestans_P13660')
    ('1', 'Phytophthora aff infestans', 'P13660')

    And the old variant without any clade at the start, e.g.

    >>> parse_fasta_entry('P._amnicola_CBS131652')
    ('', 'P. amnicola', 'CBS131652')

    And the old variant with just an accession, e.g.

    >>> parse_fasta_entry('VHS17779')
    ('', '', 'VHS17779')

    Dividing the species name into genus, species, strain etc
    is not handled here.
    """
    # The rstrip is to handle this legacy entry:
    # >PPhytophthora_fluvialis_CBS129424_JF701436_ P. fluvialis
    # See also fixing PPhytophthora later on...
    parts = text.rstrip("_").split("_")

    if len(parts) == 1:
        # Legacy variant with just an accession
        return ("", "", text)

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
    acc = parts[-1]

    if name[0] == "P.":
        name[0] = "Phytophthora"

    if clade and not clade_re.fullmatch(clade):
        raise ValueError("Clade %s not recognised from %r" % (clade, text))
    return (clade, " ".join(name), acc)


assert (parse_fasta_entry('4_P._arenaria_HQ013219')
        == ('4', 'Phytophthora arenaria', 'HQ013219'))
assert (parse_fasta_entry('1Phytophthora_aff_infestans_P13660')
        == ('1', 'Phytophthora aff infestans', 'P13660'))
assert (parse_fasta_entry('P._amnicola_CBS131652')
        == ('', 'Phytophthora amnicola', 'CBS131652'))
assert (parse_fasta_entry('ACC-ONLY') == ('', '', 'ACC-ONLY'))


def main(fasta_files, db_url, name=None, debug=True):
    """Run the script with command line arguments."""
    # Connect to the DB,
    Session = connect_to_db(db_url, echo=debug)
    session = Session()

    if not name:
        name = "Legacy import of " + " ".join(
           os.path.basename(_) for _ in fasta_files)

    db_source = DataSource(
        name=name,
        uri=":".join(fasta_files),
        notes="Imported with thapbi_pict legacy_import v%s" % __version__)
    session.add(db_source)

    seq_count = 0
    entry_count = 0
    bad_entry_count = 0
    idn_set = set()
    for filename in fasta_files:
        with open(filename) as handle:
            for title, seq in SimpleFastaParser(handle):
                if title.startswith("Control_"):
                    if debug:
                        sys.stderr.write("Ignoring control entry: %s\n"
                                         % title)
                    continue
                seq_count += 1

                # Here assume the FASTA sequence is already trimmed to the ITS1
                seq_md5 = hashlib.md5(seq.upper().encode("ascii")).hexdigest()
                its1 = ITS1(md5=seq_md5, sequence=seq)
                session.add(its1)

                # One sequence can have multiple entries
                idn = title.split(None, 1)[0]
                if idn in idn_set:
                    sys.stderr.write("WARNING: Duplicated identifier %r\n"
                                     % idn)
                idn_set.add(idn)

                entries = split_composite_entry(title.split(None, 1)[0])
                for entry in entries:
                    entry_count += 1
                    try:
                        clade, name, acc = parse_fasta_entry(entry)
                    except ValueError as e:
                        bad_entry_count += 1
                        sys.stderr.write("WARNING: %s - Can't parse: %r\n"
                                         % (e, idn))
                        continue
                    # Load into the DB
                    # Store "Phytophthora aff infestans" as
                    # genus "Phytophthora", species "aff infestans"
                    genus, species = name.split(None, 1) if name else ("", "")
                    assert genus != "P.", title
                    taxid = 0
                    # Note we use the original FASTA identifier for traceablity
                    # but means the multi-entries get the same source accession
                    record_entry = SequenceSource(source_accession=idn,
                                                  source=db_source,
                                                  its1=its1,
                                                  sequence=seq,
                                                  original_clade=clade,
                                                  original_taxid=taxid,
                                                  original_genus=genus,
                                                  original_species=species,
                                                  current_clade=clade,
                                                  current_taxid=taxid,
                                                  current_genus=genus,
                                                  current_species=species,
                                                  seq_strategy=0,
                                                  seq_platform=0,
                                                  curated_trust=0)
                    session.add(record_entry)
                    # print(clade, species, acc)
    session.commit()
    sys.stderr.write("%i sequences, %i entries including %i bad\n"
                     % (seq_count, entry_count, bad_entry_count))
