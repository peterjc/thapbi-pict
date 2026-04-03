#!/usr/bin/env python
"""Rename entries in OomyceteDB style FASTA file (via stdin/stdout).

Outputs a FASTA file in the ObiTools style to capture both the
species name and the NCBI taxid.

NOTE: Does not preserve sequence case or any line breaks, gets converted to
upper case as a single line.
"""

import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser

genus_prefices = (
    "Phytophthora_",
    "Phytopythium_",
    "Pythium_",
    "Globisporangium_",
    "Pilasporangia_",
)
suffices = ("_rps10", "_rps10plus", "_rps10_plus", "_rps1", "_(reversed)", "_-")

for title, seq in SimpleFastaParser(sys.stdin):
    idn = title.split(None, 1)[0].strip("_")
    if idn.startswith(genus_prefices) or idn.endswith(suffices):
        # Ad hoc style
        species = idn.replace("'", "").split("|", 1)[0]
        while species.endswith(suffices):
            for _ in suffices:
                if species.endswith(_):
                    species = species[: -len(_)].rstrip("_")
        species = species.replace("_", " ")
        words = species.split()
        if words[1] in ("aff.", "cf.", "sp.", "taxon", "x"):
            species = " ".join(words[:3])
        else:
            species = " ".join(words[:2])
        sys.stdout.write(f">{idn} species_name={species};\n{seq}\n")
        continue
    # Database style?
    species = title.rstrip().rsplit(";")[-1].replace("_", " ").replace("  ", " ")
    species = species.replace(".", ". ").replace("  ", " ").strip()
    try:
        idn = title[title.index("|oodb_id=") :].split("|", 2)[1]
    except ValueError:
        idn = None
    try:
        taxid = title[title.index("|ncbi_taxid=") + 12 :].split("|", 1)[0]
    except ValueError:
        taxid = None
    if idn and taxid:
        sys.stdout.write(f">{idn} species_name={species}; taxid={taxid};\n{seq}\n")
    else:
        sys.stderr.write(f">{title}\n")
