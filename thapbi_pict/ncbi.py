r"""Code for THAPBI PICT to deal with importing NCBI FASTA files.

For example, you might perform this search via the NCBI Entrez
website, against the nucleotide database, then send it to file
in FASTA format::

   its1 AND Phytophthora[Organism] AND 150:800[Sequence Length]

Or, using the NCBI Entrez Direct command line tools::

    esearch -db nucleotide \
            -query "its1 AND Phytophthora[Organism] AND 150:800[Sequence Length]" \
            | efetch -format fasta > example.fasta

Then, import this into our ITS DB using::

    thapbi_pict ncbi-import -d example.sqlite example.fasta

"""  # noqa: E501

import sys

from .import_fasta import import_fasta_file


def parse_fasta_entry(text):
    """Split an entry of Accession_Genus_Species_name_Description.

    Note we can't infer the clade without looking up the species,
    so for now this returns an empty clade.

    >>> parse_fasta_entry('LC159493.1 Phytophthora drechsleri genes for ITS1, 5.8S rRNA, ITS2, partial and complete sequence, isolate: PhWa20140918-2')
    ('', 'Phytophthora drechsleri', 'LC159493.1')

    Dividing the species name into genus, species, strain etc
    is not handled here.
    """  # noqa: E501
    parts = text.rstrip().split()
    clade = ''
    acc = parts[0]
    name = parts[1:3]  # assumes "Genus species" only (2 words)
    if len(name[0]) > 2 and name[0].startswith("P."):
        # Special case, but can we assume these are Phytophthora?
        # e.g. Y08654.1 P.cambivora ribosomal internal transcribed spacer, ITS1
        sys.stderr.write(
            "WARNING: Assuming %s from %s is Phytophthora\n" % (name[0], acc))
        name = ["Phytophthora", name[0][2:]]
    if name[0] == "Sequence":
        # Another special case
        # e.g. A57915.1 Sequence 20 from Patent EP0751227
        name = []
    return (clade, " ".join(name), acc)


assert parse_fasta_entry('LC159493.1 Phytophthora drechsleri genes') == \
    ('', 'Phytophthora drechsleri', 'LC159493.1')


def main(fasta_file, db_url, name=None, debug=True):
    """Run the script with command line arguments."""
    return(import_fasta_file(fasta_file, db_url, name=name, debug=debug,
                             # fasta_split_fn=split_composite_entry,
                             fasta_parse_fn=parse_fasta_entry))
