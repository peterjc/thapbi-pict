#!/usr/bin/env python3
"""Extract references for unknown FASTA."""
import argparse
import sys
from collections import defaultdict

from Bio.SeqIO.FastaIO import SimpleFastaParser

from thapbi_pict.db_import import parse_ncbi_fasta_entry

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.2")
    sys.exit(0)

# Parse Command Line
usage = """\
The input file should be a FASTA file of observed amplicon sequences (e.g. the
all reads output, or your unknowns), with a reference FASTA file of published
untrimmed sequences expected to contain the amplicon, and optionally another
FASTA file which is a subset of those already trimmed:

$ ./unknowns.py -i thapbi-pict.ITS1.all_reads.identity.tsv \
                -f thapbi-pict.ITS1.all_reads.fasta -o unknowns.fasta

$ ./missed_refs.py -i unknowns.fasta \
                   -f 2022-07-05_ITS1_Oomycota_34111.fasta \
                   -x 2022-07-05_ITS1_Oomycota_w32.fasta \
                   -o 2022-07-05_ITS1_Oomycota_obs.fasta

"""

# TODO - do we want to use sallseqid?

parser = argparse.ArgumentParser(
    prog="unknown_refs.py",
    description="Generate reference FASTA file from BLASTN of unknowns.",
    epilog=usage,
)
parser.add_argument(
    "-i",
    "--input",
    default="/dev/stdin",
    metavar="FILE",
    help="Input amplicons FASTA file. Required, default stdin.",
)
parser.add_argument(
    "-f",
    "--fasta",
    required=True,
    metavar="FILE",
    help="Input untrimmed reference FASTA file. Required.",
)
parser.add_argument(
    "-x",
    "--exclude",
    metavar="FILE",
    help="FASTA file of trimmed references to exclude. Optional.",
)
parser.add_argument(
    "-o",
    "--output",
    dest="output",
    default="/dev/stdout",
    metavar="FASTA",
    help="Output FASTA filename, defaults to stdout.",
)
parser.add_argument(
    "-l",
    "--left",
    default="TTTCCGTAGGTGAACCTGCGGAAGGATCATTA",
    metavar="SEQ",
    help="String to assume at the start of the amplicon, will be added to "
    "partial BLAST matches. Default TTTCCGTAGGTGAACCTGCGGAAGGATCATTA for ITS1.",
)

if len(sys.argv) == 1:
    sys.exit("ERROR: Invalid command line, try -h or --help.")
options = parser.parse_args()


def reject_title(title):
    """Check if a species name be rejected."""
    return title.split(None, 2)[1].lower() in (
        "uncultured",
        "unidentified",
        "environmental",
    )


def species_heuristics(text):
    """Apply NCBI species name heuristics."""
    # Ignore the (zero) taxid this returns
    return parse_ncbi_fasta_entry(text, [])[1]


def generate_references(
    input_fasta, reference_fasta, output_fasta, exclude_fasta, left, sep=";"
):
    """Extract FASTA file of reference sequences."""
    references = defaultdict(set)
    with open(reference_fasta) as handle:
        for title, seq in SimpleFastaParser(handle):
            if reject_title(title):
                continue
            # Take accession.X genus species (hopefully)
            references[seq.upper()].add(
                title.split(".", 1)[0] + ".X " + species_heuristics(title)
            )
    sys.stderr.write(f"Loaded {len(references)} references\n")
    exclude = set()
    if exclude_fasta:
        with open(exclude_fasta) as handle:
            for title, seq in SimpleFastaParser(handle):
                if reject_title(title):
                    continue
                exclude.add(seq.upper())
        sys.stderr.write(f"Will exclude {len(exclude)} trimmed references\n")

    with open(input_fasta) as handle:
        with open(output_fasta, "w") as out_handle:
            for _, seq in SimpleFastaParser(handle):
                seq = seq.upper()
                if seq == (
                    "TGAACCTGCGGAAGGATCATTACCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAA"
                    "CCCTTTTATTGGGGGCTTCTGTCTGGTCTGGCTTCGGCTGGATTGGGTGGCGGCTCTATC"
                    "ATGGCGACCGCTCTGAGCTTCGGCCTGGAGCTAGTAGCCCACTTTTTAAACCCATTCTTA"
                    "ATTACTGAACAAACT"
                ):
                    # Seems to have been truncated...
                    continue
                # sys.stderr.write(f"{title}\n")
                target = seq[len(left) :] if seq.startswith(left) else seq
                # Will try looking for matches which can be extended
                seq_dict = defaultdict(set)
                for ref in references:
                    if target not in ref:
                        continue
                    common = seq
                    while common not in ref and len(left) + len(common) > len(seq):
                        common = common[1:]
                    assert common in ref
                    assert target in common
                    assert len(left) + len(common) >= len(seq)
                    masked = left[: len(seq) - len(common)].lower() + common
                    assert (
                        masked.upper() == seq
                    ), f"{left[:len(seq) - len(common)]} + {common} != {seq} from {ref}"
                    seq_dict[masked].update(references[ref])
                    # sys.stderr.write(f"{title} {references[ref]}\n")

                for seq in sorted(seq_dict, key=lambda s: (s.upper(), s)):
                    out_handle.write(
                        ">%s\n%s\n"
                        % (
                            sep.join(sorted(seq_dict[seq])),
                            seq,
                        )
                    )


generate_references(
    options.input, options.fasta, options.output, options.exclude, options.left
)
