"""Summarise ITS1 classification results at folder (plate) level.

This implements the ``thapbi_pict plate-summary ...`` command.
"""

import os
import sys

from collections import Counter

from Bio.SeqIO.FastaIO import SimpleFastaParser

from .utils import find_requested_files
from .utils import load_metadata
from .utils import parse_species_tsv
from .utils import sample_sort
from .utils import split_read_name_abundance


def color_bands(meta, meta_names, sample_color_bands, debug=False):
    """Return a list for formats, one for each sample."""
    min_groups = 5
    max_groups = 0.8 * len(meta)
    if min_groups > max_groups:
        if debug:
            sys.stderr.write("DEBUG: Not enough colums to bother with coloring\n")
    else:
        for first_word in (False, True):
            for i, name in enumerate(meta_names):
                values = [_[i] for _ in meta]
                # if debug:
                #    sys.stderr.write("DEBUG: %s %r\n" % (name, values))
                if len(set(values)) == 1:
                    # All the same not helpful for banding
                    continue
                if first_word:
                    # Try massaging the data, take first 'word' only
                    values = [
                        _.replace("_", " ").replace("-", " ").split(None, 1)[0]
                        if _
                        else ""
                        for _ in values
                    ]
                if len(set(values)) < min_groups or max_groups < len(set(values)):
                    # (Almost) all same or (almost) unique not helpful
                    continue
                bands = []
                for s, value in enumerate(values):
                    if s == 0:
                        bands.append(0)
                    elif value == values[s - 1]:
                        # Same
                        bands.append(bands[-1])
                    else:
                        # Different
                        if value in values[:s]:
                            # Metadata values are not grouped, can't use for banding
                            # (might be able to use with a color key?)
                            bands = None
                            if debug:
                                sys.stderr.write(
                                    "DEBUG: Could not use %s for coloring, "
                                    "not grouped nicely\n" % name
                                )
                            break
                        bands.append(max(bands) + 1)
                if bands:
                    assert len(set(bands)) == max(bands) + 1
                    if debug:
                        sys.stderr.write(
                            "DEBUG: Using %s%s for coloring columns (%i groups)\n"
                            % (
                                name,
                                " (first word)" if first_word else "",
                                len(set(bands)),
                            )
                        )
                    return [
                        sample_color_bands[_ % len(sample_color_bands)] for _ in bands
                    ]
    if debug:
        sys.stderr.write(
            "DEBUG: No metadata deemed suitable to use for coloring sample columns\n"
        )
    return [None] * len(meta)


def main(
    inputs,
    output,
    method,
    min_abundance=1,
    excel=None,
    metadata_file=None,
    metadata_cols=None,
    metadata_fieldnames=None,
    metadata_index=None,
    debug=False,
):
    """Implement the thapbi_pict plate-summary command.

    The expectation is that the inputs represent all the samples
    from one (96 well) plate, or some other meaningful batch.
    """
    assert isinstance(inputs, list)

    if not output:
        sys.exit("ERROR: No output file specified.\n")

    samples = set()
    md5_abundance = Counter()
    abundance_by_samples = {}
    md5_species = {}
    md5_to_seq = {}

    if debug:
        sys.stderr.write("Loading FASTA sequences and abundances\n")
    for fasta_file in find_requested_files(
        [_ for _ in inputs if not _.endswith(".tsv")], ".fasta", debug
    ):
        sample = os.path.basename(fasta_file).rsplit(".", 1)[0]
        samples.add(sample)
        with open(fasta_file) as handle:
            for title, seq in SimpleFastaParser(handle):
                md5, abundance = split_read_name_abundance(title.split(None, 1)[0])
                if min_abundance > 1 and abundance < min_abundance:
                    continue
                abundance_by_samples[md5, sample] = abundance
                md5_abundance[md5] += abundance
                md5_to_seq[md5] = seq
                md5_species[md5] = set()
    samples = sample_sort(samples)
    (
        metadata_rows,
        metadata_samples,
        meta_names,
        meta_default,
        missing_meta,
    ) = load_metadata(
        metadata_file,
        metadata_cols,
        metadata_fieldnames,
        metadata_index,
        sequenced_samples=samples,
        metadata_sort=True,
        debug=debug,
    )
    # Turn row-centric metadata into a dictionary keyed on sequenced sample name
    metadata = {}
    new = []
    # Already sorted rows of the metadata values, discarded the order in the table
    for row, r_samples in zip(metadata_rows, metadata_samples):
        for sample in r_samples:
            if sample in samples:
                # print(sample, row)
                metadata[sample] = row
                assert sample not in new, sample
                new.append(sample)
    for sample in missing_meta:
        # print("Missing metadata for %s" % sample)
        assert sample not in new, sample
        new.append(sample)
    assert set(samples) == set(new)
    assert len(samples) == len(new)
    samples = new
    del metadata_rows, metadata_samples, new

    methods = method.split(",")
    for method in methods:
        if debug:
            sys.stderr.write("Loading predictions for %s\n" % method)
        tsv_files = find_requested_files(
            [_ for _ in inputs if not _.endswith(".fasta")], ".%s.tsv" % method, debug
        )
        if len(samples) != len(tsv_files):
            sys.exit(
                "ERROR: Identified %i samples from FASTA files, but %i TSV files for %s"
                % (len(samples), len(tsv_files), method)
            )
        for predicted_file in tsv_files:
            sample = os.path.basename(predicted_file).rsplit(".", 2)[0]
            assert sample in samples, predicted_file
            # TODO: Look at taxid here?
            for name, _, sp in parse_species_tsv(predicted_file, min_abundance):
                md5, abundance = split_read_name_abundance(name)
                if min_abundance > 1 and abundance < min_abundance:
                    continue
                assert abundance_by_samples[md5, sample] == abundance, name
                # Combining over all methods!
                if sp:
                    md5_species[md5].update(sp.split(";"))

    if output == "-":
        if debug:
            sys.stderr.write("DEBUG: Output to stdout...\n")
        handle = sys.stdout
    else:
        handle = open(output, "w")

    if excel:
        import xlsxwriter

        workbook = xlsxwriter.Workbook(excel)
        worksheet = workbook.add_worksheet("Sequence vs samples")
        cell_rightalign_format = workbook.add_format({"align": "right"})
        red_conditional_format = workbook.add_format(
            {"bg_color": "#FFC7CE", "font_color": "#000000"}
        )
        sample_color_bands = [
            # Simple rolling rainbow pastel pallet
            workbook.add_format({"bg_color": c.upper(), "font_color": "#000000"})
            for c in [
                "#ffb3ba",  # pink
                "#ffdfba",  # orange
                "#ffffba",  # yellow
                "#baffc9",  # geen
                "#bae1ff",  # blue
            ]
        ]

        # If there are lots of samples, set narrow column widths
        if len(samples) > 50:
            # Set column width to 2
            worksheet.set_column(5, 5 + len(samples), 2)
        elif len(samples) > 20:
            # Set column width to 4
            worksheet.set_column(5, 5 + len(samples), 4)
    else:
        workbook = None
        worksheet = None

    current_row = 0
    first_data_row = 0
    sample_formats = [None] * len(samples)
    if metadata:
        # Insert extra header rows at start for sample meta-data
        # Make a single metadata call for each sample
        meta = [metadata.get(sample, meta_default) for sample in samples]
        for i, name in enumerate(meta_names):
            handle.write("#\t\t\t\t%s\t%s\n" % (name, "\t".join(_[i] for _ in meta)))
        if worksheet:
            sample_formats = color_bands(
                meta, meta_names, sample_color_bands, debug=debug
            )
            for i, name in enumerate(meta_names):
                worksheet.write_string(i, 4, name, cell_rightalign_format)
                for s, sample in enumerate(samples):
                    worksheet.write_string(
                        i,
                        5 + s,
                        metadata.get(sample, meta_default)[i],
                        sample_formats[s],
                    )
            current_row += len(meta_names)
    handle.write(
        "#ITS1-MD5\t%s-predictions\tSequence\tSample-count\tTotal-abundance\t%s\n"
        % (",".join(methods), "\t".join(samples))
    )
    handle.write(
        "TOTAL\t-\t-\t%i\t%i\t%s\n"
        % (
            sum(
                1
                for md5 in md5_to_seq
                for sample in samples
                if (md5, sample) in abundance_by_samples
            ),
            sum(md5_abundance.values()),
            "\t".join(
                str(
                    sum(
                        abundance_by_samples.get((md5, sample), 0) for md5 in md5_to_seq
                    )
                )
                for sample in samples
            ),
        )
    )
    if worksheet:
        worksheet.write_string(current_row, 0, "ITS1-MD5")
        worksheet.write_string(current_row, 1, "%s-predictions" % ",".join(methods))
        worksheet.write_string(current_row, 2, "Sequence")
        worksheet.write_string(current_row, 3, "Sample-count")
        worksheet.write_string(current_row, 4, "Total-abundance")
        for s, sample in enumerate(samples):
            worksheet.write_string(current_row, 5 + s, sample, sample_formats[s])
        current_row += 1
        first_data_row = current_row
        worksheet.write_string(current_row, 0, "TOTAL")
        worksheet.write_string(current_row, 1, "-")
        worksheet.write_string(current_row, 2, "-")
        worksheet.write_number(
            current_row,
            3,
            sum(
                1
                for md5 in md5_to_seq
                for sample in samples
                if (md5, sample) in abundance_by_samples
            ),
        )
        worksheet.write_number(current_row, 4, sum(md5_abundance.values()))
        for s, sample in enumerate(samples):
            worksheet.write_number(
                current_row,
                5 + s,
                sum(abundance_by_samples.get((md5, sample), 0) for md5 in md5_to_seq),
                sample_formats[s],
            )
        current_row += 1
        worksheet.freeze_panes(current_row, 5)  # keep total line in view plus headers
    for total_abundance, md5 in reversed(
        sorted((v, k) for (k, v) in md5_abundance.items())
    ):
        md5_in_xxx_samples = sum(1 for _ in samples if (md5, _) in abundance_by_samples)
        handle.write(
            "%s\t%s\t%s\t%i\t%i\t%s\n"
            % (
                md5,
                ";".join(sorted(md5_species[md5])),
                md5_to_seq[md5],
                md5_in_xxx_samples,
                total_abundance,
                "\t".join(str(abundance_by_samples.get((md5, _), 0)) for _ in samples),
            )
        )
        if worksheet:
            worksheet.write_string(current_row, 0, md5)
            worksheet.write_string(current_row, 1, ";".join(sorted(md5_species[md5])))
            worksheet.write_string(current_row, 2, md5_to_seq[md5])
            worksheet.write_number(current_row, 3, md5_in_xxx_samples)
            worksheet.write_number(current_row, 4, total_abundance)
            for s, sample in enumerate(samples):
                worksheet.write_number(
                    current_row,
                    5 + s,
                    abundance_by_samples.get((md5, sample), 0),
                    sample_formats[s],
                )
            current_row += 1
            worksheet.conditional_format(
                first_data_row,
                5,
                current_row,
                5 + len(samples),
                {
                    "type": "cell",
                    "criteria": "greater than",
                    "value": 0,
                    "format": red_conditional_format,
                },
            )

    if output != "-":
        handle.close()
    if workbook:
        workbook.close()

    try:
        sys.stdout.flush()
    except BrokenPipeError:
        pass
    try:
        sys.stderr.flush()
    except BrokenPipeError:
        pass

    return 0
