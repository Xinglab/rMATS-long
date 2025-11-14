import argparse

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description="Cleanup gene IDs from ESPRESSO's updated.gtf")
    parser.add_argument('--in-gtf',
                        required=True,
                        help='the updated.gtf from ESPRESSO')
    parser.add_argument('--out-gtf',
                        required=True,
                        help='where to write the cleaned gtf')

    return parser.parse_args()


def parse_gtf_attributes(attributes_str):
    attributes = dict()
    attribute_pairs = attributes_str.split(';')
    for attribute_pair in attribute_pairs:
        attribute_pair = attribute_pair.strip()
        first_space = attribute_pair.find(' ')
        if first_space <= 0:
            continue

        key = attribute_pair[:first_space]
        value = attribute_pair[first_space + 1:]
        if value[0] == '"' and value[-1] == '"':
            # remove quotes
            value = value[1:-1]

        existing_attributes = attributes.get(key)
        if existing_attributes:
            if isinstance(existing_attributes, list):
                existing_attributes.append(value)
            else:
                attributes[key] = [existing_attributes, value]
        else:
            attributes[key] = value

    return attributes


def clean_espresso_gtf(in_handle, out_handle):
    for line in in_handle:
        if line.startswith('#'):
            out_handle.write(line)
            continue

        columns = rmats_long_utils.read_tsv_line(line)
        attribute_string = columns[8]
        attributes = parse_gtf_attributes(attribute_string)
        transcript = attributes.get('transcript_id')
        if not transcript:
            continue

        gene_string = attributes.get('gene_id')
        if not gene_string or gene_string == 'NA':
            continue

        genes = gene_string.split(',')
        if len(genes) == 1:
            out_handle.write(line)
            continue

        for gene_i, gene in enumerate(genes):
            new_transcript_id = '{}_{}'.format(transcript, gene_i)
            attributes['transcript_id'] = new_transcript_id
            attributes['gene_id'] = gene
            columns[8] = rmats_long_utils.format_gtf_attributes(attributes)
            rmats_long_utils.write_tsv_line(out_handle, columns)


def main():
    args = parse_args()
    with open(args.in_gtf, 'rt') as in_handle:
        with open(args.out_gtf, 'wt') as out_handle:
            clean_espresso_gtf(in_handle, out_handle)


if __name__ == '__main__':
    main()
