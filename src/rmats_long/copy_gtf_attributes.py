import argparse

from rmats_long import rmats_long_utils

CANONICAL_TAG = 'Ensembl_canonical'


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Copy gene_name and tag Ensembl_canonical attributes'
                     ' to a new .gtf'))
    parser.add_argument('--gencode-gtf',
                        required=True,
                        help='The .gtf with the attributes to copy')
    parser.add_argument('--other-gtf',
                        required=True,
                        help='The .gtf to add attributes to')
    parser.add_argument('--out-gtf', required=True, help='The output .gtf')

    args = parser.parse_args()
    return args


def get_attribute_mappings(gtf):
    gene_id_to_name = dict()
    canonical_transcripts = set()
    with open(gtf, 'rt') as handle:
        for line in handle:
            parsed = rmats_long_utils.parse_gtf_line(line)
            if not parsed:
                continue

            attributes = parsed['attributes']
            gene_id = attributes.get('gene_id')
            gene_name = attributes.get('gene_name')
            transcript_id = attributes.get('transcript_id')
            tags = attributes.get('tag')
            if gene_id and gene_name:
                gene_id_to_name[gene_id] = gene_name

            if not transcript_id:
                continue

            if isinstance(tags, list):
                if CANONICAL_TAG not in tags:
                    continue
            elif tags != CANONICAL_TAG:
                continue

            canonical_transcripts.add(transcript_id)

    return {
        'gene_id_to_name': gene_id_to_name,
        'canonical_transcripts': canonical_transcripts
    }


def add_attributes(gene_id_to_name, canonical_transcripts, in_gtf, out_gtf):
    with open(in_gtf, 'rt') as in_handle:
        with open(out_gtf, 'wt') as out_handle:
            for line in in_handle:
                parsed = rmats_long_utils.parse_gtf_line(line)
                if not parsed:
                    out_handle.write(line)
                    continue

                chr_name = parsed['chr']
                source = parsed['source']
                feature = parsed['feature']
                start = parsed['start']
                end = parsed['end']
                score = parsed['score']
                strand = parsed['strand']
                frame = parsed['frame']
                attributes = parsed['attributes']

                gene_id = attributes.get('gene_id')
                gene_name = gene_id_to_name.get(gene_id)
                if gene_name:
                    attributes['gene_name'] = gene_name

                transcript_id = attributes.get('transcript_id')
                if transcript_id in canonical_transcripts:
                    tags = attributes.get('tag')
                    if not tags:
                        attributes['tag'] = CANONICAL_TAG
                    elif isinstance(tags, list):
                        if CANONICAL_TAG not in tags:
                            tags.append(CANONICAL_TAG)
                    elif tags != CANONICAL_TAG:
                        attributes['tag'] = [tags, CANONICAL_TAG]

                rmats_long_utils.write_gtf_line(out_handle,
                                                chr_name,
                                                feature,
                                                start,
                                                end,
                                                strand,
                                                attributes,
                                                source=source,
                                                score=score,
                                                frame=frame)


def copy_gtf_attributes(gencode_gtf, other_gtf, out_gtf):
    mappings = get_attribute_mappings(gencode_gtf)
    gene_id_to_name = mappings['gene_id_to_name']
    canonical_transcripts = mappings['canonical_transcripts']
    add_attributes(gene_id_to_name, canonical_transcripts, other_gtf, out_gtf)


def main():
    args = parse_args()
    copy_gtf_attributes(args.gencode_gtf, args.other_gtf, args.out_gtf)

    print('copy_gtf_attributes.py finished')


if __name__ == '__main__':
    main()
