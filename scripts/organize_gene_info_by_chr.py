import argparse
import os
import os.path

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Create 1 file per chr with transcript info by gene')
    parser.add_argument('--gtf',
                        required=True,
                        help='The path to a gtf with transcript info')
    parser.add_argument(
        '--out-dir',
        required=True,
        help='The directory to create and where chr files will be written')

    return parser.parse_args()


def parse_gtf_gene(start, end, attributes, line, by_gene):
    gene_id = attributes.get('gene_id')
    if gene_id is None:
        print('gene with no gene_id: {}'.format(line))
        return

    gene_info = rmats_long_utils.try_get_or_set_default(
        by_gene, gene_id, dict())
    gene_name = attributes.get('gene_name')
    gene_info['name'] = gene_name
    gene_info['start'] = start
    gene_info['end'] = end


def parse_gtf_transcript(start, end, strand, attributes, line, by_gene):
    gene_id = attributes.get('gene_id')
    if gene_id is None:
        print('transcript with no gene_id: {}'.format(line))
        return

    transcript_id = attributes.get('transcript_id')
    if transcript_id is None:
        print('transcript with no transcript_id: {}'.format(line))
        return

    gene_info = rmats_long_utils.try_get_or_set_default(
        by_gene, gene_id, dict())
    by_transcript = rmats_long_utils.try_get_or_set_default(
        gene_info, 'transcripts', dict())
    transcript_info = rmats_long_utils.try_get_or_set_default(
        by_transcript, transcript_id, dict())
    transcript_name = attributes.get('transcript_name')
    transcript_info['name'] = transcript_name
    transcript_info['strand'] = strand
    transcript_info['start'] = start
    transcript_info['end'] = end


def parse_gtf_exon(start, end, attributes, line, by_gene):
    gene_id = attributes.get('gene_id')
    if gene_id is None:
        print('exon with no gene_id: {}'.format(line))
        return

    transcript_id = attributes.get('transcript_id')
    if transcript_id is None:
        print('exon with no transcript_id: {}'.format(line))
        return

    gene_info = rmats_long_utils.try_get_or_set_default(
        by_gene, gene_id, dict())
    by_transcript = rmats_long_utils.try_get_or_set_default(
        gene_info, 'transcripts', dict())
    transcript_info = rmats_long_utils.try_get_or_set_default(
        by_transcript, transcript_id, dict())
    exons = rmats_long_utils.try_get_or_set_default(transcript_info, 'exons',
                                                    list())
    exons.append((start, end))


def parse_gtf(path):
    transcripts_by_gene_by_chr = dict()
    with open(path, 'rt') as handle:
        for line in handle:
            parsed = rmats_long_utils.parse_gtf_line(line)
            if parsed is None:
                continue

            chr_name = parsed['chr']
            by_gene = rmats_long_utils.try_get_or_set_default(
                transcripts_by_gene_by_chr, chr_name, dict())
            start = parsed['start']
            end = parsed['end']
            strand = parsed['strand']
            attributes = parsed['attributes']
            feature = parsed['feature']
            if feature == 'gene':
                parse_gtf_gene(start, end, attributes, line, by_gene)
            elif feature == 'transcript':
                parse_gtf_transcript(start, end, strand, attributes, line,
                                     by_gene)
            elif feature == 'exon':
                parse_gtf_exon(start, end, attributes, line, by_gene)

    return transcripts_by_gene_by_chr


# The mapping avoids any issues with special characters in chr names
def write_chr_name_mapping(names, out_dir):
    mapping_path = rmats_long_utils.get_chr_id_mapping_path(out_dir)
    chr_to_id = dict()
    with open(mapping_path, 'wt') as handle:
        rmats_long_utils.write_tsv_line(handle, ['name', 'id'])
        for name_i, name in enumerate(names):
            id_str = str(name_i)
            rmats_long_utils.write_tsv_line(handle, [name, id_str])
            chr_to_id[name] = id_str

    return chr_to_id


def create_sjs_string(exons):
    if len(exons) < 2:
        return None

    exons.sort()
    sjs = list()
    prev_end = exons[0][1]
    for start, end in exons[1:]:
        sj = '{}-{}'.format(prev_end, start)
        sjs.append(sj)
        prev_end = end

    return ','.join(sjs)


def write_transcript_output_for_a_gene(by_transcript, handle):
    to_sort = list()
    for transcript_id, transcript_info in by_transcript.items():
        start = transcript_info['start']
        end = transcript_info['end']
        to_sort.append((start, end, transcript_id))

    to_sort.sort()
    for start, end, transcript_id in to_sort:
        transcript_info = by_transcript[transcript_id]
        name = transcript_info.get('name', '')
        if name is None:
            name = ''

        strand = transcript_info['strand']
        exons = transcript_info.get('exons')
        if exons is None:
            continue

        sjs_string = create_sjs_string(exons)
        if sjs_string is None:
            sjs_string = ''

        rmats_long_utils.write_tsv_line(
            handle,
            [str(start),
             str(end), sjs_string, strand, transcript_id, name])


def get_gene_start_and_end_from_transcripts(transcripts):
    min_start = None
    max_end = None
    for transcript in transcripts.values():
        start = transcript['start']
        end = transcript['end']
        if min_start is None:
            min_start = start
            max_end = end
            continue

        if start < min_start:
            min_start = start

        if end > max_end:
            max_end = end

    return min_start, max_end


def write_output_for_chr(by_gene, handle):
    to_sort = list()
    for gene_id, gene_info in by_gene.items():
        start = gene_info.get('start')
        end = gene_info.get('end')
        if start is None or end is None:
            transcripts = gene_info.get('transcripts')
            if transcripts is None:
                continue

            start, end = get_gene_start_and_end_from_transcripts(transcripts)

        to_sort.append((start, end, gene_id))

    to_sort.sort()
    for start, end, gene_id in to_sort:
        gene_info = by_gene[gene_id]
        transcripts = gene_info.get('transcripts')
        if transcripts is None:
            continue

        name = gene_info.get('name', '')
        if name is None:
            name = ''

        rmats_long_utils.write_tsv_line(
            handle, [str(start), str(end), gene_id, name])
        write_transcript_output_for_a_gene(transcripts, handle)


def write_output(transcripts_by_gene_by_chr, out_dir):
    chr_names = sorted(transcripts_by_gene_by_chr.keys())
    chr_to_id = write_chr_name_mapping(chr_names, out_dir)
    for chr_name in chr_names:
        by_gene = transcripts_by_gene_by_chr[chr_name]
        chr_id = chr_to_id[chr_name]
        out_path = rmats_long_utils.get_chr_file_path(out_dir, chr_id)
        with open(out_path, 'wt') as handle:
            write_output_for_chr(by_gene, handle)


def main():
    args = parse_args()
    out_dir = os.path.abspath(args.out_dir)
    rmats_long_utils.create_output_dir(out_dir, check_empty=True)
    transcripts_by_gene_by_chr = parse_gtf(args.gtf)
    write_output(transcripts_by_gene_by_chr, out_dir)
    print('organize_gene_info_by_chr.py finished')


if __name__ == '__main__':
    main()
