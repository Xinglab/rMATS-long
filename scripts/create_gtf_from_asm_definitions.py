import argparse
import os
import os.path

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Create a .gtf file with ASM definitions')
    parser.add_argument(
        '--event-dir',
        required=True,
        help='The output directory from detect_splicing_events.py')
    parser.add_argument('--out-gtf',
                        required=True,
                        help='The output .gtf file to create')

    args = parser.parse_args()
    return args


def find_asm_files(event_dir):
    files_by_chr = dict()
    chr_mapping = rmats_long_utils.parse_chr_mapping(event_dir)
    id_to_chr = chr_mapping['id_to_chr']
    for chr_id, chr_name in id_to_chr.items():
        chr_file_path = rmats_long_utils.get_chr_file_path(event_dir, chr_id)
        if os.path.exists(chr_file_path):
            files_by_chr[chr_name] = chr_file_path

    return files_by_chr


def write_transcript_line(handle, chr_name, exons, strand, gene_id, asm_id,
                          isoform_id):
    feature = 'transcript'
    start = exons[0][0]
    end = exons[-1][1]
    attributes = dict()
    attributes['gene_id'] = gene_id
    attributes['asm_id'] = asm_id
    attributes['transcript_id'] = isoform_id
    rmats_long_utils.write_gtf_line(handle, chr_name, feature, start, end,
                                    strand, attributes)


def write_exon_line(handle, chr_name, exon, strand, gene_id, asm_id,
                    isoform_id):
    feature = 'exon'
    start = exon[0]
    end = exon[1]
    attributes = dict()
    attributes['gene_id'] = gene_id
    attributes['asm_id'] = asm_id
    attributes['transcript_id'] = isoform_id
    rmats_long_utils.write_gtf_line(handle, chr_name, feature, start, end,
                                    strand, attributes)


def write_gtf_lines_with_chr_file_handle(chr_name, chr_handle, gtf_handle):
    gene_id = None
    gene_i = None
    for row in rmats_long_utils.row_iterator_for_tsv_with_header(chr_handle):
        row_gene_id = row['gene_id']
        if row_gene_id != '':
            gene_id = row_gene_id
            gene_i = row['gene_i']
        else:
            row_gene_i = row['gene_i']
            if row_gene_i != gene_i:
                raise Exception('Unexpected gene_i: {}'.format(row))

        asm_id = row['asm_id']
        strand = row['strand']
        isoform_ids_str = row['isoform_ids']
        isoform_ids = rmats_long_utils.parse_isoform_ids_str(isoform_ids_str)
        isoforms = row['isoforms']
        parsed_isoforms = rmats_long_utils.parse_isoforms_str(isoforms)
        for isoform_i, isoform in enumerate(parsed_isoforms):
            isoform_id = isoform_ids[isoform_i]
            is_start = True
            exons = list()
            for coord in isoform:
                if is_start:
                    exons.append([coord, None])
                    is_start = False
                else:
                    exons[-1][1] = coord
                    is_start = True

            write_transcript_line(gtf_handle, chr_name, exons, strand, gene_id,
                                  asm_id, isoform_id)
            for exon in exons:
                write_exon_line(gtf_handle, chr_name, exon, strand, gene_id,
                                asm_id, isoform_id)


def write_gtf_lines_using_chr_files(files_by_chr, gtf_handle):
    for chr_name, path in files_by_chr.items():
        with open(path, 'rt') as chr_handle:
            write_gtf_lines_with_chr_file_handle(chr_name, chr_handle,
                                                 gtf_handle)


def create_gtf_from_asm_definitions(event_dir, out_gtf):
    files_by_chr = find_asm_files(event_dir)
    with open(out_gtf, 'wt') as gtf_handle:
        write_gtf_lines_using_chr_files(files_by_chr, gtf_handle)


def main():
    args = parse_args()
    create_gtf_from_asm_definitions(args.event_dir, args.out_gtf)
    print('create_gtf_from_asm_definitions.py finished')


if __name__ == '__main__':
    main()
