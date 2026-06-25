import argparse
import os
import os.path

from rmats_long import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Create alignment files with a read corresponding to'
                     ' each isoform'))
    parser.add_argument(
        '--event-dir',
        required=True,
        help='The output directory from detect_splicing_events.py')
    parser.add_argument(
        '--out-dir',
        required=True,
        help='The directory to create and where new files will be written')

    return parser.parse_args()


def is_isoform_start(event_start, coord):
    return (event_start == 'source') or (str(coord) == event_start)


def is_isoform_end(event_end, coord, coord_i, final_i):
    if event_end == 'sink':
        if coord_i == final_i:
            return True
    elif str(coord) == event_end:
        return True

    return False


def get_isoform_start_end_and_sjs(strand, event_start, event_end, isoform):
    start = None
    end = None
    sjs = list()

    if strand == '-':
        isoform = list(reversed(isoform))

    final_i = len(isoform) - 1
    at_exon_start = False
    for coord_i, coord in enumerate(isoform):
        at_exon_start = not at_exon_start
        if start is None:
            if is_isoform_start(event_start, coord):
                start = coord
                if is_isoform_end(event_end, coord, coord_i, final_i):
                    end = coord
                    break

                if not at_exon_start:
                    sjs.append([coord, None])

            continue

        if is_isoform_end(event_end, coord, coord_i, final_i):
            end = coord
            if at_exon_start:
                sjs[-1][-1] = coord

            break

        if at_exon_start:
            sjs[-1][-1] = coord
        else:
            sjs.append([coord, None])

    if strand == '-':
        start, end = end, start
        reversed_sjs = list()
        for sj_start, sj_end in reversed(sjs):
            reversed_sjs.append([sj_end, sj_start])

        sjs = reversed_sjs

    return {'start': start, 'end': end, 'sjs': sjs}


def create_alignments_from_events_with_handles(event_handle, out_handle,
                                               index_handle):
    sample = ''  # leave the sample column empty
    prev_gene_i = None
    for row in rmats_long_utils.row_iterator_for_tsv_with_header(event_handle):
        gene_i = row['gene_i']
        gene_id = row['gene_id']
        strand = row['strand']
        event_start = row['start']
        event_end = row['end']
        isoforms_str = row['isoforms']
        isoform_ids_str = row['isoform_ids']
        isoforms = rmats_long_utils.parse_isoforms_str(isoforms_str)
        isoform_ids = rmats_long_utils.parse_isoform_ids_str(isoform_ids_str)
        if prev_gene_i != gene_i:
            prev_gene_i = gene_i
            offset = out_handle.tell()
            rmats_long_utils.write_tsv_line(index_handle,
                                            [gene_id, str(offset)])
        else:
            gene_id = ''

        for isoform_i, isoform in enumerate(isoforms):
            isoform_id = isoform_ids[isoform_i]
            isoform_details = get_isoform_start_end_and_sjs(
                strand, event_start, event_end, isoform)
            start = isoform_details['start']
            end = isoform_details['end']
            sjs = isoform_details['sjs']
            sjs_string = rmats_long_utils.create_sjs_string(sjs)
            out_columns = [
                gene_i, gene_id, sample,
                str(start),
                str(end), sjs_string, strand, isoform_id
            ]
            rmats_long_utils.write_tsv_line(out_handle, out_columns)
            # only write the gene_id once
            gene_id = ''


def create_alignments_from_events(event_dir, out_dir):
    chr_id_mapping = rmats_long_utils.parse_chr_mapping(event_dir)
    chr_to_id = chr_id_mapping['chr_to_id']
    rmats_long_utils.copy_chr_name_mapping(event_dir, out_dir)
    for chr_id in chr_to_id.values():
        event_path = rmats_long_utils.get_chr_file_path(event_dir, chr_id)
        out_path = rmats_long_utils.get_chr_file_path(out_dir, chr_id)
        index_path = '{}.index'.format(out_path)
        with open(event_path, 'rt') as event_handle:
            with open(out_path, 'wt') as out_handle:
                with open(index_path, 'wt') as index_handle:
                    create_alignments_from_events_with_handles(
                        event_handle, out_handle, index_handle)


def main():
    args = parse_args()
    out_dir = os.path.abspath(args.out_dir)
    rmats_long_utils.create_output_dir(out_dir, check_empty=True)
    create_alignments_from_events(args.event_dir, out_dir)
    print('create_alignments_from_events.py finished')


if __name__ == '__main__':
    main()
