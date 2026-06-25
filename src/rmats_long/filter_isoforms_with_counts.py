import argparse
import multiprocessing
import os
import os.path
import shutil

from rmats_long import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Filter ASM isoforms based on compatible reads')
    parser.add_argument(
        '--event-dir',
        required=True,
        help='The output directory from detect_splicing_events.py')
    parser.add_argument(
        '--asm-counts-dir',
        required=True,
        help='The output directory from count_reads_for_asms.py')
    parser.add_argument(
        '--out-event-dir',
        required=True,
        help='The path for the output filtered event directory')
    parser.add_argument(
        '--out-counts-dir',
        required=True,
        help='The path for the output filtered read assignments')
    parser.add_argument(
        '--all-sjs-count',
        type=int,
        default=1,
        help=('The minimum number of compatible reads that use all splice'
              ' junctions of the isoform (default %(default)s)'))
    parser.add_argument(
        '--num-threads',
        type=int,
        default=1,
        help='The number of threads to use (default %(default)s)')

    args = parser.parse_args()
    if args.all_sjs_count <= 0:
        raise Exception('--all-sjs-count must be > 0: {}'.format(
            args.all_sjs_count))

    return args


def filter_isoforms_for_row(row, kept_isoforms):
    isoforms_str = row['isoforms']
    isoforms = rmats_long_utils.parse_isoforms_str(isoforms_str)
    isoform_ids_str = row['isoform_ids']
    isoform_ids = rmats_long_utils.parse_isoform_ids_str(isoform_ids_str)

    filtered_isoforms = list()
    filtered_isoform_ids = list()
    for isoform_i, isoform in enumerate(isoforms):
        isoform_id = isoform_ids[isoform_i]
        if isoform_id not in kept_isoforms:
            continue

        filtered_isoforms.append(isoform)
        filtered_isoform_ids.append(isoform_id)

    return {'isoforms': filtered_isoforms, 'isoform_ids': filtered_isoform_ids}


def filter_event_rows_for_asm(asm_id, kept_isoforms, event_headers,
                              event_state, in_event_handle, out_event_handle):
    parsed_asm_id = rmats_long_utils.parse_asm_id(asm_id)
    event_i = parsed_asm_id['event_i']
    while event_state['next_line']:
        columns = rmats_long_utils.read_tsv_line(event_state['next_line'])
        row = dict(zip(event_headers, columns))
        gene_id = row['gene_id']
        if gene_id:
            event_state['gene_id_to_write'] = gene_id

        found_asm_id = row['asm_id']
        parsed_found_asm_id = rmats_long_utils.parse_asm_id(found_asm_id)
        found_event_i = parsed_found_asm_id['event_i']
        if found_event_i > event_i:
            break

        event_state['next_line'] = in_event_handle.readline()
        if found_event_i < event_i:
            continue

        filtered_results = filter_isoforms_for_row(row, kept_isoforms)
        filtered_isoforms = filtered_results['isoforms']
        filtered_isoform_ids = filtered_results['isoform_ids']
        if not filtered_isoforms:
            continue

        isoforms_str = rmats_long_utils.format_isoforms_str(filtered_isoforms)
        isoform_ids_str = rmats_long_utils.format_isoform_ids_str(
            filtered_isoform_ids)
        out_columns = list()
        for header in event_headers:
            if header == 'isoforms':
                out_columns.append(isoforms_str)
            elif header == 'isoform_ids':
                out_columns.append(isoform_ids_str)
            elif header == 'gene_id':
                out_columns.append(event_state['gene_id_to_write'])
                event_state['gene_id_to_write'] = ''
            else:
                value = row[header]
                out_columns.append(value)

        rmats_long_utils.write_tsv_line(out_event_handle, out_columns)


def filter_count_rows_for_asm(asm_id, kept_isoforms, count_headers,
                              count_state, in_count_handle, out_count_handle):
    parsed_asm_id = rmats_long_utils.parse_asm_id(asm_id)
    event_i = parsed_asm_id['event_i']
    while count_state['next_line']:
        columns = rmats_long_utils.read_tsv_line(count_state['next_line'])
        row = dict(zip(count_headers, columns))
        gene_id = row['gene_id']
        if gene_id:
            count_state['gene_id_to_write'] = gene_id

        found_asm_id = row['asm_id']
        parsed_found_asm_id = rmats_long_utils.parse_asm_id(found_asm_id)
        found_event_i = parsed_found_asm_id['event_i']
        if found_event_i > event_i:
            break

        count_state['next_line'] = in_count_handle.readline()
        if found_event_i < event_i:
            continue

        isoform_id = row['isoform_id']
        if isoform_id not in kept_isoforms:
            continue

        out_columns = list()
        for header in count_headers:
            if header == 'gene_id':
                out_columns.append(count_state['gene_id_to_write'])
                count_state['gene_id_to_write'] = ''
            else:
                value = row[header]
                out_columns.append(value)

        rmats_long_utils.write_tsv_line(out_count_handle, out_columns)


def filter_isoforms_with_counts_with_handles(min_all_sjs_count,
                                             in_event_handle, in_count_handle,
                                             all_sjs_handle, out_event_handle,
                                             out_count_handle):
    event_header_line = in_event_handle.readline()
    event_headers = rmats_long_utils.read_tsv_line(event_header_line)
    rmats_long_utils.write_tsv_line(out_event_handle, event_headers)
    next_event_line = in_event_handle.readline()

    count_header_line = in_count_handle.readline()
    count_headers = rmats_long_utils.read_tsv_line(count_header_line)
    next_count_line = in_count_handle.readline()
    rmats_long_utils.write_tsv_line(out_count_handle, count_headers)

    prev_asm_id = None
    kept_isoforms = set()
    # The gene_id is only written once.
    # Additional rows with the same gene_id have it set to ''
    event_state = {'next_line': next_event_line, 'gene_id_to_write': ''}
    count_state = {'next_line': next_count_line, 'gene_id_to_write': ''}
    for row in rmats_long_utils.row_iterator_for_tsv_with_header(
            all_sjs_handle):
        asm_id = row['asm_id']
        if asm_id != prev_asm_id:
            if prev_asm_id is not None:
                filter_event_rows_for_asm(prev_asm_id, kept_isoforms,
                                          event_headers, event_state,
                                          in_event_handle, out_event_handle)
                filter_count_rows_for_asm(prev_asm_id, kept_isoforms,
                                          count_headers, count_state,
                                          in_count_handle, out_count_handle)

            prev_asm_id = asm_id
            kept_isoforms = set()

        isoform_id = row['isoform_id']
        all_sjs_count_str = row['all_sjs_matched_count']
        all_sjs_count = int(all_sjs_count_str)
        if all_sjs_count < min_all_sjs_count:
            continue

        kept_isoforms.add(isoform_id)

    if prev_asm_id is not None:
        filter_event_rows_for_asm(prev_asm_id, kept_isoforms, event_headers,
                                  event_state, in_event_handle,
                                  out_event_handle)
        filter_count_rows_for_asm(prev_asm_id, kept_isoforms, count_headers,
                                  count_state, in_count_handle,
                                  out_count_handle)


def copy_header_line(in_path, out_path):
    with open(in_path, 'rt') as in_handle:
        header_line = in_handle.readline()

    with open(out_path, 'wt') as out_handle:
        out_handle.write(header_line)


def filter_isoforms_with_counts_thread(min_all_sjs_count, input_queue):
    while True:
        input_data = input_queue.get()
        if input_data is None:
            return

        in_event_path = input_data['in_event_path']
        out_event_path = input_data['out_event_path']
        in_graph_path = input_data['in_graph_path']
        out_graph_path = input_data['out_graph_path']
        in_count_path = input_data['in_count_path']
        out_count_path = input_data['out_count_path']
        all_sjs_path = input_data['all_sjs_path']

        if not os.path.exists(in_count_path):
            if os.path.exists(in_event_path):
                copy_header_line(in_event_path, out_event_path)

            continue

        shutil.copy(in_graph_path, out_graph_path)
        with open(in_event_path, 'rt') as in_event_handle:
            with open(in_count_path, 'rt') as in_count_handle:
                with open(all_sjs_path, 'rt') as all_sjs_handle:
                    with open(out_event_path, 'wt') as out_event_handle:
                        with open(out_count_path, 'wt') as out_count_handle:
                            filter_isoforms_with_counts_with_handles(
                                min_all_sjs_count, in_event_handle,
                                in_count_handle, all_sjs_handle,
                                out_event_handle, out_count_handle)


def filter_isoforms_with_counts(event_dir, counts_dir, out_event_dir,
                                out_counts_dir, min_all_sjs_count,
                                num_threads):
    chr_id_mapping = rmats_long_utils.parse_chr_mapping(event_dir)
    chr_to_id = chr_id_mapping['chr_to_id']
    rmats_long_utils.copy_chr_name_mapping(event_dir, out_event_dir)
    rmats_long_utils.copy_chr_name_mapping(counts_dir, out_counts_dir)

    num_chr_names = len(chr_to_id)
    num_threads = min(num_threads, num_chr_names)
    input_queue = multiprocessing.Queue(num_chr_names + num_threads)
    threads = list()
    for _ in range(num_threads):
        thread = multiprocessing.Process(
            target=filter_isoforms_with_counts_thread,
            args=(min_all_sjs_count, input_queue))
        threads.append(thread)
        thread.start()

    for chr_id in chr_to_id.values():
        in_event_path = rmats_long_utils.get_chr_file_path(event_dir, chr_id)
        out_event_path = rmats_long_utils.get_chr_file_path(
            out_event_dir, chr_id)
        in_graph_path = rmats_long_utils.get_graph_file_path(event_dir, chr_id)
        out_graph_path = rmats_long_utils.get_graph_file_path(
            out_event_dir, chr_id)
        in_count_path = rmats_long_utils.get_chr_file_path(counts_dir, chr_id)
        out_count_path = rmats_long_utils.get_chr_file_path(
            out_counts_dir, chr_id)
        all_sjs_path = rmats_long_utils.get_all_sjs_path(counts_dir, chr_id)
        input_queue.put({
            'in_event_path': in_event_path,
            'out_event_path': out_event_path,
            'in_graph_path': in_graph_path,
            'out_graph_path': out_graph_path,
            'in_count_path': in_count_path,
            'out_count_path': out_count_path,
            'all_sjs_path': all_sjs_path
        })

    # Signal the end of inputs
    for _ in range(num_threads):
        input_queue.put(None)

    rmats_long_utils.join_threads_and_raise_if_error(threads)


def main():
    args = parse_args()
    out_event_dir = os.path.abspath(args.out_event_dir)
    rmats_long_utils.create_output_dir(out_event_dir, check_empty=True)
    out_counts_dir = os.path.abspath(args.out_counts_dir)
    rmats_long_utils.create_output_dir(out_counts_dir, check_empty=True)
    filter_isoforms_with_counts(args.event_dir, args.asm_counts_dir,
                                out_event_dir, out_counts_dir,
                                args.all_sjs_count, args.num_threads)

    print('filter_isoforms_with_counts.py finished')


if __name__ == '__main__':
    main()
