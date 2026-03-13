import argparse
import multiprocessing
import os
import os.path

from rmats_long import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Aggregate read to isoform compatibility for'
        ' pairs of ASMs in a gene')
    parser.add_argument(
        '--samples-tsv',
        required=True,
        help='The path to a file where the first column is the sample name')
    parser.add_argument(
        '--counts-dir',
        required=True,
        help='The directory of read isoform compatibility files')
    parser.add_argument(
        '--event-dir',
        required=True,
        help='The output directory from detect_splicing_events.py')
    parser.add_argument(
        '--out-dir',
        required=True,
        help='The output directory for aggregated read count files')
    parser.add_argument('--num-threads',
                        default=1,
                        type=int,
                        help='how many threads to use')
    parser.add_argument('--max-isoforms-per-asm',
                        type=int,
                        help='how many threads to use')

    args = parser.parse_args()
    return args


def read_samples_tsv(samples_tsv):
    samples = list()
    with open(samples_tsv, 'rt') as handle:
        for line in handle:
            columns = rmats_long_utils.read_tsv_line(line)
            sample = columns[0]
            samples.append(sample)

    return samples


def write_out_header(samples, out_handle):
    columns = ['asm_pair', 'asm_1_isoforms', 'asm_2_isoforms']
    columns.extend(samples)
    rmats_long_utils.write_tsv_line(out_handle, columns)


def get_isoforms_for_asm_read(gene_file_headers, line, asm, in_handle):
    isoforms = list()
    read_id = None
    sample = None
    if not line:
        line = in_handle.readline()

    while line:
        columns = rmats_long_utils.read_tsv_line(line)
        row = dict(zip(gene_file_headers, columns))
        found_asm = row['asm']
        if found_asm != asm:
            return {
                'line': line,
                'read': read_id,
                'isoforms': isoforms,
                'sample': sample
            }

        found_read = int(row['read'])
        if read_id is None:
            read_id = found_read

        if found_read != read_id:
            return {
                'line': line,
                'read': read_id,
                'isoforms': isoforms,
                'sample': sample
            }

        sample = row['sample']
        isoform = row['isoform']
        isoforms.append(isoform)
        line = in_handle.readline()

    return {
        'line': None,
        'read': read_id,
        'isoforms': isoforms,
        'sample': sample
    }


def update_combo_counts(sample, isoforms_1, isoforms_2, counts_by_sample_combo,
                        all_combos):
    counts_by_combo = rmats_long_utils.try_get_or_set_default(
        counts_by_sample_combo, sample, dict())

    isoforms_1.sort()
    isoforms_2.sort()
    combo = (tuple(isoforms_1), tuple(isoforms_2))
    all_combos.add(combo)
    old_count = counts_by_combo.get(combo, 0)
    counts_by_combo[combo] = old_count + 1


def output_combo_counts_for_asm_pair(asm_1, asm_2, samples,
                                     counts_by_sample_combo, all_combos,
                                     out_handle):
    pair_string = '{},{}'.format(asm_1, asm_2)
    sorted_combos = sorted(list(all_combos))
    for combo in sorted_combos:
        asm_1_isoforms, asm_2_isoforms = combo
        comma_sep_1 = ','.join(asm_1_isoforms)
        comma_sep_2 = ','.join(asm_2_isoforms)
        columns = [pair_string, comma_sep_1, comma_sep_2]
        for sample in samples:
            counts_by_combo = counts_by_sample_combo.get(sample, dict())
            count = counts_by_combo.get(combo, 0)
            columns.append(str(count))

        rmats_long_utils.write_tsv_line(out_handle, columns)


def check_if_asms_overlap(asm_1_details, asm_2_details):
    start_1 = asm_1_details['start']
    end_1 = asm_1_details['end']
    if asm_1_details['strand'] == '-':
        start_1, end_1 = end_1, start_1

    start_2 = asm_2_details['start']
    end_2 = asm_2_details['end']
    if asm_2_details['strand'] == '-':
        start_2, end_2 = end_2, start_2

    if (start_1 is not None) and (end_2 is not None) and (start_1 > end_2):
        return False

    if (start_2 is not None) and (end_1 is not None) and (start_2 > end_1):
        return False

    return True


def get_and_output_combo_counts_for_asm_pair(gene_file_headers, samples, asm_1,
                                             asm_2, asm_1_offset, asm_2_offset,
                                             asm_1_details, asm_2_details,
                                             in_handle_1, in_handle_2,
                                             out_handle):
    asms_overlap = check_if_asms_overlap(asm_1_details, asm_2_details)
    if asms_overlap:
        return

    counts_by_sample_combo = dict()
    all_combos = set()
    asm_1_result = None
    asm_2_result = None
    in_handle_1.seek(asm_1_offset)
    in_handle_2.seek(asm_2_offset)
    asm_1_line = in_handle_1.readline()
    asm_2_line = in_handle_2.readline()
    while asm_1_line and asm_2_line:
        if asm_1_result is None:
            asm_1_result = get_isoforms_for_asm_read(gene_file_headers,
                                                     asm_1_line, asm_1,
                                                     in_handle_1)
        if asm_2_result is None:
            asm_2_result = get_isoforms_for_asm_read(gene_file_headers,
                                                     asm_2_line, asm_2,
                                                     in_handle_2)

        asm_1_read = asm_1_result['read']
        asm_2_read = asm_2_result['read']
        if asm_1_read is None:
            asm_1_line = None
            asm_1_result = None
            break

        if asm_2_read is None:
            asm_2_line = None
            asm_2_result = None
            break

        if asm_1_read < asm_2_read:
            update_combo_counts(asm_1_result['sample'],
                                asm_1_result['isoforms'], list(),
                                counts_by_sample_combo, all_combos)
            asm_1_line = asm_1_result['line']
            asm_1_result = None
        elif asm_2_read < asm_1_read:
            update_combo_counts(asm_2_result['sample'], list(),
                                asm_2_result['isoforms'],
                                counts_by_sample_combo, all_combos)
            asm_2_line = asm_2_result['line']
            asm_2_result = None
        else:
            update_combo_counts(asm_1_result['sample'],
                                asm_1_result['isoforms'],
                                asm_2_result['isoforms'],
                                counts_by_sample_combo, all_combos)
            asm_1_line = asm_1_result['line']
            asm_2_line = asm_2_result['line']
            asm_1_result = None
            asm_2_result = None

    while asm_1_line:
        if asm_1_result is None:
            asm_1_result = get_isoforms_for_asm_read(gene_file_headers,
                                                     asm_1_line, asm_1,
                                                     in_handle_1)

        asm_1_read = asm_1_result['read']
        if asm_1_read is None:
            asm_1_line = None
            asm_1_result = None
            break

        update_combo_counts(asm_1_result['sample'], asm_1_result['isoforms'],
                            list(), counts_by_sample_combo, all_combos)
        asm_1_line = asm_1_result['line']
        asm_1_result = None

    while asm_2_line:
        if asm_2_result is None:
            asm_2_result = get_isoforms_for_asm_read(gene_file_headers,
                                                     asm_2_line, asm_2,
                                                     in_handle_2)

        asm_2_read = asm_2_result['read']
        if asm_2_read is None:
            asm_2_line = None
            asm_2_result = None
            break

        update_combo_counts(asm_2_result['sample'], list(),
                            asm_2_result['isoforms'], counts_by_sample_combo,
                            all_combos)
        asm_2_line = asm_2_result['line']
        asm_2_result = None

    output_combo_counts_for_asm_pair(asm_1, asm_2, samples,
                                     counts_by_sample_combo, all_combos,
                                     out_handle)


def update_asm_details(event_headers, asm, event_handle_and_line, asm_details):
    target_asm = rmats_long_utils.parse_asm_id(asm)
    target_event_i = target_asm['event_i']
    line = event_handle_and_line.read_line()
    while line:
        columns = rmats_long_utils.read_tsv_line(line)
        row = dict(zip(event_headers, columns))
        found_asm = rmats_long_utils.parse_asm_id(row['asm_id'])
        found_event_i = found_asm['event_i']
        if found_event_i > target_event_i:
            break

        line = event_handle_and_line.clear_and_read_line()
        if found_event_i < target_event_i:
            continue

        if row['start'] == 'source':
            start = None
        else:
            start = int(row['start'])

        if row['end'] == 'sink':
            end = None
        else:
            end = int(row['end'])

        num_isoforms = len(row['isoform_ids'].split(';'))
        details = {
            'start': start,
            'end': end,
            'strand': row['strand'],
            'num_isoforms': num_isoforms,
        }
        asm_details[asm] = details


def process_gene_reads(gene_file_headers, event_headers, gene_file_path,
                       samples, max_isoforms_per_asm, event_handle_and_line,
                       out_handle):
    asms = list()
    asm_offsets = list()
    asm_details = dict()
    prev_asm = None
    # Open two handles, one for each asm.
    # Separate handles may have better performance than constantly seeking.
    with open(gene_file_path, 'rt') as in_handle_1:
        with open(gene_file_path, 'rt') as in_handle_2:
            offset = in_handle_1.tell()
            line = in_handle_1.readline()
            while line:
                columns = rmats_long_utils.read_tsv_line(line)
                row = dict(zip(gene_file_headers, columns))
                asm = row['asm']
                if (prev_asm is None) or (asm != prev_asm):
                    update_asm_details(event_headers, asm,
                                       event_handle_and_line, asm_details)
                    details = asm_details[asm]
                    over_max = (max_isoforms_per_asm
                                and (details['num_isoforms']
                                     > max_isoforms_per_asm))
                    has_source_or_sink = None in [
                        details['start'], details['end']
                    ]
                    if not (over_max or has_source_or_sink):
                        asms.append(asm)
                        asm_offsets.append(offset)

                prev_asm = asm
                offset = in_handle_1.tell()
                line = in_handle_1.readline()

            num_asms = len(asms)
            for asm_1_i, asm_1 in enumerate(asms):
                asm_1_offset = asm_offsets[asm_1_i]
                asm_1_details = asm_details[asm_1]
                for asm_2_i in range(asm_1_i + 1, num_asms):
                    asm_2 = asms[asm_2_i]
                    asm_2_offset = asm_offsets[asm_2_i]
                    asm_2_details = asm_details[asm_2]
                    get_and_output_combo_counts_for_asm_pair(
                        gene_file_headers, samples, asm_1, asm_2, asm_1_offset,
                        asm_2_offset, asm_1_details, asm_2_details,
                        in_handle_1, in_handle_2, out_handle)


def count_reads_for_asm_pairs_with_handles(samples, max_isoforms_per_asm,
                                           counts_handle, event_handle,
                                           tmp_path, out_handle):
    tmp_headers = ['asm', 'read', 'sample', 'isoform']
    tmp_handle = None
    event_handle_and_line = rmats_long_utils.HandleAndNextLine(event_handle)
    event_header_line = event_handle_and_line.read_line()
    event_handle_and_line.clear_and_read_line()  # move past the header line
    event_headers = rmats_long_utils.read_tsv_line(event_header_line)
    try:
        prev_gene_id = None
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(
                counts_handle):
            asm_id = row['asm_id']
            gene_id = row['gene_id']
            read_id = row['read_id']
            sample_id = row['sample_id']
            isoform_id = row['isoform_id']
            if not gene_id:
                gene_id = prev_gene_id

            if gene_id != prev_gene_id:
                if tmp_handle is not None:
                    tmp_handle.close()
                    tmp_handle = None
                    process_gene_reads(tmp_headers, event_headers, tmp_path,
                                       samples, max_isoforms_per_asm,
                                       event_handle_and_line, out_handle)

                tmp_handle = open(tmp_path, 'wt')

            columns = [asm_id, read_id, sample_id, isoform_id]
            rmats_long_utils.write_tsv_line(tmp_handle, columns)
            prev_gene_id = gene_id

        if tmp_handle is not None:
            tmp_handle.close()
            tmp_handle = None
            process_gene_reads(tmp_headers, event_headers, tmp_path, samples,
                               max_isoforms_per_asm, event_handle_and_line,
                               out_handle)
    finally:
        if tmp_handle is not None:
            tmp_handle.close()
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def get_inputs_for_chrs(counts_dir, event_dir, out_dir):
    chr_infos = list()
    counts_names = os.listdir(counts_dir)
    for name in counts_names:
        if not name.startswith('chr_id_'):
            continue

        counts_path = os.path.join(counts_dir, name)
        event_path = os.path.join(event_dir, name)
        out_path = os.path.join(out_dir, name)
        counts_size = rmats_long_utils.get_file_size(counts_path)
        chr_info = {
            'counts_path': counts_path,
            'event_path': event_path,
            'out_path': out_path,
            'counts_size': counts_size,
        }
        chr_infos.append(chr_info)

    chr_infos.sort(key=lambda x: x['counts_size'], reverse=True)
    return chr_infos


def count_reads_for_asm_pairs_for_chr(samples, max_isoforms_per_asm,
                                      counts_path, event_path, out_path):
    tmp_path = '{}.tmp'.format(out_path)
    with open(out_path, 'wt') as out_handle:
        write_out_header(samples, out_handle)
        with open(counts_path, 'rt') as counts_handle:
            with open(event_path, 'rt') as event_handle:
                count_reads_for_asm_pairs_with_handles(samples,
                                                       max_isoforms_per_asm,
                                                       counts_handle,
                                                       event_handle, tmp_path,
                                                       out_handle)


def count_reads_for_asm_pairs_thread(samples, max_isoforms_per_asm,
                                     thread_inputs):
    while True:
        arguments = rmats_long_utils.try_get_from_queue_with_short_wait(
            thread_inputs)
        if arguments is None:
            return

        counts_path = arguments['counts_path']
        event_path = arguments['event_path']
        out_path = arguments['out_path']
        count_reads_for_asm_pairs_for_chr(samples, max_isoforms_per_asm,
                                          counts_path, event_path, out_path)


def count_reads_for_asm_pairs(num_threads, samples, max_isoforms_per_asm,
                              counts_dir, event_dir, out_dir):
    chr_infos = get_inputs_for_chrs(counts_dir, event_dir, out_dir)
    num_jobs = len(chr_infos)
    threads = list()
    thread_inputs = multiprocessing.Queue(num_jobs)
    for chr_info in chr_infos:
        arguments = {
            'counts_path': chr_info['counts_path'],
            'event_path': chr_info['event_path'],
            'out_path': chr_info['out_path'],
        }
        thread_inputs.put(arguments)

    for _ in range(num_threads):
        thread = multiprocessing.Process(
            target=count_reads_for_asm_pairs_thread,
            args=(samples, max_isoforms_per_asm, thread_inputs))
        threads.append(thread)
        thread.start()

    rmats_long_utils.join_threads_and_raise_if_error(threads)


def main():
    args = parse_args()
    rmats_long_utils.create_output_dir(args.out_dir, check_empty=True)
    samples = read_samples_tsv(args.samples_tsv)
    count_reads_for_asm_pairs(args.num_threads, samples,
                              args.max_isoforms_per_asm, args.counts_dir,
                              args.event_dir, args.out_dir)

    print('count_reads_for_asm_pairs.py finished')


if __name__ == '__main__':
    main()
