import argparse
import multiprocessing
import os
import os.path
import shutil
import subprocess
import tempfile

from rmats_long import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Create 1 file per chr with read info by gene')
    parser.add_argument(
        '--gtf-dir',
        required=True,
        help='The output directory from organize_gene_info_by_chr.py')
    parser.add_argument(
        '--out-dir',
        required=True,
        help='The directory to create and where new files will be written')
    parser.add_argument(
        '--samples-tsv',
        required=True,
        help=('The path to a file where each line has 2 tab separated columns:'
              ' sample name, then a path to an output file from'
              ' simplify_alignment_info.py. A sample name can have multiple'
              ' lines if it has multiple input files.'))
    parser.add_argument('--sort-buffer-size',
                        default='2G',
                        help=('Used for the --buffer-size argument of sort.'
                              ' Default: %(default)s'))
    parser.add_argument(
        '--num-threads',
        type=int,
        default=1,
        help='The number of threads to use (default %(default)s)')

    return parser.parse_args()


def parse_samples_tsv(samples_path):
    sample_files = list()
    with open(samples_path, 'rt') as handle:
        for line in handle:
            columns = rmats_long_utils.read_tsv_line(line)
            if columns == ['']:
                continue  # Skip empty lines

            sample_name = columns[0]
            path = columns[1]
            index = '{}.index'.format(path)
            details = {'sample': sample_name, 'path': path, 'index': index}
            sample_files.append(details)

    return sample_files


def load_offsets_by_chr_from_index(index):
    offsets_by_chr = dict()
    with open(index, 'rt') as handle:
        for line in handle:
            columns = rmats_long_utils.read_tsv_line(line)
            chr_name = columns[0]
            offset_str = columns[1]
            offset = int(offset_str)
            offsets_by_chr[chr_name] = offset

    return offsets_by_chr


def load_indices(sample_files):
    indices_by_path = dict()
    for details in sample_files:
        path = details['path']
        index = details['index']
        offset_by_chr = load_offsets_by_chr_from_index(index)
        indices_by_path[path] = offset_by_chr

    return indices_by_path


def get_sorted_samples(sample_files):
    samples = list()
    for details in sample_files:
        sample = details['sample']
        samples.append(sample)

    samples.sort()
    return samples


def parse_read_line(line):
    parsed = dict()
    columns = rmats_long_utils.read_tsv_line(line)
    parsed['columns'] = columns
    parsed['chr'] = columns[0]
    start_str = columns[1]
    parsed['start'] = int(start_str)
    end_str = columns[2]
    parsed['end'] = int(end_str)
    # sjs_string = columns[3]
    # strand = columns[4]

    return parsed


def initialize_next_read_offsets(sample_files, chr_name, offsets_for_chr):
    for details in sample_files:
        path = details['path']
        offset = offsets_for_chr.get(path)
        details['next_offset'] = offset
        if offset is None:
            continue

        with open(path, 'rt') as handle:
            handle.seek(offset)
            line = handle.readline()
            if not line:
                details['next_offset'] = None
                continue

            parsed = parse_read_line(line)
            found_chr = parsed['chr']
            if found_chr != chr_name:
                raise Exception('Expected {} at {} but found {}: {}'.format(
                    chr_name, offset, found_chr, path))

            details['next_coord'] = parsed['start']


def select_next_sample_file(sample_files):
    next_sample = None
    for details in sample_files:
        if details['next_offset'] is None:
            continue

        if (((next_sample is None)
             or (details['next_coord'] < next_sample['next_coord']))):
            next_sample = details

    return next_sample


# pending_genes is sorted by end coordinate
def drop_old_genes(start, pending_genes):
    keep_i = None
    for gene_i, gene in enumerate(pending_genes['genes']):
        if gene['end'] < start:
            continue

        keep_i = gene_i
        break

    dropped_genes = list()
    if keep_i is None:
        dropped_genes = pending_genes['genes']
        pending_genes['genes'] = list()
    elif keep_i > 0:
        dropped_genes = pending_genes['genes'][:keep_i]
        pending_genes['genes'] = pending_genes['genes'][keep_i:]

    return dropped_genes


def load_next_gene(handle, gene_i):
    line = handle.readline()
    while line:
        columns = rmats_long_utils.read_tsv_line(line)
        # Gene lines have 4 columns, transcript lines have 6
        if len(columns) != 4:
            line = handle.readline()
            continue

        start_str = columns[0]
        start = int(start_str)
        end_str = columns[1]
        end = int(end_str)
        gene_id = columns[2]
        return {
            'start': start,
            'end': end,
            'gene_i': gene_i,
            'gene_id': gene_id
        }

    return None


def add_new_pending_genes(end, pending_genes, gtf_handle):
    new_genes = list()
    next_gene = pending_genes['next_gene']
    if (next_gene is not None) and (next_gene['start'] > end):
        return new_genes

    gene_i = 0
    if next_gene is not None:
        gene_i = next_gene['gene_i'] + 1

    while True:
        if next_gene is None:
            next_gene = load_next_gene(gtf_handle, gene_i)
            pending_genes['next_gene'] = next_gene
            if next_gene is None:
                break

            gene_i += 1

        if next_gene['start'] > end:
            break

        new_genes.append(next_gene)
        next_gene = None

    if new_genes:
        combined = pending_genes['genes'] + new_genes
        pending_genes['genes'] = sorted(combined, key=lambda x: x['end'])

    return new_genes


def check_for_gene_overlap(start, end, pending_genes, gtf_handle):
    result = dict()
    overlapped_genes = list()
    result['genes'] = overlapped_genes
    past_some_gene = False
    for gene in pending_genes['genes']:
        if start > gene['end']:
            past_some_gene = True
        elif end >= gene['start']:
            overlapped_genes.append(gene)

    result['past_some_gene'] = past_some_gene
    new_genes = add_new_pending_genes(end, pending_genes, gtf_handle)
    for gene in new_genes:
        # All new_genes have gene['start'] <= end
        if start <= gene['end']:
            overlapped_genes.append(gene)

    return result


def process_reads(chr_name, sample_file, sorted_samples, total_by_sample,
                  total_by_gene_by_sample, pending_genes, gtf_handle,
                  out_handle, gene_count_handle):
    sample = sample_file['sample']
    read_total = total_by_sample.get(sample, 0)
    with open(sample_file['path'], 'rt') as read_handle:
        offset = sample_file['next_offset']
        if offset != 0:
            read_handle.seek(offset)

        line = read_handle.readline()
        if line:
            parsed = parse_read_line(line)
            # Since this sample_file was selected based on start coordinate, this
            # start can be used to see which genes are no longer needed.
            dropped_genes = drop_old_genes(parsed['start'], pending_genes)
            write_gene_count_lines(sorted_samples, dropped_genes,
                                   total_by_gene_by_sample, gene_count_handle)

        while line:
            if parsed['chr'] != chr_name:
                sample_file['next_offset'] = None
                break

            read_total += 1
            overlap_result = check_for_gene_overlap(parsed['start'],
                                                    parsed['end'],
                                                    pending_genes, gtf_handle)
            # The first 2 columns are gene_i and gene_id which
            # are set for each gene.
            out_columns = [None, None, sample]
            # Add alignment info, but skip columns[0] which is chr_name
            out_columns += parsed['columns'][1:]
            for gene in overlap_result['genes']:
                gene_id = gene['gene_id']
                out_columns[0] = str(gene['gene_i'])
                out_columns[1] = gene_id
                gene_total_by_sample = rmats_long_utils.try_get_or_set_default(
                    total_by_gene_by_sample, gene_id, dict())
                old_gene_total = gene_total_by_sample.get(sample, 0)
                gene_total_by_sample[sample] = old_gene_total + 1
                rmats_long_utils.write_tsv_line(out_handle, out_columns)

            offset = read_handle.tell()
            line = read_handle.readline()
            if line:
                parsed = parse_read_line(line)

            # Move on to a new file. When the other files are also
            # past the gene then it can be removed from pending_genes.
            if overlap_result['past_some_gene']:
                break

        if line and parsed['chr'] == chr_name:
            sample_file['next_coord'] = parsed['start']
            sample_file['next_offset'] = offset
        else:
            sample_file['next_offset'] = None

    total_by_sample[sample] = read_total


def sort_by_gene_i(path, buffer_size):
    dir_path = os.path.dirname(os.path.abspath(path))
    with tempfile.TemporaryDirectory(suffix='_tmp',
                                     prefix='sort_by_gene',
                                     dir=dir_path) as temp_dir:
        tmp_path = os.path.join(temp_dir, 'sorted.tmp')
        gene_key_arg = '-k1,1g'
        env = {'LC_ALL': 'C'}  # to ensure sort order
        command = [
            'sort', '--buffer-size', buffer_size, '--temporary-directory',
            temp_dir, gene_key_arg, '--output', tmp_path, path
        ]
        subprocess.run(command, env=env, check=True)
        shutil.move(tmp_path, path)


def remove_gene_id_from_extra_lines_and_create_index(align_path):
    dir_path = os.path.dirname(os.path.abspath(align_path))
    index_path = '{}.index'.format(align_path)
    with tempfile.TemporaryDirectory(suffix='_tmp',
                                     prefix='remove_gene_id',
                                     dir=dir_path) as temp_dir:
        tmp_path = os.path.join(temp_dir, 'sorted.tmp')
        current_gene_i = None
        with open(index_path, 'wt') as index_handle:
            with open(tmp_path, 'wt') as out_handle:
                with open(align_path, 'rt') as in_handle:
                    for line in in_handle:
                        columns = rmats_long_utils.read_tsv_line(line)
                        gene_i = columns[0]
                        if (((current_gene_i is None)
                             or (gene_i != current_gene_i))):
                            current_gene_i = gene_i
                            gene_id = columns[1]
                            offset = out_handle.tell()
                            out_handle.write(line)
                            out_columns = [gene_id, str(offset)]
                            rmats_long_utils.write_tsv_line(
                                index_handle, out_columns)
                            continue

                        columns[1] = ''
                        rmats_long_utils.write_tsv_line(out_handle, columns)

        shutil.move(tmp_path, align_path)


def write_gene_count_lines(sorted_samples, dropped_genes,
                           total_by_gene_by_sample, gene_count_handle):
    for gene in dropped_genes:
        gene_i = gene['gene_i']
        gene_id = gene['gene_id']
        total_by_sample = total_by_gene_by_sample.pop(gene_id, dict())
        columns = [str(gene_i), gene_id]
        for sample in sorted_samples:
            total = total_by_sample.get(sample, 0)
            columns.append(str(total))

        rmats_long_utils.write_tsv_line(gene_count_handle, columns)


def update_total_by_sample(partial_by_sample, total_by_sample):
    for sample, value in partial_by_sample.items():
        old_total = total_by_sample.get(sample, 0)
        total_by_sample[sample] = old_total + value


def write_total_by_sample(sorted_samples, total_by_sample, out_dir):
    out_path = os.path.join(out_dir, 'sample_read_totals.tsv')
    headers = ['sample', 'total']
    with open(out_path, 'wt') as out_handle:
        rmats_long_utils.write_tsv_line(out_handle, headers)
        for sample in sorted_samples:
            total = total_by_sample[sample]
            columns = [sample, str(total)]
            rmats_long_utils.write_tsv_line(out_handle, columns)


def calculate_cpm_multiplier_by_sample(total_by_sample):
    cpm_mult_by_sample = dict()
    for sample, total in total_by_sample.items():
        if total == 0:
            mult = 0
        else:
            mult = 1e6 / total

        cpm_mult_by_sample[sample] = mult

    return cpm_mult_by_sample


def combine_gene_count_files_to_cpm(chr_to_id, sorted_samples, total_by_sample,
                                    gene_count_tmp_dir, out_dir):
    chr_ids = sorted(chr_to_id.values())
    cpm_mult_by_sample = calculate_cpm_multiplier_by_sample(total_by_sample)
    combined_path = os.path.join(out_dir, 'sample_gene_cpm.tsv')
    with open(combined_path, 'wt') as out_handle:
        headers = ['gene_id']
        headers.extend(sorted_samples)
        rmats_long_utils.write_tsv_line(out_handle, headers)

        for chr_id in chr_ids:
            gene_count_path = rmats_long_utils.get_chr_file_path(
                gene_count_tmp_dir, chr_id)
            if not os.path.exists(gene_count_path):
                continue

            with open(gene_count_path, 'rt') as count_handle:
                for line in count_handle:
                    columns = rmats_long_utils.read_tsv_line(line)
                    # gene_i = columns[0]
                    gene_id = columns[1]
                    cpm_columns = [gene_id]
                    for sample_i, sample in enumerate(sorted_samples):
                        sample_count = columns[2 + sample_i]
                        sample_count = int(sample_count)
                        cpm_mult = cpm_mult_by_sample.get(sample, 0)
                        sample_cpm = sample_count * cpm_mult
                        formatted = rmats_long_utils.format_float(sample_cpm)
                        cpm_columns.append(formatted)

                    rmats_long_utils.write_tsv_line(out_handle, cpm_columns)


def organize_alignments_thread(buffer_size, sample_files, sorted_samples,
                               input_queue, output_queue):
    total_by_sample = dict()
    while True:
        input_data = input_queue.get()
        if input_data is None:
            break

        chr_name = input_data['chr_name']
        offsets_for_chr = input_data['offsets_for_chr']
        out_path = input_data['out_path']
        gtf_path = input_data['gtf_path']
        gene_count_path = input_data['gene_count_path']

        total_by_gene_by_sample = dict()
        initialize_next_read_offsets(sample_files, chr_name, offsets_for_chr)
        pending_genes = {'genes': list(), 'next_gene': None}
        with open(out_path, 'wt') as out_handle:
            with open(gtf_path, 'rt') as gtf_handle:
                with open(gene_count_path, 'wt') as gene_count_handle:
                    while True:
                        sample_file = select_next_sample_file(sample_files)
                        if sample_file is None:
                            break

                        process_reads(chr_name, sample_file, sorted_samples,
                                      total_by_sample, total_by_gene_by_sample,
                                      pending_genes, gtf_handle, out_handle,
                                      gene_count_handle)

                    write_gene_count_lines(sorted_samples,
                                           pending_genes['genes'],
                                           total_by_gene_by_sample,
                                           gene_count_handle)

        sort_by_gene_i(out_path, buffer_size)
        sort_by_gene_i(gene_count_path, buffer_size)
        remove_gene_id_from_extra_lines_and_create_index(out_path)

    output_queue.put(total_by_sample)


def organize_alignments_and_write_files(gtf_dir, buffer_size, num_threads,
                                        sample_files, sorted_samples, out_dir,
                                        gene_count_tmp_dir):
    chr_id_mapping = rmats_long_utils.parse_chr_mapping(gtf_dir)
    chr_to_id = chr_id_mapping['chr_to_id']
    rmats_long_utils.copy_chr_name_mapping(gtf_dir, out_dir)

    indices_by_path = load_indices(sample_files)
    num_chr_names = len(chr_to_id)
    num_threads = min(num_threads, num_chr_names)
    input_queue = multiprocessing.Queue(num_chr_names + num_threads)
    output_queue = multiprocessing.Queue(num_threads)
    threads = list()
    for _ in range(num_threads):
        thread = multiprocessing.Process(target=organize_alignments_thread,
                                         args=(buffer_size, sample_files,
                                               sorted_samples, input_queue,
                                               output_queue))
        threads.append(thread)
        thread.start()

    for chr_name, chr_id in chr_to_id.items():
        offsets_for_chr = dict()
        for path, offsets in indices_by_path.items():
            offset = offsets.get(chr_name)
            offsets_for_chr[path] = offset

        out_path = rmats_long_utils.get_chr_file_path(out_dir, chr_id)
        gtf_path = rmats_long_utils.get_chr_file_path(gtf_dir, chr_id)
        gene_count_path = rmats_long_utils.get_chr_file_path(
            gene_count_tmp_dir, chr_id)
        input_queue.put({
            'chr_name': chr_name,
            'offsets_for_chr': offsets_for_chr,
            'out_path': out_path,
            'gtf_path': gtf_path,
            'gene_count_path': gene_count_path
        })

    # Signal the end of inputs
    for _ in range(num_threads):
        input_queue.put(None)

    rmats_long_utils.join_threads_and_raise_if_error(threads)

    total_by_sample = dict()
    for _ in range(num_threads):
        thread_total_by_sample = output_queue.get()
        update_total_by_sample(thread_total_by_sample, total_by_sample)

    write_total_by_sample(sorted_samples, total_by_sample, out_dir)
    combine_gene_count_files_to_cpm(chr_to_id, sorted_samples, total_by_sample,
                                    gene_count_tmp_dir, out_dir)


def remove_empty_out_files(out_dir):
    file_names = os.listdir(out_dir)
    for name in file_names:
        path = os.path.join(out_dir, name)
        size_in_bytes = rmats_long_utils.get_file_size(path)
        if size_in_bytes == 0:
            os.remove(path)


def main():
    args = parse_args()
    out_dir = os.path.abspath(args.out_dir)
    rmats_long_utils.create_output_dir(out_dir, check_empty=True)
    sample_files = parse_samples_tsv(args.samples_tsv)
    sorted_samples = get_sorted_samples(sample_files)
    with tempfile.TemporaryDirectory(suffix='_tmp',
                                     prefix='gene_count_',
                                     dir=out_dir) as gene_count_tmp_dir:
        organize_alignments_and_write_files(args.gtf_dir,
                                            args.sort_buffer_size,
                                            args.num_threads, sample_files,
                                            sorted_samples, out_dir,
                                            gene_count_tmp_dir)

    remove_empty_out_files(out_dir)
    print('organize_alignment_info_by_gene_and_chr.py finished')


if __name__ == '__main__':
    main()
