import datetime
import importlib
import importlib.resources
import math
import os
import os.path
import queue
import shutil
import subprocess
import sys


def get_rmats_long_executable():
    return 'rmats-long'


def get_python_executable():
    python_executable = sys.executable
    if not python_executable:
        python_executable = 'python'

    return python_executable


def get_script_command_line():
    py_exe = get_python_executable()
    command_parts = [py_exe]
    command_parts.extend(sys.argv)
    return ' '.join(command_parts)


def get_script_dir():
    py_script_rel_path = sys.argv[0]
    py_script_abs_path = os.path.abspath(py_script_rel_path)
    script_dir = os.path.dirname(py_script_abs_path)
    return script_dir


def get_source_code_commit():
    file_name = 'source_code_commit.txt'
    try:
        with importlib.resources.path('rmats_long', file_name) as path_object:
            path = os.fspath(path_object)
            with open(path, 'rt') as handle:
                commit_string = handle.read()
                return commit_string.strip()
    except:
        return 'unknown'


def create_output_dir(dir_path, check_empty=False):
    if os.path.exists(dir_path):
        if not os.path.isdir(dir_path):
            raise Exception('{} already exists and is a file'.format(dir_path))

        if check_empty:
            files = os.listdir(dir_path)
            if len(files) != 0:
                raise Exception(
                    '{} exists and is not an empty directory'.format(dir_path))

        return

    os.makedirs(dir_path)


def get_chr_id_mapping_path(dir_path):
    return os.path.join(dir_path, 'chr_name_id_mapping.tsv')


def get_chr_file_path(dir_path, chr_id):
    return os.path.join(dir_path, 'chr_id_{}.tsv'.format(chr_id))


def get_graph_file_path(dir_path, chr_id):
    return os.path.join(dir_path, 'graph_{}.txt'.format(chr_id))


def all_sjs_prefix():
    return 'all_sjs_'


def get_all_sjs_path(dir_path, chr_id):
    prefix = all_sjs_prefix()
    return os.path.join(dir_path, '{}{}.tsv'.format(prefix, chr_id))


def find_gene_graph_in_graph_file(gene_id, path):
    found_gene = False
    graph_lines = list()
    if not (path and os.path.exists(path)):
        return graph_lines

    with open(path, 'rt') as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not found_gene:
                if 'digraph' not in line:
                    continue

                line_gene = line.split()[2]
                if line_gene != gene_id:
                    continue

                found_gene = True
                graph_lines.append(line)
                continue

            graph_lines.append(line)
            if line == '}':
                break

    return graph_lines


def get_chr_id_from_path(path):
    # name is like chr_id_{}.tsv
    prefix = 'chr_id_'
    suffix = '.tsv'
    return get_chr_id_from_path_with_prefix_and_suffix(path, prefix, suffix)


def get_chr_id_from_all_sjs_path(path):
    prefix = all_sjs_prefix()
    suffix = '.tsv'
    return get_chr_id_from_path_with_prefix_and_suffix(path, prefix, suffix)


def get_chr_id_from_path_with_prefix_and_suffix(path, prefix, suffix):
    name = os.path.basename(path)
    if not (name.startswith(prefix) and name.endswith(suffix)):
        raise Exception('Unexpected format for {}'.format(path))

    chr_id = name[len(prefix):-len(suffix)]
    return chr_id


def copy_chr_name_mapping(orig_dir, new_dir):
    orig_path = get_chr_id_mapping_path(orig_dir)
    new_path = get_chr_id_mapping_path(new_dir)
    shutil.copy(orig_path, new_path)


def parse_chr_mapping(dir_path):
    chr_to_id = dict()
    id_to_chr = dict()
    path = get_chr_id_mapping_path(dir_path)
    with open(path, 'rt') as handle:
        for row in row_iterator_for_tsv_with_header(handle):
            chr_name = row['name']
            chr_id = row['id']
            chr_id = int(chr_id)
            chr_to_id[chr_name] = chr_id
            id_to_chr[chr_id] = chr_name

    return {'chr_to_id': chr_to_id, 'id_to_chr': id_to_chr}


def parse_asm_id(asm_id):
    parts = asm_id.split('_')
    if len(parts) != 2:
        raise Exception('could not parse asm_id: {}'.format(asm_id))

    chr_id = parts[0]
    event_i = int(parts[1])
    return {'chr': chr_id, 'event_i': event_i}


def parse_asm_isoform_id(isoform_id):
    parts = isoform_id.split('_')
    if len(parts) != 3:
        raise Exception(
            'could not parse ASM isoform_id: {}'.format(isoform_id))

    asm_id = '_'.join(parts[:2])
    isoform_i = int(parts[2])
    return {'asm': asm_id, 'isoform_i': isoform_i}


def create_sjs_string(sjs):
    strings = list()
    for start, end in sjs:
        sj = '{}-{}'.format(start, end)
        strings.append(sj)

    return ','.join(strings)


def parse_sjs_string(sjs_string):
    sjs = list()
    if not sjs_string:
        return sjs

    parts = sjs_string.split(',')
    for part in parts:
        start_str, end_str = part.split('-')
        start = int(start_str)
        end = int(end_str)
        sjs.append((start, end))

    return sjs


def format_isoforms_str_from_exons(isoforms):
    isoform_strs = list()
    for isoform in isoforms:
        coord_strs = list()
        for exon in isoform:
            start, end = exon
            coord_strs.append(str(start))
            coord_strs.append(str(end))

        isoform_str = ','.join(coord_strs)
        isoform_strs.append(isoform_str)

    return ';'.join(isoform_strs)


def format_isoforms_str(isoforms):
    isoform_strs = list()
    for isoform in isoforms:
        coord_strs = list()
        for coord in isoform:
            coord_strs.append(str(coord))

        isoform_str = ','.join(coord_strs)
        isoform_strs.append(isoform_str)

    return ';'.join(isoform_strs)


def parse_isoforms_str(isoforms_str):
    isoforms = list()
    parts = isoforms_str.split(';')
    for part in parts:
        coord_strings = part.split(',')
        isoform = list()
        for coord_str in coord_strings:
            coord = int(coord_str)
            isoform.append(coord)

        isoforms.append(isoform)

    return isoforms


def format_isoform_ids_str(isoform_id_strs):
    return ';'.join(isoform_id_strs)


def parse_isoform_ids_str(isoform_ids_str):
    return isoform_ids_str.split(';')


def parse_group_file(file_path):
    groups = None
    with open(file_path, 'rt') as handle:
        for line_i, line in enumerate(handle):
            if line and line_i != 0:
                print('expected only 1 line in {}'.format(file_path))

            groups = line.rstrip('\n').split(',')

    if not groups:
        raise Exception('could not parse groups from: {}'.format(file_path))

    return groups


def parse_gtf_line(line):
    if line.strip().startswith('#'):
        # skip comment lines
        return None

    columns = read_tsv_line(line)
    chr_name = columns[0]
    source = columns[1]
    feature = columns[2]
    start_str = columns[3]
    start = int(start_str)
    end_str = columns[4]
    end = int(end_str)
    score = columns[5]
    strand = columns[6]
    frame = columns[7]
    attributes_str = columns[8]
    attributes = parse_gtf_attributes(attributes_str)
    return {
        'chr': chr_name,
        'source': source,
        'feature': feature,
        'start': start,
        'end': end,
        'score': score,
        'strand': strand,
        'frame': frame,
        'attributes': attributes,
        'line': line
    }


def write_gtf_line(handle,
                   chr_name,
                   feature,
                   start,
                   end,
                   strand,
                   attributes,
                   source='.',
                   score='.',
                   frame='.'):
    columns = [
        chr_name, source, feature,
        str(start),
        str(end), score, strand, frame
    ]
    formatted_attrs = format_gtf_attributes(attributes)
    columns.append(formatted_attrs)
    write_tsv_line(handle, columns)


def format_gtf_attributes(attributes):
    key_value_strings = list()
    for key, value in attributes.items():
        key_value_strings.append('{} "{}";'.format(key, value))

    return ' '.join(key_value_strings)


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


def select_significant_transcripts(gene_id, diff_transcripts_path):
    if not diff_transcripts_path:
        return list()

    most_significant_positive = None
    most_significant_negative = None
    with open(diff_transcripts_path, 'rt') as transcripts_handle:
        for row in row_iterator_for_tsv_with_header(transcripts_handle):
            if row['gene_id'] != gene_id:
                continue

            transcript = row['feature_id']
            pvalue = parse_float(row['adj_pvalue'])
            delta = parse_float(row['delta_isoform_proportion'])
            if math.isnan(pvalue) or math.isnan(delta):
                continue

            parsed = {'id': transcript, 'pvalue': pvalue, 'delta': delta}
            if delta > 0:
                if ((most_significant_positive is None)
                        or (pvalue < most_significant_positive['pvalue'])):
                    most_significant_positive = parsed
            if delta < 0:
                if ((most_significant_negative is None)
                        or (pvalue < most_significant_negative['pvalue'])):
                    most_significant_negative = parsed

    if most_significant_positive and most_significant_negative:
        pos_pvalue = most_significant_positive['pvalue']
        neg_pvalue = most_significant_negative['pvalue']
        if pos_pvalue < neg_pvalue:
            return [
                most_significant_positive['id'],
                most_significant_negative['id']
            ]

        return [
            most_significant_negative['id'], most_significant_positive['id']
        ]
    if most_significant_positive:
        return [most_significant_positive['id']]
    if most_significant_negative:
        return [most_significant_negative['id']]

    return list()


def row_iterator_for_tsv_with_header(in_handle):
    for line_i, line in enumerate(in_handle):
        columns = read_tsv_line(line)
        if line_i == 0:
            headers = columns
            continue

        yield dict(zip(headers, columns))


def parse_abundance_file(abundance_path):
    counts_by_gene_by_transcript_by_sample = dict()
    total_by_gene_by_sample = dict()
    total_by_sample = dict()
    expected_headers = ['transcript_ID', 'transcript_name', 'gene_ID']
    with open(abundance_path, 'rt') as in_handle:
        for line_i, line in enumerate(in_handle):
            in_columns = read_tsv_line(line)
            if line_i == 0:
                in_headers = in_columns
                first_headers = in_headers[:len(expected_headers)]
                if first_headers != expected_headers:
                    raise Exception('expected headers in {} to start with {}'
                                    ' but found {}'.format(
                                        abundance_path, expected_headers,
                                        first_headers))

                sample_names = in_headers[len(expected_headers):]
                continue

            row = dict(zip(in_headers, in_columns))
            gene_str = row['gene_ID']
            genes = gene_str.split(',')
            transcript = row['transcript_ID']
            if genes == ['NA']:
                # Only update the total read count for each sample
                for sample in sample_names:
                    count = parse_float(row[sample])
                    old_total = total_by_sample.get(sample, 0)
                    total_by_sample[sample] = old_total + count

                continue

            for gene in genes:
                gene_total_by_sample = total_by_gene_by_sample.get(gene)
                if not gene_total_by_sample:
                    gene_total_by_sample = dict()
                    total_by_gene_by_sample[gene] = gene_total_by_sample

                counts_by_transcript_by_sample = counts_by_gene_by_transcript_by_sample.get(
                    gene)
                if not counts_by_transcript_by_sample:
                    counts_by_transcript_by_sample = dict()
                    counts_by_gene_by_transcript_by_sample[gene] = (
                        counts_by_transcript_by_sample)

                counts_by_sample = counts_by_transcript_by_sample.get(
                    transcript)
                if not counts_by_sample:
                    counts_by_sample = dict()
                    counts_by_transcript_by_sample[transcript] = (
                        counts_by_sample)

                for sample in sample_names:
                    count = parse_float(row[sample])
                    old_total = total_by_sample.get(sample, 0)
                    total_by_sample[sample] = old_total + count
                    old_gene_total = gene_total_by_sample.get(sample, 0)
                    gene_total_by_sample[sample] = old_gene_total + count
                    counts_by_sample[sample] = count

    return {
        'counts_by_gene_by_transcript_by_sample':
            counts_by_gene_by_transcript_by_sample,
        'total_by_gene_by_sample':
            total_by_gene_by_sample,
        'total_by_sample':
            total_by_sample,
        'sample_names':
            sample_names
    }


def parse_asm_counts_file(asm_counts):
    counts_by_asm_by_isoform_by_sample = dict()
    total_by_asm_by_sample = dict()
    total_by_sample = dict()
    expected_headers = [
        'asm_id', 'gene_id', 'sample_id', 'isoform_id', 'count', 'prop'
    ]
    with open(asm_counts, 'rt') as in_handle:
        for line_i, line in enumerate(in_handle):
            in_columns = read_tsv_line(line)
            if line_i == 0:
                in_headers = in_columns
                if in_headers != expected_headers:
                    raise Exception('expected headers in {} to be {}'
                                    ' but found {}'.format(
                                        asm_counts, expected_headers,
                                        in_headers))
                continue

            row = dict(zip(in_headers, in_columns))
            asm_id = row['asm_id']
            sample = row['sample_id']
            isoform_id = row['isoform_id']
            count = float(row['count'])
            asm_total_by_sample = try_get_or_set_default(
                total_by_asm_by_sample, asm_id, dict())
            counts_by_isoform_by_sample = try_get_or_set_default(
                counts_by_asm_by_isoform_by_sample, asm_id, dict())
            counts_by_sample = try_get_or_set_default(
                counts_by_isoform_by_sample, isoform_id, dict())
            old_total = total_by_sample.get(sample, 0)
            total_by_sample[sample] = old_total + count
            old_asm_total = asm_total_by_sample.get(sample, 0)
            asm_total_by_sample[sample] = old_asm_total + count
            counts_by_sample[sample] = count

    sample_names = sorted(total_by_sample.keys())
    return {
        'counts_by_asm_by_isoform_by_sample':
            counts_by_asm_by_isoform_by_sample,
        'total_by_asm_by_sample':
            total_by_asm_by_sample,
        'total_by_sample':
            total_by_sample,
        'sample_names':
            sample_names
    }


def parse_sample_totals_file(total_path):
    total_by_sample = dict()
    with open(total_path, 'rt') as handle:
        rows = row_iterator_for_tsv_with_header(handle)
        for row in rows:
            sample = row['sample']
            total = int(row['total'])
            total_by_sample[sample] = total

    return total_by_sample


def read_tsv_line(line):
    return line.rstrip('\n').split('\t')


def write_tsv_line(handle, columns):
    handle.write('{}\n'.format('\t'.join(columns)))


def parse_float(string):
    if string == 'NA':
        return math.nan

    return float(string)


def format_float(value, digits=4):
    if value is None or math.isnan(value):
        return 'NA'

    format_str = '{{:.{}g}}'.format(digits)
    formatted = format_str.format(value)
    return formatted


def round_float_string(string):
    value = parse_float(string)
    return format_float(value)


def parse_bool(string):
    if string == 'True':
        return True

    if string == 'False':
        return False

    raise Exception('Unexpected bool string: {}'.format(string))


def get_sample_cpm_from_row(sample, row):
    cpm_column = '{}_cpm'.format(sample)
    sample_cpm = row[cpm_column]
    sample_cpm = float(sample_cpm)
    return sample_cpm


def get_sample_count_from_row(sample, row):
    count_column = '{}_count'.format(sample)
    sample_count = row[count_column]
    sample_count = float(sample_count)
    return sample_count


def run_command(command):
    print('running: {}'.format(command))
    subprocess.run(command, check=True)


def run_package_r_script(script_name, arguments):
    # importlib.resources.path was deprecated but then un-deprecated:
    # https://github.com/python/cpython/issues/116608
    with importlib.resources.path('rmats_long', script_name) as path_object:
        script_path = os.fspath(path_object)
        abs_script_path = os.path.abspath(script_path)
        command = ['Rscript', abs_script_path]
        command.extend(arguments)
        run_command(command)


# import the module so that a compiled Cython version can be used.
def run_package_py_script(script_name, arguments):
    parts = script_name.split('.')
    if (len(parts) != 2) and parts[1] != 'py':
        raise Exception('expected {} to have the format {{script_name}}.py')

    module_name = parts[0]
    full_module_name = 'rmats_long.{}'.format(module_name)
    module = importlib.import_module(full_module_name)
    argv_backup = sys.argv
    new_argv = [script_name]
    new_argv.extend(arguments)
    sys.argv = new_argv
    module.main()
    sys.argv = argv_backup


def run_package_script(script_name, arguments):
    if script_name.endswith('.R'):
        run_package_r_script(script_name, arguments)
    else:
        run_package_py_script(script_name, arguments)


def try_get_or_set_default(container, key, default):
    value = container.get(key)
    if value is None:
        container[key] = default
        value = default

    return value


def try_get_from_queue_with_short_wait(in_queue):
    try:
        # Wait at most 1 second for value
        value = in_queue.get(True, 1)
    except queue.Empty:
        return None

    return value


def try_get_from_queue_without_wait(in_queue):
    try:
        value = in_queue.get(False)
    except queue.Empty:
        return None

    return value


def try_put_to_queue_with_short_wait(out_queue, value):
    try:
        out_queue.put(value, True, 1)
    except queue.Full:
        return False

    return True


def try_put_to_queue_without_wait(out_queue, value):
    try:
        out_queue.put(value, False)
    except queue.Full:
        return False

    return True


def drain_queue(to_drain):
    # This is a best effort to make sure that the queue is empty.
    # It's possible for the queue to raise queue.Empty due to locking
    # rather than actually being empty.
    while True:
        value = try_get_from_queue_with_short_wait(to_drain)
        if value is None:
            break


def raise_exception_if_thread_exited_early(threads):
    for thread in threads:
        if not thread.is_alive():
            exit_code = thread.exitcode
            raise Exception(
                'thread exited early with value: {}'.format(exit_code))


def try_join_thread(thread):
    thread.join(1)  # wait at most 1 second
    exit_code = thread.exitcode
    if exit_code is None:
        return 'running'
    if exit_code != 0:
        return 'error'

    return 'joined'


def join_threads_and_raise_if_error(threads):
    any_error = False
    while threads and (not any_error):
        running_threads = list()
        for thread in threads:
            join_result = try_join_thread(thread)
            if join_result == 'running':
                running_threads.append(thread)
                continue

            if join_result == 'joined':
                continue

            any_error = True

        threads = running_threads

    if any_error:
        for thread in threads:
            thread.terminate()

        message = 'A thread exited with an error.'
        if threads:
            message = '{} Terminated {} threads'.format(message, len(threads))

        raise Exception(message)


class HandleAndNextLine:
    def __init__(self, handle):
        self.handle = handle
        self.next_line = None
        self.offset = self.handle.tell()

    def read_line(self):
        if self.next_line is not None:
            return self.next_line

        self.next_line = self.handle.readline()
        return self.next_line

    def clear_line(self):
        self.offset = self.handle.tell()
        self.next_line = None

    def clear_and_read_line(self):
        self.clear_line()
        return self.read_line()


def get_file_size(path):
    stat_result = os.stat(path)
    size_in_bytes = stat_result.st_size
    return size_in_bytes


def get_current_time_string():
    format_str = '%Y-%m-%dT%H:%M:%S'
    now = datetime.datetime.now()
    time_str = now.strftime(format_str)
    return time_str
