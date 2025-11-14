import argparse
import os
import os.path
import shutil
import tempfile

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Detect differential isoform expression using a'
                     ' multinomial model'))
    parser.add_argument(
        '--abundance', help='The path to the abundance.esp file from ESPRESSO')
    parser.add_argument('--asm-counts',
                        help='The output file from count_reads_for_asms.py')
    parser.add_argument('--out-dir',
                        required=True,
                        help='The path to use as the output directory')
    parser.add_argument(
        '--group-1',
        required=True,
        help=('The path to a file listing the sample names for group 1. The'
              ' file should have a single line with the sample names as a'
              ' comma separated list. The sample names should match with the'
              ' ESPRESSO abundance column names.'))
    parser.add_argument(
        '--group-2',
        required=True,
        help='The path to a file listing the sample names for group 2.')
    parser.add_argument(
        '--num-threads',
        type=int,
        default=1,
        help='The number of threads to use (default: %(default)s)')
    parser.add_argument(
        '--adj-pvalue',
        type=float,
        default=0.05,
        help='The cutoff for adjusted p-value (default: %(default)s)')
    parser.add_argument('--use-unadjusted-pvalue',
                        action='store_true',
                        help='Use pvalue instead of adj_pvalue for the cutoff')
    parser.add_argument(
        '--delta-proportion',
        type=float,
        default=0.05,
        help='The cutoff for delta isoform proportion (default: %(default)s)')
    parser.add_argument(
        '--covar-tsv',
        help=('A .tsv with 1 line per sample. The first line has the column'
              ' names. The first column is sample_id. Each additional column'
              ' is a covariate.'))
    parser.add_argument('--use-drimseq',
                        action='store_true',
                        help='Use DRIMSeq instead of run_stat_model.py')
    parser.add_argument(
        '--min-cpm-per-asm',
        type=float,
        default=0,
        help=('Only consider ASMs where at least 1 sample has at least'
              ' this CPM of reads assigned to the ASM.'
              ' (default: %(default)s)'))
    parser.add_argument('--sample-read-total-tsv',
                        required=True,
                        help=('A .tsv file with two columns: sample and total.'
                              ' The 1st line is the header'))
    parser.add_argument(
        '--min-isoform-reads',
        type=int,
        default=1,
        help=('Only consider isoforms with at least this many reads'
              ' (default: %(default)s)'))
    parser.add_argument(
        '--limit-asm-to-top-n-isoforms',
        type=int,
        default=50,
        help=('Only consider the top N isoforms with the highest'
              ' total proportion across samples for each ASM'
              ' (default: %(default)s)'))
    parser.add_argument(
        '--average-reads-per-group',
        type=float,
        default=10,
        help=('For each sample group require the average read count to be at'
              ' least this value in order to be significant'
              ' (default: %(default)s)'))
    parser.add_argument(
        '--average-cpm-per-group',
        type=float,
        help=('For each sample group require the average CPM to be at least'
              ' this value in order to be significant'))
    parser.add_argument(
        '--min-cpm-per-group',
        type=float,
        help=('For each sample group require the min CPM to be at least'
              ' this value in order to be significant'))
    parser.add_argument(
        '--asm-proportion-of-gene',
        type=float,
        default=0.05,
        help=('Require the ASM CPM to be at least this proportion of the'
              ' gene CPM in at least 1 sample in order to be significant'
              ' (default: %(default)s)'))
    parser.add_argument(
        '--gene-cpm-tsv',
        required=True,
        help=('A .tsv file with gene_id as the 1st header and then an'
              ' additional header for each sample name'))

    args = parser.parse_args()
    has_abun = args.abundance is not None
    has_asm = args.asm_counts is not None
    if (has_abun and has_asm) or not (has_abun or has_asm):
        parser.error('Exactly 1 of --abundance or --asm-counts is required')
    if args.use_drimseq and has_asm:
        parser.error('DRIMSeq is not implemented for --asm-counts')

    return args


def sort_tsv_by_columns(tsv_path, column_names):
    rows_with_sort_key = list()
    with open(tsv_path, 'rt') as handle:
        for line_i, line in enumerate(handle):
            columns = rmats_long_utils.read_tsv_line(line)
            if line_i == 0:
                header_line = line
                header = columns
                column_indices = list()
                for name in column_names:
                    column_indices.append(header.index(name))

                continue

            sort_key = list()
            for column_i in column_indices:
                sort_key.append(columns[column_i])

            rows_with_sort_key.append((line, sort_key))

    rows_with_sort_key.sort(key=lambda p: p[1])
    with open(tsv_path, 'wt') as handle:
        handle.write(header_line)
        for row, key in rows_with_sort_key:
            handle.write(row)


def detect_isoforms(script_dir, abundance_path, asm_counts, is_asm,
                    use_drimseq, out_dir, num_threads, group_1_path,
                    group_2_path, covar_tsv, min_isoform_reads,
                    min_cpm_per_asm, limit_asm_to_top_n_isoforms,
                    sample_read_total_tsv):
    if use_drimseq:
        script_path = os.path.join(script_dir,
                                   'detect_differential_isoforms.R')
        command = [
            'Rscript', script_path, abundance_path, out_dir,
            str(num_threads), group_1_path, group_2_path
        ]
        if covar_tsv:
            command.append(covar_tsv)
    else:
        script_path = os.path.join(script_dir, 'run_stat_model.py')
        python_executable = rmats_long_utils.get_python_executable()
        command = [
            python_executable, script_path, '--out-dir', out_dir,
            '--num-threads',
            str(num_threads), '--group-1', group_1_path, '--group-2',
            group_2_path, '--min-isoform-reads',
            str(min_isoform_reads), '--min-cpm-per-asm',
            str(min_cpm_per_asm), '--sample-read-total-tsv',
            sample_read_total_tsv, '--limit-asm-to-top-n-isoforms',
            str(limit_asm_to_top_n_isoforms)
        ]
        if is_asm:
            command.extend(['--counts-tsv', asm_counts])
        else:
            command.extend(['--abundance', abundance_path])

        if covar_tsv:
            command.extend(['--covar-tsv', covar_tsv])

    rmats_long_utils.run_command(command)


def rename_and_change_headers(in_path, out_path, header_changes):
    with open(in_path, 'rt') as in_handle:
        with open(out_path, 'wt') as out_handle:
            columns_to_delete = list()
            for line_i, line in enumerate(in_handle):
                columns = rmats_long_utils.read_tsv_line(line)
                if line_i == 0:
                    new_headers = list()
                    for header_i, header in enumerate(columns):
                        new_header = header_changes.get(header, header)
                        if new_header is None:
                            columns_to_delete.append(header_i)
                        else:
                            new_headers.append(new_header)

                    columns_to_delete.reverse()
                    rmats_long_utils.write_tsv_line(out_handle, new_headers)
                    continue

                # columns_to_delete is in descending order so that
                # the elements are removed starting from the back of the list.
                for col_i in columns_to_delete:
                    columns.pop(col_i)

                rmats_long_utils.write_tsv_line(out_handle, columns)

    os.remove(in_path)


def sort_output(out_dir, is_asm, use_drimseq):
    if is_asm:
        return  # run_stat_model.py already sorted the files

    gene_tsv = os.path.join(out_dir, 'differential_genes.tsv')
    transcript_tsv = os.path.join(out_dir, 'differential_transcripts.tsv')
    if not use_drimseq:
        # Change file names and column names
        asm_tsv = os.path.join(out_dir, 'differential_asms.tsv')
        isoform_tsv = os.path.join(out_dir, 'differential_isoforms.tsv')
        rename_and_change_headers(asm_tsv, gene_tsv, {'asm_id': None})
        rename_and_change_headers(isoform_tsv, transcript_tsv, {
            'asm_id': None,
            'isoform_id': 'feature_id'
        })

    sort_tsv_by_columns(gene_tsv, ['gene_id'])
    sort_tsv_by_columns(transcript_tsv, ['gene_id', 'feature_id'])


def calculate_isoform_proportion(script_dir, out_dir, abundance_path, is_asm,
                                 group_1_path, group_2_path,
                                 sample_read_total_tsv, python_executable):
    proportion_script_path = os.path.join(script_dir,
                                          'calculate_isoform_proportion.py')

    if is_asm:
        count_path = os.path.join(out_dir, 'count.tsv')
        diff_transcripts = os.path.join(out_dir, 'differential_isoforms.tsv')
        diff_asms = os.path.join(out_dir, 'differential_asms.tsv')
    else:
        diff_transcripts = os.path.join(out_dir,
                                        'differential_transcripts.tsv')
        diff_asms = os.path.join(out_dir, 'differential_genes.tsv')

    tmp_file = '{}.tmp'.format(diff_transcripts)
    command = [
        python_executable, proportion_script_path, '--diff-transcripts',
        diff_transcripts, '--group-1', group_1_path, '--group-2', group_2_path,
        '--tmp-file', tmp_file, '--sample-read-total-tsv',
        sample_read_total_tsv
    ]
    if is_asm:
        command.extend(['--asm-counts', count_path])
    else:
        command.extend(['--abundance', abundance_path])

    rmats_long_utils.run_command(command)

    return {'transcripts': diff_transcripts, 'asms': diff_asms}


class GeneCpmReader:
    def __init__(self, handle):
        self.handle = handle
        self.headers = None
        self.samples = list()
        self.current_gene = None
        self.current_cpm_by_sample = None

    def read_header_line(self):
        if not self.handle:
            return

        line = self.handle.readline()
        self.headers = rmats_long_utils.read_tsv_line(line)
        for header_i, header in enumerate(self.headers):
            if header_i == 0:
                expected_header = 'gene_id'
                if header != expected_header:
                    raise Exception(
                        'Expected {} as the first column but found {}'.format(
                            expected_header, header))

                continue

            self.samples.append(header)

    def get_sample_names(self):
        return self.samples[:]

    def get_cpm_for_gene(self, gene):
        if gene == self.current_gene:
            return self.current_cpm_by_sample

        if self.headers is None:
            self.read_header_line()

        while True:
            line = self.handle.readline()
            if not line:
                return None

            columns = rmats_long_utils.read_tsv_line(line)
            row = dict(zip(self.headers, columns))
            found_gene = row['gene_id']
            if found_gene == gene:
                cpm_by_sample = dict()
                for sample in self.samples:
                    cpm = float(row[sample])
                    cpm_by_sample[sample] = cpm

                self.current_gene = gene
                self.current_cpm_by_sample = cpm_by_sample
                return cpm_by_sample


def get_sample_names_from_transcript_headers(transcript_headers):
    sample_names = list()
    cpm_suffix = '_cpm'
    for header in transcript_headers:
        if header.endswith(cpm_suffix):
            name = header[:-len(cpm_suffix)]
            sample_names.append(name)

    return sample_names


def write_asm_header_line(sample_names, orig_headers, asms_out_handle):
    new_cpm_headers = list()
    new_prop_headers = list()
    by_sample = dict()
    for sample in sample_names:
        cpm_header = '{}_cpm'.format(sample)
        prop_header = '{}_proportion_of_gene_cpm'.format(sample)
        new_cpm_headers.append(cpm_header)
        new_prop_headers.append(prop_header)
        by_sample[sample] = {'cpm': cpm_header, 'prop': prop_header}

    new_headers = new_cpm_headers + new_prop_headers
    out_columns = orig_headers[:]
    out_columns.extend(new_headers)
    rmats_long_utils.write_tsv_line(asms_out_handle, out_columns)
    return {'by_sample': by_sample, 'all': new_headers}


def write_asm_row(sample_names, orig_asm_headers, new_headers, asm_row,
                  gene_cpm_by_sample, asm_cpm_by_sample, out_handle):
    headers_by_sample = new_headers['by_sample']
    all_new_headers = new_headers['all']
    for sample in sample_names:
        sample_headers = headers_by_sample[sample]
        cpm_header = sample_headers['cpm']
        prop_header = sample_headers['prop']
        gene_cpm = gene_cpm_by_sample[sample]
        asm_cpm = asm_cpm_by_sample[sample]
        asm_row[cpm_header] = asm_cpm
        if not gene_cpm:
            prop_of_gene = 0
        else:
            prop_of_gene = asm_cpm / gene_cpm

        asm_row[prop_header] = prop_of_gene

    out_columns = list()
    for header in orig_asm_headers:
        value = asm_row[header]
        out_columns.append(value)

    for header in all_new_headers:
        value = asm_row[header]
        formatted = rmats_long_utils.format_float(value)
        out_columns.append(formatted)

    rmats_long_utils.write_tsv_line(out_handle, out_columns)


def add_asm_cpm_and_gene_ratio_columns_with_handles(transcripts_handle,
                                                    asms_in_handle,
                                                    gene_cpm_reader, is_asm,
                                                    asms_out_handle):
    asm_header_line = asms_in_handle.readline()
    if not asm_header_line:
        raise Exception('No headers in ASM tsv file')

    orig_asm_headers = rmats_long_utils.read_tsv_line(asm_header_line)
    if is_asm:
        asm_id_column = 'asm_id'
    else:
        # If starting from isoform abundance, treat each gene as an ASM
        asm_id_column = 'gene_id'

    transcript_header_line = transcripts_handle.readline()
    if not transcript_header_line:
        raise Exception('No headers in transcript tsv file')

    transcript_headers = rmats_long_utils.read_tsv_line(transcript_header_line)
    sample_names = get_sample_names_from_transcript_headers(transcript_headers)
    new_asm_headers = write_asm_header_line(sample_names, orig_asm_headers,
                                            asms_out_handle)

    asm_line = asms_in_handle.readline()
    transcript_line = None
    while asm_line:
        orig_asm_columns = rmats_long_utils.read_tsv_line(asm_line)
        asm_row = dict(zip(orig_asm_headers, orig_asm_columns))
        asm_id = asm_row[asm_id_column]
        gene_id = asm_row['gene_id']
        gene_cpm_by_sample = gene_cpm_reader.get_cpm_for_gene(gene_id)
        asm_cpm_by_sample = dict()
        for sample in sample_names:
            asm_cpm_by_sample[sample] = 0

        # After the first ASM, the transcript line will already be set
        # since the transcript loop breaks when a new ASM is found.
        if transcript_line is None:
            transcript_line = transcripts_handle.readline()

        while transcript_line:
            transcript_columns = rmats_long_utils.read_tsv_line(
                transcript_line)
            transcript_row = dict(zip(transcript_headers, transcript_columns))
            found_id = transcript_row[asm_id_column]
            if found_id != asm_id:
                write_asm_row(sample_names, orig_asm_headers, new_asm_headers,
                              asm_row, gene_cpm_by_sample, asm_cpm_by_sample,
                              asms_out_handle)
                break

            for sample in sample_names:
                sample_cpm = rmats_long_utils.get_sample_cpm_from_row(
                    sample, transcript_row)
                asm_cpm_by_sample[sample] += sample_cpm

            transcript_line = transcripts_handle.readline()

        # If the the end of the transcript file is reached then
        # write out the new ASM line.
        if not transcript_line:
            write_asm_row(sample_names, orig_asm_headers, new_asm_headers,
                          asm_row, gene_cpm_by_sample, asm_cpm_by_sample,
                          asms_out_handle)

        asm_line = asms_in_handle.readline()


def add_asm_cpm_and_gene_ratio_columns(diff_transcripts_path, diff_asms_path,
                                       gene_cpm_tsv, is_asm):
    asms_basename = os.path.basename(diff_asms_path)
    out_dir = os.path.dirname(diff_asms_path)
    with tempfile.TemporaryDirectory(
            suffix='_tmp',
            prefix='add_asm_cpm_and_gene_ratio_columns_',
            dir=out_dir) as temp_dir:
        tmp_asms_path = os.path.join(temp_dir, asms_basename)
        with open(diff_transcripts_path, 'rt') as transcripts_handle:
            with open(diff_asms_path, 'rt') as asms_in_handle:
                with open(gene_cpm_tsv, 'rt') as gene_cpm_handle:
                    with open(tmp_asms_path, 'wt') as asms_out_handle:
                        gene_cpm_reader = GeneCpmReader(gene_cpm_handle)
                        gene_cpm_reader.read_header_line()
                        add_asm_cpm_and_gene_ratio_columns_with_handles(
                            transcripts_handle, asms_in_handle,
                            gene_cpm_reader, is_asm, asms_out_handle)

        shutil.move(tmp_asms_path, diff_asms_path)


def count_significant_isoforms(script_dir, out_dir, python_executable,
                               diff_transcripts, diff_asms, adj_pvalue,
                               use_unadjusted_pvalue, delta_proportion,
                               average_reads_per_group, average_cpm_per_group,
                               min_cpm_per_group, asm_proportion_of_gene,
                               group_1_path, group_2_path, is_asm):
    count_script_path = os.path.join(script_dir,
                                     'count_significant_isoforms.py')
    if is_asm:
        filtered_transcripts = os.path.join(
            out_dir, 'differential_isoforms_filtered.tsv')
    else:
        filtered_transcripts = os.path.join(
            out_dir, 'differential_transcripts_filtered.tsv')

    command = [
        python_executable, count_script_path, '--diff-transcripts',
        diff_transcripts, '--diff-asms', diff_asms, '--out-tsv',
        filtered_transcripts, '--adj-pvalue',
        str(adj_pvalue), '--delta-proportion',
        str(delta_proportion), '--average-reads-per-group',
        str(average_reads_per_group), '--group-1', group_1_path, '--group-2',
        group_2_path, '--asm-proportion-of-gene',
        str(asm_proportion_of_gene)
    ]
    if is_asm:
        command.append('--is-asm')

    if use_unadjusted_pvalue:
        command.append('--use-unadjusted-pvalue')

    if average_cpm_per_group is not None:
        command.extend(['--average-cpm-per-group', str(average_cpm_per_group)])

    if min_cpm_per_group is not None:
        command.extend(['--min-cpm-per-group', str(min_cpm_per_group)])

    rmats_long_utils.run_command(command)


def detect_differential_isoforms(args):
    script_dir = rmats_long_utils.get_script_dir()
    python_executable = rmats_long_utils.get_python_executable()

    is_asm = args.asm_counts is not None
    detect_isoforms(script_dir, args.abundance, args.asm_counts, is_asm,
                    args.use_drimseq, args.out_dir, args.num_threads,
                    args.group_1, args.group_2, args.covar_tsv,
                    args.min_isoform_reads, args.min_cpm_per_asm,
                    args.limit_asm_to_top_n_isoforms,
                    args.sample_read_total_tsv)
    sort_output(args.out_dir, is_asm, args.use_drimseq)
    diff_paths = calculate_isoform_proportion(script_dir, args.out_dir,
                                              args.abundance, is_asm,
                                              args.group_1, args.group_2,
                                              args.sample_read_total_tsv,
                                              python_executable)
    diff_transcripts_path = diff_paths['transcripts']
    diff_asms_path = diff_paths['asms']

    add_asm_cpm_and_gene_ratio_columns(diff_transcripts_path, diff_asms_path,
                                       args.gene_cpm_tsv, is_asm)
    count_significant_isoforms(
        script_dir, args.out_dir, python_executable, diff_transcripts_path,
        diff_asms_path, args.adj_pvalue, args.use_unadjusted_pvalue,
        args.delta_proportion, args.average_reads_per_group,
        args.average_cpm_per_group, args.min_cpm_per_group,
        args.asm_proportion_of_gene, args.group_1, args.group_2, is_asm)


def main():
    args = parse_args()
    detect_differential_isoforms(args)


if __name__ == '__main__':
    main()
