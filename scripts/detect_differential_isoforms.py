import argparse
import os
import os.path
import subprocess
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Detect differential isoform expression using DRIMSeq'))
    parser.add_argument('--abundance',
                        required=True,
                        help='The path to abundance.esp file from ESPRESSO')
    parser.add_argument('--out-dir',
                        required=True,
                        help='The path to use as output directory')
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
        required=False,
        help='The number of threads to use (default: %(default)s)')

    return parser.parse_args()


def sort_tsv_by_columns(tsv_path, column_names):
    rows_with_sort_key = list()
    with open(tsv_path, 'rt') as handle:
        for line_i, line in enumerate(handle):
            columns = line.rstrip('\n').split('\t')
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


def run_command(command):
    print('running: {}'.format(command))
    subprocess.run(command, check=True)


def detect_isoforms(script_dir, abundance_path, out_dir, num_threads,
                    group_1_path, group_2_path):
    r_script_path = os.path.join(script_dir, 'detect_differential_isoforms.R')
    command = [
        'Rscript', r_script_path, abundance_path, out_dir,
        str(num_threads), group_1_path, group_2_path
    ]
    run_command(command)


def sort_output(out_dir):
    gene_tsv = os.path.join(out_dir, 'differential_genes.tsv')
    transcript_tsv = os.path.join(out_dir, 'differential_transcripts.tsv')
    sort_tsv_by_columns(gene_tsv, ['gene_id'])
    sort_tsv_by_columns(transcript_tsv, ['gene_id', 'feature_id'])


def calculate_isoform_proportion(script_dir, out_dir, abundance_path,
                                 group_1_path, group_2_path,
                                 python_executable):
    proportion_script_path = os.path.join(script_dir,
                                          'calculate_isoform_proportion.py')
    diff_transcripts = os.path.join(out_dir, 'differential_transcripts.tsv')
    tmp_file = '{}.tmp'.format(diff_transcripts)
    command = [
        python_executable, proportion_script_path, '--diff-transcripts',
        diff_transcripts, '--abundance', abundance_path, '--group-1',
        group_1_path, '--group-2', group_2_path, '--tmp-file', tmp_file
    ]
    run_command(command)


def count_significant_isoforms(script_dir, out_dir, python_executable):
    count_script_path = os.path.join(script_dir,
                                     'count_significant_isoforms.py')
    diff_transcripts = os.path.join(out_dir, 'differential_transcripts.tsv')
    command = [
        python_executable, count_script_path, '--diff-transcripts',
        diff_transcripts
    ]
    run_command(command)


def get_python_executable():
    python_executable = sys.executable
    if not python_executable:
        python_executable = 'python'

    return python_executable


def detect_differential_isoforms(args):
    py_script_rel_path = sys.argv[0]
    py_script_abs_path = os.path.abspath(py_script_rel_path)
    script_dir = os.path.dirname(py_script_abs_path)
    python_executable = get_python_executable()

    detect_isoforms(script_dir, args.abundance, args.out_dir, args.num_threads,
                    args.group_1, args.group_2)
    sort_output(args.out_dir)
    calculate_isoform_proportion(script_dir, args.out_dir, args.abundance,
                                 args.group_1, args.group_2, python_executable)
    count_significant_isoforms(script_dir, args.out_dir, python_executable)


def main():
    args = parse_args()
    detect_differential_isoforms(args)


if __name__ == '__main__':
    main()
