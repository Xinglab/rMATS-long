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


def sort_tsv_by_column(tsv_path, column_name):
    rows_with_sort_key = list()
    na_rows = list()
    with open(tsv_path, 'rt') as handle:
        for line_i, line in enumerate(handle):
            columns = line.rstrip('\n').split('\t')
            if line_i == 0:
                header_line = line
                header = columns
                column_i = header.index(column_name)
                continue

            sort_string = columns[column_i]
            if sort_string == 'NA':
                na_rows.append(line)
            else:
                sort_value = float(sort_string)
                rows_with_sort_key.append((line, sort_value))

    rows_with_sort_key.sort(key=lambda p: p[1])
    with open(tsv_path, 'wt') as handle:
        handle.write(header_line)
        for row, key in rows_with_sort_key:
            handle.write(row)

        for row in na_rows:
            handle.write(row)


def detect_differential_isoforms(args):
    py_script_rel_path = sys.argv[0]
    py_script_abs_path = os.path.abspath(py_script_rel_path)
    script_dir = os.path.dirname(py_script_abs_path)
    r_script_path = os.path.join(script_dir, 'detect_differential_isoforms.R')

    command = [
        'Rscript', r_script_path, args.abundance, args.out_dir,
        str(args.num_threads), args.group_1, args.group_2
    ]
    print('running: {}'.format(command))
    subprocess.run(command, check=True)
    transcript_tsv = os.path.join(args.out_dir, 'differential_transcripts.tsv')
    gene_tsv = os.path.join(args.out_dir, 'differential_genes.tsv')
    sort_tsv_by_column(transcript_tsv, 'adj_pvalue')
    sort_tsv_by_column(gene_tsv, 'adj_pvalue')


def main():
    args = parse_args()
    detect_differential_isoforms(args)


if __name__ == '__main__':
    main()
