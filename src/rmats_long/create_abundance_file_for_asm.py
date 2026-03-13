import argparse
import os
import os.path

from rmats_long import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Create an abundance.esp file based on ASM files')
    parser.add_argument('--asm-id',
                        required=True,
                        help='The ID of the ASM to create output for')
    parser.add_argument(
        '--event-dir',
        required=True,
        help='The output directory from detect_splicing_events.py')
    parser.add_argument(
        '--sample-totals',
        required=True,
        help='The path to a .tsv file with read totals for each sample')
    parser.add_argument(
        '--count-tsv',
        required=True,
        help=('The path to a .tsv file with read counts by sample for'
              ' each isoform'))
    parser.add_argument('--out-file',
                        required=True,
                        help='The file for the output abundance')

    args = parser.parse_args()
    return args


def get_isoforms_for_asm(asm_id, event_dir):
    isoform_ids = list()
    asm_id_parts = asm_id.split('_')
    chr_i = asm_id_parts[0]
    tsv_name = 'chr_id_{}.tsv'.format(chr_i)
    tsv_path = os.path.join(event_dir, tsv_name)
    with open(tsv_path, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            found_asm = row['asm_id']
            if found_asm != asm_id:
                continue

            id_strings = row['isoform_ids']
            isoform_ids = id_strings.split(';')
            break

    return isoform_ids


def read_sample_totals(sample_totals_path):
    total_by_sample = dict()
    with open(sample_totals_path, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            sample = row['sample']
            total_str = row['total']
            total = float(total_str)
            total_by_sample[sample] = total

    return total_by_sample


def read_isoform_counts(isoform_ids, count_tsv):
    counts_by_isoform_by_sample = dict()
    for isoform_id in isoform_ids:
        counts_by_isoform_by_sample[isoform_id] = dict()

    with open(count_tsv, 'rt') as in_handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(
                in_handle):
            sample = row['sample_id']
            isoform = row['isoform_id']
            count = float(row['count'])
            counts_by_sample = counts_by_isoform_by_sample.get(isoform)
            if counts_by_sample is None:
                continue

            counts_by_sample[sample] = count

    return counts_by_isoform_by_sample


def create_abundance_file_for_asm(asm_id, event_dir, sample_totals_path,
                                  count_tsv, out_file):
    isoform_ids = get_isoforms_for_asm(asm_id, event_dir)
    total_by_sample = read_sample_totals(sample_totals_path)
    counts_by_isoform_by_sample = read_isoform_counts(isoform_ids, count_tsv)

    samples = sorted(total_by_sample)
    headers = ['transcript_ID', 'transcript_name', 'gene_ID']
    headers.extend(samples)
    out_columns = [None, 'NA', asm_id]
    remaining_for_samples = list()
    for sample in samples:
        sample_total = total_by_sample[sample]
        out_columns.append('0')
        remaining_for_samples.append(sample_total)

    with open(out_file, 'wt') as out_handle:
        rmats_long_utils.write_tsv_line(out_handle, headers)
        for isoform, counts_by_sample in counts_by_isoform_by_sample.items():
            out_columns[0] = isoform
            for sample_i, sample in enumerate(samples):
                value = counts_by_sample.get(sample, 0)
                formatted = rmats_long_utils.format_float(value, digits=2)
                out_columns[sample_i + 3] = formatted
                remaining_for_samples[sample_i] -= value

            rmats_long_utils.write_tsv_line(out_handle, out_columns)

        other_abundance_columns = ['NA', 'NA', 'NA']
        for remaining in remaining_for_samples:
            other_abundance_columns.append(str(remaining))

        rmats_long_utils.write_tsv_line(out_handle, other_abundance_columns)


def main():
    args = parse_args()
    create_abundance_file_for_asm(args.asm_id, args.event_dir,
                                  args.sample_totals, args.count_tsv,
                                  args.out_file)
    print('create_abundance_file_for_asm.py finished')


if __name__ == '__main__':
    main()
