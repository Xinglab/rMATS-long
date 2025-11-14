import argparse
import shutil

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Calculate isoform proportion for each sample')
    parser.add_argument(
        '--diff-transcripts',
        required=True,
        help=('The path to differential_isoforms.tsv from run_stat_model.py.'
              ' Additional columns will be appended.'))
    parser.add_argument('--abundance',
                        help='The path to abundance.esp file from ESPRESSO')
    parser.add_argument(
        '--asm-counts',
        help='The count.tsv output file from run_stat_model.py')
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
    parser.add_argument('--tmp-file',
                        required=True,
                        help='A path to use for temporary data')
    parser.add_argument('--sample-read-total-tsv',
                        required=True,
                        help=('A .tsv file with two columns: sample and total.'
                              ' The 1st line is the header'))

    args = parser.parse_args()
    has_abun = args.abundance is not None
    has_asm = args.asm_counts is not None
    if (has_abun and has_asm) or not (has_abun or has_asm):
        parser.error('Exactly 1 of --abundance or --asm-counts is required')

    return args


def parse_abundance(abundance_path):
    parsed_details = rmats_long_utils.parse_abundance_file(abundance_path)
    counts_by_gene_by_transcript_by_sample = (
        parsed_details['counts_by_gene_by_transcript_by_sample'])
    counts_by_isoform_by_sample = dict()
    for counts_by_transcript_by_sample in (
            counts_by_gene_by_transcript_by_sample.values()):
        for transcript, counts_by_sample in (
                counts_by_transcript_by_sample.items()):
            new_counts_by_sample = dict()
            counts_by_isoform_by_sample[transcript] = new_counts_by_sample
            for sample, count in counts_by_sample.items():
                new_counts_by_sample[sample] = count
    return {
        'total_by_gene_by_sample': parsed_details['total_by_gene_by_sample'],
        'counts_by_isoform_by_sample': counts_by_isoform_by_sample,
    }


def parse_asm_counts(asm_counts):
    parsed_details = rmats_long_utils.parse_asm_counts_file(asm_counts)
    counts_by_gene_by_transcript_by_sample = (
        parsed_details['counts_by_asm_by_isoform_by_sample'])
    counts_by_isoform_by_sample = dict()
    for counts_by_transcript_by_sample in (
            counts_by_gene_by_transcript_by_sample.values()):
        for transcript, counts_by_sample in (
                counts_by_transcript_by_sample.items()):
            new_counts_by_sample = dict()
            counts_by_isoform_by_sample[transcript] = new_counts_by_sample
            for sample, count in counts_by_sample.items():
                new_counts_by_sample[sample] = count
    return {
        'total_by_gene_by_sample': parsed_details['total_by_asm_by_sample'],
        'counts_by_isoform_by_sample': counts_by_isoform_by_sample,
    }


def calculate_group_values(samples, total_by_sample, counts_by_sample,
                           cpm_mult_by_sample):
    group_props = list()
    group_counts = list()
    group_cpm_values = list()
    for sample in samples:
        cpm_mult = cpm_mult_by_sample.get(sample, 0)
        gene_count = total_by_sample.get(sample, 0)
        isoform_count = counts_by_sample.get(sample, 0)
        if gene_count == 0:
            prop = 'NA'
        else:
            prop = isoform_count / gene_count

        group_props.append(prop)
        group_counts.append(isoform_count)
        cpm_value = isoform_count * cpm_mult
        group_cpm_values.append(cpm_value)

    return {
        'props': group_props,
        'counts': group_counts,
        'cpm_values': group_cpm_values
    }


def calc_average_without_na(values):
    total = 0
    count = 0
    for value in values:
        if value == 'NA':
            continue

        total += value
        count += 1

    if count == 0:
        return 'NA'

    return total / count


def calc_delta(a, b):
    if 'NA' in [a, b]:
        return 'NA'

    return a - b


def calculate_cpm_mult_by_sample(sample_read_total_tsv):
    total_by_sample = rmats_long_utils.parse_sample_totals_file(
        sample_read_total_tsv)
    cpm_mult_by_sample = dict()
    for sample, total in total_by_sample.items():
        cpm_mult = 1e6 / total
        cpm_mult_by_sample[sample] = cpm_mult

    return cpm_mult_by_sample


def add_sample_value_headers(all_samples, headers):
    sample_prop_headers = list()
    sample_count_headers = list()
    sample_cpm_headers = list()
    for sample in all_samples:
        sample_prop_headers.append('{}_proportion'.format(sample))
        sample_count_headers.append('{}_count'.format(sample))
        sample_cpm_headers.append('{}_cpm'.format(sample))

    headers.extend(sample_prop_headers)
    headers.append('group_1_average_proportion')
    headers.append('group_2_average_proportion')
    headers.append('delta_isoform_proportion')
    headers.extend(sample_count_headers)
    headers.extend(sample_cpm_headers)


def add_sample_values(group_1_props, group_2_props, group_1_avg, group_2_avg,
                      delta, group_1_values, group_2_values, columns):
    new_values = list()
    new_values.extend(group_1_props)
    new_values.extend(group_2_props)
    new_values.append(group_1_avg)
    new_values.append(group_2_avg)
    new_values.append(delta)
    new_values.extend(group_1_values['counts'])
    new_values.extend(group_2_values['counts'])
    new_values.extend(group_1_values['cpm_values'])
    new_values.extend(group_2_values['cpm_values'])
    for value in new_values:
        if value == 'NA':
            columns.append(value)
        else:
            columns.append(rmats_long_utils.format_float(value))


def append_proportion_columns(diff_transcripts_path, abundance_details, is_asm,
                              group_1_samples, group_2_samples,
                              sample_read_total_tsv, tmp_path):
    cpm_mult_by_sample = calculate_cpm_mult_by_sample(sample_read_total_tsv)
    total_by_gene_by_sample = abundance_details['total_by_gene_by_sample']
    counts_by_isoform_by_sample = (
        abundance_details['counts_by_isoform_by_sample'])
    all_samples = group_1_samples + group_2_samples
    with open(tmp_path, 'wt') as out_handle:
        with open(diff_transcripts_path, 'rt') as in_handle:
            for line_i, line in enumerate(in_handle):
                columns = rmats_long_utils.read_tsv_line(line)
                if line_i == 0:
                    orig_headers = columns
                    new_headers = orig_headers[:]
                    add_sample_value_headers(all_samples, new_headers)
                    rmats_long_utils.write_tsv_line(out_handle, new_headers)
                    continue

                row = dict(zip(orig_headers, columns))
                if is_asm:
                    gene = row['asm_id']
                    isoform = row['isoform_id']
                else:
                    gene = row['gene_id']
                    isoform = row['feature_id']

                total_by_sample = total_by_gene_by_sample[gene]
                counts_by_sample = counts_by_isoform_by_sample[isoform]
                group_1_values = calculate_group_values(
                    group_1_samples, total_by_sample, counts_by_sample,
                    cpm_mult_by_sample)
                group_2_values = calculate_group_values(
                    group_2_samples, total_by_sample, counts_by_sample,
                    cpm_mult_by_sample)
                group_1_props = group_1_values['props']
                group_2_props = group_2_values['props']
                group_1_avg = calc_average_without_na(group_1_props)
                group_2_avg = calc_average_without_na(group_2_props)
                delta = calc_delta(group_1_avg, group_2_avg)
                add_sample_values(group_1_props, group_2_props, group_1_avg,
                                  group_2_avg, delta, group_1_values,
                                  group_2_values, columns)
                rmats_long_utils.write_tsv_line(out_handle, columns)

    shutil.move(tmp_path, diff_transcripts_path)


def calculate_isoform_proportion(diff_transcripts_path, abundance_path,
                                 asm_counts, group_1_path, group_2_path,
                                 sample_read_total_tsv, tmp_path):
    group_1_samples = rmats_long_utils.parse_group_file(group_1_path)
    group_2_samples = rmats_long_utils.parse_group_file(group_2_path)
    if abundance_path:
        abundance_details = parse_abundance(abundance_path)
        is_asm = False
    else:
        abundance_details = parse_asm_counts(asm_counts)
        is_asm = True

    append_proportion_columns(diff_transcripts_path, abundance_details, is_asm,
                              group_1_samples, group_2_samples,
                              sample_read_total_tsv, tmp_path)


def main():
    args = parse_args()
    calculate_isoform_proportion(args.diff_transcripts, args.abundance,
                                 args.asm_counts, args.group_1, args.group_2,
                                 args.sample_read_total_tsv, args.tmp_file)


if __name__ == '__main__':
    main()
