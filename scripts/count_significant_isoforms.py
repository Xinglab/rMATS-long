import argparse
import math

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Count isoforms that meet the cutoff values')
    parser.add_argument('--diff-transcripts',
                        required=True,
                        help='The path to the differential transcript results')
    parser.add_argument('--diff-asms',
                        required=True,
                        help='The path to the differential asms results')
    parser.add_argument(
        '--out-tsv',
        required=True,
        help='The path to write transcripts that meet the cutoff values')
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
    parser.add_argument('--is-asm',
                        action='store_true',
                        help='Use if running with ASM output')
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
        '--group-1',
        required=True,
        help=('The path to a file listing the sample names for group 1. The'
              ' file should have a single line with the sample names as a'
              ' comma separated list.'))
    parser.add_argument(
        '--group-2',
        required=True,
        help='The path to a file listing the sample names for group 2')

    return parser.parse_args()


def count_significant_isoforms(transcript_path, asm_path, out_path, adj_pvalue,
                               use_unadjusted_pvalue, delta_proportion,
                               average_reads_per_group, average_cpm_per_group,
                               min_cpm_per_group, asm_proportion_of_gene,
                               group_1_path, group_2_path, is_asm):
    group_1 = list()
    group_2 = list()
    if group_1_path:
        group_1 = rmats_long_utils.parse_group_file(group_1_path)

    if group_2_path:
        group_2 = rmats_long_utils.parse_group_file(group_2_path)

    with open(transcript_path, 'rt') as transcript_handle:
        with open(asm_path, 'rt') as asm_handle:
            with open(out_path, 'wt') as out_handle:
                count_significant_isoforms_with_handles(
                    transcript_handle, asm_handle, out_handle, adj_pvalue,
                    use_unadjusted_pvalue, delta_proportion,
                    average_reads_per_group, average_cpm_per_group,
                    min_cpm_per_group, asm_proportion_of_gene, group_1,
                    group_2, is_asm)


def check_if_any_sample_meets_gene_cpm_proportion(asm_proportion_of_gene,
                                                  asm_row):
    prop_suffix = '_proportion_of_gene_cpm'
    for name, value in asm_row.items():
        if not name.endswith(prop_suffix):
            continue

        prop = float(value)
        if prop >= asm_proportion_of_gene:
            return True

    return False


def check_if_meets_average_cpm_or_count_per_group(
        key, average_threshold, sample_group_cpm_and_count_info):
    group_1_by_sample = sample_group_cpm_and_count_info.get(1, dict())
    group_2_by_sample = sample_group_cpm_and_count_info.get(2, dict())
    num_group_1_samples = len(group_1_by_sample)
    num_group_2_samples = len(group_2_by_sample)
    if num_group_1_samples == 0:
        raise Exception('No samples in group 1')

    if num_group_2_samples == 0:
        raise Exception('No samples in group 2')

    group_1_sum = 0
    for cpm_and_count in group_1_by_sample.values():
        group_1_sum += cpm_and_count[key]

    group_2_sum = 0
    for cpm_and_count in group_2_by_sample.values():
        group_2_sum += cpm_and_count[key]

    group_1_avg = group_1_sum / num_group_1_samples
    group_2_avg = group_2_sum / num_group_2_samples
    return ((group_1_avg >= average_threshold)
            and (group_2_avg >= average_threshold))


def check_if_meets_average_reads_per_group(average_reads_per_group,
                                           sample_group_cpm_and_count_info):
    return check_if_meets_average_cpm_or_count_per_group(
        'count', average_reads_per_group, sample_group_cpm_and_count_info)


def check_if_meets_average_cpm_per_group(average_cpm_per_group,
                                         sample_group_cpm_and_count_info):
    return check_if_meets_average_cpm_or_count_per_group(
        'cpm', average_cpm_per_group, sample_group_cpm_and_count_info)


def check_if_meets_min_cpm_per_group(min_cpm_per_group,
                                     sample_group_cpm_and_count_info):
    group_1_by_sample = sample_group_cpm_and_count_info.get(1, dict())
    group_2_by_sample = sample_group_cpm_and_count_info.get(2, dict())
    group_1_min_cpm = None
    for cpm_and_count in group_1_by_sample.values():
        cpm = cpm_and_count['cpm']
        if group_1_min_cpm is None:
            group_1_min_cpm = cpm
        else:
            group_1_min_cpm = min(group_1_min_cpm, cpm)

    group_2_min_cpm = None
    for cpm_and_count in group_2_by_sample.values():
        cpm = cpm_and_count['cpm']
        if group_2_min_cpm is None:
            group_2_min_cpm = cpm
        else:
            group_2_min_cpm = min(group_2_min_cpm, cpm)

    return ((group_1_min_cpm >= min_cpm_per_group)
            and (group_2_min_cpm >= min_cpm_per_group))


def get_sample_cpm_and_count_info_by_group(rows, group_1, group_2):
    details = {1: dict(), 2: dict()}

    group_and_sample = list()
    for sample in group_1:
        group_and_sample.append((1, sample))

    for sample in group_2:
        group_and_sample.append((2, sample))

    for row in rows:
        for group, sample in group_and_sample:
            sample_cpm = rmats_long_utils.get_sample_cpm_from_row(sample, row)
            sample_count = rmats_long_utils.get_sample_count_from_row(
                sample, row)
            by_sample = details[group]
            cpm_and_count = by_sample.get(sample)
            if not cpm_and_count:
                cpm_and_count = {'cpm': 0, 'count': 0}
                by_sample[sample] = cpm_and_count

            cpm_and_count['cpm'] += sample_cpm
            cpm_and_count['count'] += sample_count

    return details


def process_asm_rows(rows, asm_row, asm, gene, pvalue_col_name,
                     delta_prop_col_name, adj_pvalue, delta_proportion,
                     average_reads_per_group, average_cpm_per_group,
                     min_cpm_per_group, asm_proportion_of_gene, group_1,
                     group_2, asms, genes, out_handle):
    isoform_count = 0
    sample_group_cpm_and_count_info = get_sample_cpm_and_count_info_by_group(
        rows, group_1, group_2)
    if asm_proportion_of_gene:
        any_had_cpm_proportion = check_if_any_sample_meets_gene_cpm_proportion(
            asm_proportion_of_gene, asm_row)
        if not any_had_cpm_proportion:
            return isoform_count

    if average_reads_per_group:
        meets_avg_reads = check_if_meets_average_reads_per_group(
            average_reads_per_group, sample_group_cpm_and_count_info)
        if not meets_avg_reads:
            return isoform_count

    if average_cpm_per_group:
        meets_avg_cpm = check_if_meets_average_cpm_per_group(
            average_cpm_per_group, sample_group_cpm_and_count_info)
        if not meets_avg_cpm:
            return isoform_count

    if min_cpm_per_group:
        meets_min_cpm = check_if_meets_min_cpm_per_group(
            min_cpm_per_group, sample_group_cpm_and_count_info)
        if not meets_min_cpm:
            return isoform_count

    for row in rows:
        found_pvalue = rmats_long_utils.parse_float(row[pvalue_col_name])
        found_prop = rmats_long_utils.parse_float(row[delta_prop_col_name])
        if math.isnan(found_pvalue) or math.isnan(found_prop):
            continue

        if (((found_pvalue <= adj_pvalue)
             and (abs(found_prop) >= delta_proportion))):
            asms.add(asm)
            genes.add(gene)
            isoform_count += 1
            out_handle.write(row['original_line'])

    return isoform_count


def print_message(isoform_count, asms, genes, pvalue_col_name, adj_pvalue,
                  delta_prop_col_name, delta_proportion,
                  average_reads_per_group, average_cpm_per_group,
                  min_cpm_per_group, asm_proportion_of_gene, is_asm):
    message_parts = list()
    message_parts.append('found {} isoforms'.format(isoform_count))
    if is_asm:
        asm_count = len(asms)
        message_parts.append('from {} ASMs'.format(asm_count))

    gene_count = len(genes)
    message_parts.append('from {} genes'.format(gene_count))
    message_parts.append('with')
    message_parts.append('{} <= {}'.format(pvalue_col_name, adj_pvalue))
    message_parts.append('and')
    message_parts.append('abs({}) >= {}'.format(delta_prop_col_name,
                                                delta_proportion))

    if average_reads_per_group:
        message_parts.append('and')
        message_parts.append(
            'average reads per group >= {}'.format(average_reads_per_group))

    if average_cpm_per_group:
        message_parts.append('and')
        message_parts.append(
            'average CPM per group >= {}'.format(average_cpm_per_group))

    if min_cpm_per_group:
        message_parts.append('and')
        message_parts.append(
            'min CPM per group >= {}'.format(min_cpm_per_group))

    if asm_proportion_of_gene:
        message_parts.append('and')
        message_parts.append('ASM CPM >= {}% of gene CPM'.format(
            asm_proportion_of_gene * 100))

    print(' '.join(message_parts))


def count_significant_isoforms_with_handles(
        transcript_handle, asm_handle, out_handle, adj_pvalue,
        use_unadjusted_pvalue, delta_proportion, average_reads_per_group,
        average_cpm_per_group, min_cpm_per_group, asm_proportion_of_gene,
        group_1, group_2, is_asm):
    asms = set()
    genes = set()
    isoform_count = 0
    rows_for_asm = list()
    current_asm = None
    current_gene = None
    pvalue_col_name = 'adj_pvalue'
    if use_unadjusted_pvalue:
        pvalue_col_name = 'pvalue'

    delta_prop_col_name = 'delta_isoform_proportion'

    asm_iterator = rmats_long_utils.row_iterator_for_tsv_with_header(
        asm_handle)
    asm_row = None
    for line_i, line in enumerate(transcript_handle):
        columns = rmats_long_utils.read_tsv_line(line)
        if line_i == 0:
            headers = columns
            out_handle.write(line)
            continue

        row = dict(zip(headers, columns))
        row['original_line'] = line
        gene = row['gene_id']
        if is_asm:
            asm = row['asm_id']
            if asm != current_asm:
                if rows_for_asm:
                    isoform_count += process_asm_rows(
                        rows_for_asm, asm_row, current_asm, current_gene,
                        pvalue_col_name, delta_prop_col_name, adj_pvalue,
                        delta_proportion, average_reads_per_group,
                        average_cpm_per_group, min_cpm_per_group,
                        asm_proportion_of_gene, group_1, group_2, asms, genes,
                        out_handle)

                rows_for_asm = list()
                current_asm = asm
                asm_row = next(asm_iterator, None)
                if (asm_row is None) or (asm_row['asm_id'] != asm):
                    raise Exception(
                        'matching ASM row not found for {}'.format(row))
        else:
            if gene != current_gene:
                if rows_for_asm:
                    isoform_count += process_asm_rows(
                        rows_for_asm, asm_row, current_asm, current_gene,
                        pvalue_col_name, delta_prop_col_name, adj_pvalue,
                        delta_proportion, average_reads_per_group,
                        average_cpm_per_group, min_cpm_per_group,
                        asm_proportion_of_gene, group_1, group_2, asms, genes,
                        out_handle)

                rows_for_asm = list()
                asm_row = next(asm_iterator, None)
                if (asm_row is None) or (asm_row['gene_id'] != gene):
                    raise Exception(
                        'matching gene row not found for {}'.format(row))

        current_gene = gene
        rows_for_asm.append(row)

    if rows_for_asm:
        isoform_count += process_asm_rows(
            rows_for_asm, asm_row, current_asm, current_gene, pvalue_col_name,
            delta_prop_col_name, adj_pvalue, delta_proportion,
            average_reads_per_group, average_cpm_per_group, min_cpm_per_group,
            asm_proportion_of_gene, group_1, group_2, asms, genes, out_handle)

    print_message(isoform_count, asms, genes, pvalue_col_name, adj_pvalue,
                  delta_prop_col_name, delta_proportion,
                  average_reads_per_group, average_cpm_per_group,
                  min_cpm_per_group, asm_proportion_of_gene, is_asm)


def main():
    args = parse_args()
    count_significant_isoforms(
        args.diff_transcripts, args.diff_asms, args.out_tsv, args.adj_pvalue,
        args.use_unadjusted_pvalue, args.delta_proportion,
        args.average_reads_per_group, args.average_cpm_per_group,
        args.min_cpm_per_group, args.asm_proportion_of_gene, args.group_1,
        args.group_2, args.is_asm)


if __name__ == '__main__':
    main()
