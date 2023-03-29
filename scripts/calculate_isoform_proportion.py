import argparse
import shutil


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Calculate isoform proportion for each sample'))
    parser.add_argument('--diff-transcripts',
                        required=True,
                        help=('The path to the DRIMSeq transcript output file.'
                              ' Additional columns will be appended.'))
    parser.add_argument('--abundance',
                        required=True,
                        help='The path to abundance.esp file from ESPRESSO')
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

    return parser.parse_args()


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


def parse_abundance(abundance_path):
    details = dict()
    total_by_gene_by_sample = dict()
    counts_by_isoform_by_sample = dict()
    details = {
        'total_by_gene_by_sample': total_by_gene_by_sample,
        'counts_by_isoform_by_sample': counts_by_isoform_by_sample,
    }
    with open(abundance_path, 'rt') as handle:
        for line_i, line in enumerate(handle):
            columns = line.rstrip('\n').split('\t')
            if line_i == 0:
                headers = columns
                samples = headers[3:]
                continue

            row = dict(zip(headers, columns))
            gene = row['gene_ID']
            isoform = row['transcript_ID']
            if gene == 'NA':
                continue

            total_by_sample = total_by_gene_by_sample.get(gene)
            if not total_by_sample:
                total_by_sample = dict()
                total_by_gene_by_sample[gene] = total_by_sample

            count_by_sample = counts_by_isoform_by_sample.get(isoform)
            if not count_by_sample:
                count_by_sample = dict()
                counts_by_isoform_by_sample[isoform] = count_by_sample

            for sample in samples:
                count = float(row[sample])
                count_by_sample[sample] = count
                old_total = total_by_sample.get(sample, 0)
                total_by_sample[sample] = old_total + count

    return details


def write_tsv_line(handle, columns):
    handle.write('{}\n'.format('\t'.join(columns)))


def calc_group_proportions(samples, total_by_sample, counts_by_sample):
    group_props = list()
    for sample in samples:
        gene_count = total_by_sample[sample]
        isoform_count = counts_by_sample[sample]
        if gene_count == 0:
            prop = 'NA'
        else:
            prop = isoform_count / gene_count

        group_props.append(prop)

    return group_props


def calc_average(values):
    if 'NA' in values:
        return 'NA'

    total = sum(values)
    return total / len(values)


def calc_delta(a, b):
    if 'NA' in [a, b]:
        return 'NA'

    return a - b


def append_proportion_columns(diff_transcripts_path, abundance_details,
                              group_1_samples, group_2_samples, tmp_path):
    total_by_gene_by_sample = abundance_details['total_by_gene_by_sample']
    counts_by_isoform_by_sample = (
        abundance_details['counts_by_isoform_by_sample'])
    with open(tmp_path, 'wt') as out_handle:
        with open(diff_transcripts_path, 'rt') as in_handle:
            for line_i, line in enumerate(in_handle):
                columns = line.rstrip('\n').split('\t')
                if line_i == 0:
                    orig_headers = columns
                    new_headers = orig_headers[:]
                    for sample in group_1_samples + group_2_samples:
                        new_headers.append('{}_proportion'.format(sample))

                    new_headers.append('group_1_average_proportion')
                    new_headers.append('group_2_average_proportion')
                    new_headers.append('delta_isoform_proportion')
                    write_tsv_line(out_handle, new_headers)
                    continue

                row = dict(zip(orig_headers, columns))
                gene = row['gene_id']
                isoform = row['feature_id']
                total_by_sample = total_by_gene_by_sample[gene]
                counts_by_sample = counts_by_isoform_by_sample[isoform]
                group_1_props = calc_group_proportions(group_1_samples,
                                                       total_by_sample,
                                                       counts_by_sample)
                group_2_props = calc_group_proportions(group_2_samples,
                                                       total_by_sample,
                                                       counts_by_sample)
                group_1_avg = calc_average(group_1_props)
                group_2_avg = calc_average(group_2_props)
                delta = calc_delta(group_1_avg, group_2_avg)
                new_values = (group_1_props + group_2_props +
                              [group_1_avg, group_2_avg, delta])
                for value in new_values:
                    if value == 'NA':
                        columns.append(value)
                    else:
                        columns.append(str(round(value, 4)))

                write_tsv_line(out_handle, columns)

    shutil.move(tmp_path, diff_transcripts_path)


def calculate_isoform_proportion(diff_transcripts_path, abundance_path,
                                 group_1_path, group_2_path, tmp_path):
    group_1_samples = parse_group_file(group_1_path)
    group_2_samples = parse_group_file(group_2_path)
    abundance_details = parse_abundance(abundance_path)
    append_proportion_columns(diff_transcripts_path, abundance_details,
                              group_1_samples, group_2_samples, tmp_path)


def main():
    args = parse_args()
    calculate_isoform_proportion(args.diff_transcripts, args.abundance,
                                 args.group_1, args.group_2, args.tmp_file)


if __name__ == '__main__':
    main()
