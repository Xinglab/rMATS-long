import argparse
import math

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Count isoforms that meet the cutoff values'))
    parser.add_argument('--diff-transcripts',
                        required=True,
                        help='The path to the differential transcript results')
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
        default=0.1,
        help='The cutoff for delta isoform proportion (default: %(default)s)')

    return parser.parse_args()


def count_significant_isoforms(transcript_path, out_path, adj_pvalue,
                               use_unadjusted_pvalue, delta_proportion):
    with open(transcript_path, 'rt') as in_handle:
        with open(out_path, 'wt') as out_handle:
            count_significant_isoforms_with_handles(in_handle, out_handle,
                                                    adj_pvalue,
                                                    use_unadjusted_pvalue,
                                                    delta_proportion)


def count_significant_isoforms_with_handles(in_handle, out_handle, adj_pvalue,
                                            use_unadjusted_pvalue,
                                            delta_proportion):
    genes = set()
    isoform_count = 0
    pvalue_col_name = 'adj_pvalue'
    if use_unadjusted_pvalue:
        pvalue_col_name = 'pvalue'

    delta_prop_col_name = 'delta_isoform_proportion'
    for line_i, line in enumerate(in_handle):
        columns = line.rstrip('\n').split('\t')
        if line_i == 0:
            headers = columns
            out_handle.write(line)
            continue

        row = dict(zip(headers, columns))
        gene = row['gene_id']
        found_pvalue = rmats_long_utils.parse_float(row[pvalue_col_name])
        found_prop = rmats_long_utils.parse_float(row[delta_prop_col_name])
        if math.isnan(found_pvalue) or math.isnan(found_prop):
            continue

        if (((found_pvalue <= adj_pvalue)
             and (abs(found_prop) >= delta_proportion))):
            genes.add(gene)
            isoform_count += 1
            out_handle.write(line)

    gene_count = len(genes)
    print('found {} isoforms from {} genes with'
          ' {} <= {} and abs({}) >= {}'.format(isoform_count, gene_count,
                                               pvalue_col_name, adj_pvalue,
                                               delta_prop_col_name,
                                               delta_proportion))


def main():
    args = parse_args()
    count_significant_isoforms(args.diff_transcripts, args.out_tsv,
                               args.adj_pvalue, args.use_unadjusted_pvalue,
                               args.delta_proportion)


if __name__ == '__main__':
    main()
