import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Count isoforms that meet cutoffs'))
    parser.add_argument('--diff-transcripts',
                        required=True,
                        help='The path to the differential transcript results')
    parser.add_argument(
        '--adj-pvalue',
        type=float,
        default=0.05,
        help='The cutoff for adjusted p-value (default: %(default)s)')
    parser.add_argument(
        '--delta-proportion',
        type=float,
        default=0.1,
        help='The cutoff for delta isoform proportion (default: %(default)s)')

    return parser.parse_args()


def count_significant_isoforms(transcript_path, adj_pvalue, delta_proportion):
    genes = set()
    isoform_count = 0
    pvalue_col_name = 'adj_pvalue'
    delta_prop_col_name = 'delta_isoform_proportion'
    with open(transcript_path, 'rt') as handle:
        for line_i, line in enumerate(handle):
            columns = line.rstrip('\n').split('\t')
            if line_i == 0:
                headers = columns
                continue

            row = dict(zip(headers, columns))
            gene = row['gene_id']
            found_pvalue_str = row[pvalue_col_name]
            found_prop_str = row[delta_prop_col_name]
            if (found_pvalue_str == 'NA') or (found_prop_str == 'NA'):
                continue

            found_pvalue = float(found_pvalue_str)
            found_prop = float(found_prop_str)
            if (((found_pvalue <= adj_pvalue)
                 and (abs(found_prop) >= delta_proportion))):
                genes.add(gene)
                isoform_count += 1

    gene_count = len(genes)
    print('found {} isoforms from {} genes with'
          ' {} <= {} and abs({}) >= {}'.format(isoform_count, gene_count,
                                               pvalue_col_name, adj_pvalue,
                                               delta_prop_col_name,
                                               delta_proportion))


def main():
    args = parse_args()
    count_significant_isoforms(args.diff_transcripts, args.adj_pvalue,
                               args.delta_proportion)


if __name__ == '__main__':
    main()
