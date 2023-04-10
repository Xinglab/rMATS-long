import argparse
import math
import os
import os.path
import tempfile

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(description=(
        'Analyze ESPRESSO results and produce plots for significant isoforms'))
    parser.add_argument(
        '--abundance',
        required=True,
        help='The path to the abundance.esp file from ESPRESSO')
    parser.add_argument('--updated-gtf',
                        required=True,
                        help='The path to the updated.gtf file from ESPRESSO')
    parser.add_argument(
        '--gencode-gtf',
        help=('The path to a gencode annotation.gtf file. Will be used to'
              ' identify the Ensemble canonical isoform and the gene name'))
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
        help='The path to a file listing the sample names for group 2')
    parser.add_argument('--group-1-name',
                        default='group 1',
                        help='A name for group 1 (default %(default)s)')
    parser.add_argument('--group-2-name',
                        default='group 2',
                        help='A name for group 2 (default %(default)s)')
    parser.add_argument('--out-dir',
                        required=True,
                        help='The path to use as the output directory')
    parser.add_argument(
        '--num-threads',
        type=int,
        default=1,
        help='The number of threads to use (default %(default)s)')
    parser.add_argument(
        '--plot-top-n',
        type=int,
        default=10,
        help=('Generate plots for the top "n" significant genes. To plot all'
              ' significant genes use --plot-top-n -1. (default %(default)s)'))
    parser.add_argument(
        '--plot-file-type',
        choices=['.pdf', '.png'],
        default='.pdf',
        help='The file type for output plots (default %(default)s))')
    parser.add_argument(
        '--diff-transcripts',
        help=('The path to the differential transcript results. If given then'
              ' skip the differential isoform calculation.'))
    parser.add_argument(
        '--adj-pvalue',
        type=float,
        default=0.05,
        help='The cutoff for adjusted p-value (default %(default)s)')
    parser.add_argument(
        '--delta-proportion',
        type=float,
        default=0.1,
        help='The cutoff for delta isoform proportion (default %(default)s)')

    return parser.parse_args()


def process_genes(genes, temp_dir, sorted_paths, other_gene_abundance_totals,
                  out_dir, group_1, group_1_name, group_2, group_2_name,
                  plot_file_type):
    with open(sorted_paths['abundance'], 'rt') as abundance_handle:
        with open(sorted_paths['updated_gtf'], 'rt') as updated_gtf_handle:
            with open(sorted_paths['diff_transcripts'],
                      'rt') as diff_transcripts_handle:
                sorted_gencode_path = sorted_paths.get('gencode_gtf')
                if sorted_gencode_path:
                    with open(sorted_gencode_path, 'rt') as gencode_gtf_handle:
                        process_genes_with_handles(
                            genes, temp_dir, abundance_handle,
                            updated_gtf_handle, diff_transcripts_handle,
                            gencode_gtf_handle, other_gene_abundance_totals,
                            out_dir, group_1, group_1_name, group_2,
                            group_2_name, plot_file_type)
                else:
                    process_genes_with_handles(
                        genes, temp_dir, abundance_handle, updated_gtf_handle,
                        diff_transcripts_handle, None,
                        other_gene_abundance_totals, out_dir, group_1,
                        group_1_name, group_2, group_2_name, plot_file_type)


def process_genes_with_handles(genes, temp_dir, abundance_handle,
                               updated_gtf_handle, diff_transcripts_handle,
                               gencode_gtf_handle, other_gene_abundance_totals,
                               out_dir, group_1, group_1_name, group_2,
                               group_2_name, plot_file_type):
    abundance_info = {
        'handle': abundance_handle,
        'next_line': None,
        'header': None
    }
    updated_gtf_info = {
        'handle': updated_gtf_handle,
        'next_line': None,
        'header': None
    }
    diff_transcripts_info = {
        'handle': diff_transcripts_handle,
        'next_line': None,
        'header': None
    }
    gencode_gtf_info = {
        'handle': gencode_gtf_handle,
        'next_line': None,
        'header': None
    }
    set_header_line_if_possible(abundance_info)
    set_header_line_if_possible(diff_transcripts_info)
    for gene in genes:
        temp_files_for_gene = write_temp_files_for_gene(
            gene, temp_dir, abundance_info, updated_gtf_info,
            diff_transcripts_info, gencode_gtf_info,
            other_gene_abundance_totals)
        process_gene(gene, temp_files_for_gene['abundance'],
                     temp_files_for_gene['updated_gtf'],
                     temp_files_for_gene['gencode_gtf'],
                     temp_files_for_gene['diff_transcripts'], out_dir, group_1,
                     group_1_name, group_2, group_2_name, plot_file_type)


def set_header_line_if_possible(info):
    if not info['handle']:
        return

    first_line = info['handle'].readline()
    orig_gene, orig_line = parse_sorted_gene_line(first_line)
    if orig_gene is not None:
        raise Exception('Error looking for header in {}'.format(info))

    info['header'] = orig_line


def write_temp_files_for_gene(gene, temp_dir, abundance_info, updated_gtf_info,
                              diff_transcripts_info, gencode_gtf_info,
                              other_gene_abundance_totals):
    temp_abundance_path = os.path.join(temp_dir,
                                       '{}_abundance.esp'.format(gene))
    write_temp_file_for_gene(gene, temp_abundance_path, abundance_info)
    add_other_gene_totals_to_abundance(other_gene_abundance_totals[gene],
                                       temp_abundance_path)
    temp_updated_gtf_path = os.path.join(temp_dir,
                                         '{}_updated.gtf'.format(gene))
    write_temp_file_for_gene(gene, temp_updated_gtf_path, updated_gtf_info)
    temp_gencode_gtf_path = os.path.join(temp_dir,
                                         '{}_gencode.gtf'.format(gene))
    write_temp_file_for_gene(gene, temp_gencode_gtf_path, gencode_gtf_info)
    temp_diff_transcripts_path = os.path.join(
        temp_dir, '{}_diff_transcripts.tsv'.format(gene))
    write_temp_file_for_gene(gene, temp_diff_transcripts_path,
                             diff_transcripts_info)
    temp_files = {
        'abundance': temp_abundance_path,
        'updated_gtf': temp_updated_gtf_path,
        'gencode_gtf': None,
        'diff_transcripts': temp_diff_transcripts_path
    }
    if gencode_gtf_info['handle']:
        temp_files['gencode_gtf'] = temp_gencode_gtf_path

    return temp_files


def write_temp_file_for_gene(gene, temp_path, in_handle_info):
    in_handle = in_handle_info['handle']
    if in_handle is None:
        return

    with open(temp_path, 'wt') as out_handle:
        if in_handle_info['header']:
            out_handle.write(in_handle_info['header'])

        if in_handle_info['next_line'] is not None:
            next_line = in_handle_info['next_line']
            orig_gene, orig_line = parse_sorted_gene_line(next_line)
            if orig_gene is None:
                in_handle_info['next_line'] = None
            elif orig_gene == gene:
                out_handle.write(orig_line)
                in_handle_info['next_line'] = None
            else:
                return

        for next_line in in_handle:
            orig_gene, orig_line = parse_sorted_gene_line(next_line)
            if orig_gene is None:
                continue
            if orig_gene == gene:
                out_handle.write(orig_line)
            else:
                in_handle_info['next_line'] = next_line
                return


def add_other_gene_totals_to_abundance(other_gene_abundance_totals, out_path):
    with open(out_path, 'at') as out_handle:
        rmats_long_utils.write_tsv_line(out_handle,
                                        other_gene_abundance_totals)


def process_gene(gene_id, abundance, updated_gtf, gencode_gtf,
                 diff_transcripts, out_dir, group_1, group_1_name, group_2,
                 group_2_name, plot_file_type):
    visualize_isoforms(gene_id, abundance, updated_gtf, gencode_gtf,
                       diff_transcripts, out_dir, group_1, group_1_name,
                       group_2, group_2_name, plot_file_type)

    significant_transcripts = rmats_long_utils.select_significant_transcripts(
        gene_id, diff_transcripts)
    most_significant = significant_transcripts[0]
    classify_isoforms(gene_id, out_dir, updated_gtf, gencode_gtf,
                      most_significant)
    if gencode_gtf:
        gene_name, canonical_transcript = (
            rmats_long_utils.get_gene_name_and_canonical_transcript_from_gtf(
                gene_id, gencode_gtf))
        if canonical_transcript:
            classify_isoforms(gene_id, out_dir, updated_gtf, gencode_gtf,
                              canonical_transcript)


def count_significant_isoforms(out_dir, diff_transcripts, adj_pvalue,
                               delta_proportion):
    python_executable = rmats_long_utils.get_python_executable()
    script_dir = rmats_long_utils.get_script_dir()
    count_script = os.path.join(script_dir, 'count_significant_isoforms.py')
    out_tsv = os.path.join(out_dir, 'differential_transcripts_filtered.tsv')
    command = [
        python_executable, count_script, '--diff-transcripts',
        diff_transcripts, '--out-tsv', out_tsv, '--adj-pvalue',
        str(adj_pvalue), '--delta-proportion',
        str(delta_proportion)
    ]
    rmats_long_utils.run_command(command)


def detect_differential_isoforms(abundance, out_dir, group_1, group_2,
                                 adj_pvalue, delta_proportion, num_threads):
    python_executable = rmats_long_utils.get_python_executable()
    script_dir = rmats_long_utils.get_script_dir()
    detect_script = os.path.join(script_dir, 'detect_differential_isoforms.py')
    command = [
        python_executable, detect_script, '--abundance', abundance,
        '--out-dir', out_dir, '--group-1', group_1, '--group-2', group_2,
        '--adj-pvalue',
        str(adj_pvalue), '--delta-proportion',
        str(delta_proportion), '--num-threads',
        str(num_threads)
    ]
    rmats_long_utils.run_command(command)


def visualize_isoforms(gene_id, abundance, updated_gtf, gencode_gtf,
                       diff_transcripts, out_dir, group_1, group_1_name,
                       group_2, group_2_name, plot_file_type):
    python_executable = rmats_long_utils.get_python_executable()
    script_dir = rmats_long_utils.get_script_dir()
    visualize_script = os.path.join(script_dir, 'visualize_isoforms.py')
    command = [
        python_executable, visualize_script, '--gene-id', gene_id,
        '--abundance', abundance, '--updated-gtf', updated_gtf,
        '--diff-transcripts', diff_transcripts, '--out-dir', out_dir,
        '--group-1', group_1, '--group-2', group_2, '--group-1-name',
        group_1_name, '--group-2-name', group_2_name, '--plot-file-type',
        plot_file_type
    ]
    if gencode_gtf:
        command.extend(['--gencode-gtf', gencode_gtf])

    rmats_long_utils.run_command(command)


def classify_isoforms(gene, out_dir, updated_gtf, gencode_gtf,
                      main_transcript_id):
    python_executable = rmats_long_utils.get_python_executable()
    script_dir = rmats_long_utils.get_script_dir()
    classify_script = os.path.join(script_dir,
                                   'classify_isoform_differences.py')
    out_tsv = os.path.join(
        out_dir,
        '{}_isoform_differences_from_{}.tsv'.format(gene, main_transcript_id))
    command = [
        python_executable, classify_script, '--updated-gtf', updated_gtf,
        '--out-tsv', out_tsv, '--main-transcript-id', main_transcript_id
    ]
    if gencode_gtf:
        command.extend(['--gencode-gtf', gencode_gtf])

    rmats_long_utils.run_command(command)


def get_top_n_genes(diff_transcripts, plot_top_n):
    genes_to_sort = dict()
    with open(diff_transcripts, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            pvalue = rmats_long_utils.parse_float(row['adj_pvalue'])
            if math.isnan(pvalue):
                pvalue = 2  # to be last in sort

            abs_delta = abs(
                rmats_long_utils.parse_float(row['delta_isoform_proportion']))
            if math.isnan(abs_delta):
                abs_delta = -1  # to be last in sort

            gene = row['gene_id']
            sort_tuple = (pvalue, -abs_delta, gene)
            existing = genes_to_sort.get(gene)
            if not existing:
                genes_to_sort[gene] = sort_tuple
                continue

            if sort_tuple < existing:
                genes_to_sort[gene] = sort_tuple

    sorted_genes = sorted(genes_to_sort.keys(),
                          key=lambda gene: genes_to_sort[gene])
    if plot_top_n < 0:
        return sorted_genes

    return sorted_genes[:plot_top_n]


def parse_sorted_gene_line(string):
    gene, line = string.split('\t', maxsplit=1)
    if gene == 'None':
        return None, line

    without_quotes = gene.strip('"')
    return without_quotes, line


def sort_file_by_genes(genes, orig_path, sorted_path, get_gene_from_line):
    # header and comment lines will be put under None in lines_by_gene
    lines_by_gene = {None: list()}
    for gene in genes:
        lines_by_gene[gene] = list()

    with open(orig_path, 'rt') as in_handle:
        for line_i, line in enumerate(in_handle):
            gene_for_line = get_gene_from_line(line_i, line)
            if gene_for_line in lines_by_gene:
                lines_by_gene[gene_for_line].append(line)

    with open(sorted_path, 'wt') as out_handle:
        # write header and comment lines first
        for gene in [None] + genes:
            for line in lines_by_gene[gene]:
                # Write the gene in quotes or None before the original line.
                # This allows identifying the gene for each line later.
                if gene is None:
                    out_handle.write('None')
                else:
                    out_handle.write('"{}"'.format(gene))

                out_handle.write('\t')
                out_handle.write(line)


def get_gene_from_line_gtf(line_i, line):
    parsed = rmats_long_utils.parse_gtf_line(line)
    if parsed is None:
        return None

    gene_id = parsed['attributes'].get('gene_id')
    return gene_id


def get_gene_from_line_abundance(line_i, line):
    if line_i == 0:
        return None

    columns = line.rstrip('\n').split('\t')
    return columns[2]


def get_gene_from_line_diff_transcripts(line_i, line):
    if line_i == 0:
        return None

    columns = line.rstrip('\n').split('\t')
    return columns[0]


# The abundance file will be filtered down to each gene when running
# scripts for that gene. The visualization needs to know the total read count
# for each sample in order to calculate CPM. Add extra lines for each gene
# with the total read count for other genes for each sample. Then when the
# file is filtered to a single gene, the CPM can still be calculated.
def get_other_gene_abundance_totals_by_gene(genes, abundance):
    rows_by_gene = dict()
    parsed = rmats_long_utils.parse_abundance_file(abundance)
    total_by_sample = parsed['total_by_sample']
    total_by_gene_by_sample = parsed['total_by_gene_by_sample']
    sample_names = parsed['sample_names']
    # ['transcript_ID', 'transcript_name', 'gene_ID']
    shared_columns = ['NA', 'NA', 'NA']
    for gene in genes:
        sample_columns = list()
        for sample in sample_names:
            total = total_by_sample[sample]
            gene_total = total_by_gene_by_sample[gene][sample]
            non_gene_total = total - gene_total
            sample_columns.append(
                rmats_long_utils.format_float(non_gene_total))

        rows_by_gene[gene] = shared_columns + sample_columns

    return rows_by_gene


def sort_files_by_genes(genes, abundance, updated_gtf, gencode_gtf,
                        diff_transcripts, temp_dir):
    sorted_paths = dict()
    sorted_abundance = os.path.join(temp_dir, 'abundance.esp')
    sorted_paths['abundance'] = sorted_abundance
    sort_file_by_genes(genes, abundance, sorted_abundance,
                       get_gene_from_line_abundance)

    sorted_updated_gtf = os.path.join(temp_dir, 'updated.gtf')
    sorted_paths['updated_gtf'] = sorted_updated_gtf
    sort_file_by_genes(genes, updated_gtf, sorted_updated_gtf,
                       get_gene_from_line_gtf)

    sorted_diff_transcripts = os.path.join(temp_dir, 'diff_transcripts.tsv')
    sorted_paths['diff_transcripts'] = sorted_diff_transcripts
    sort_file_by_genes(genes, diff_transcripts, sorted_diff_transcripts,
                       get_gene_from_line_diff_transcripts)
    if gencode_gtf:
        sorted_gencode_gtf = os.path.join(temp_dir, 'gencode.gtf')
        sorted_paths['gencode_gtf'] = sorted_gencode_gtf
        sort_file_by_genes(genes, gencode_gtf, sorted_gencode_gtf,
                           get_gene_from_line_gtf)

    return sorted_paths


def rmats_long(args):
    rmats_long_utils.create_output_dir(args.out_dir)
    if args.diff_transcripts:
        count_significant_isoforms(args.out_dir, args.diff_transcripts,
                                   args.adj_pvalue, args.delta_proportion)
        diff_transcripts_path = args.diff_transcripts
    else:
        detect_differential_isoforms(args.abundance, args.out_dir,
                                     args.group_1, args.group_2,
                                     args.adj_pvalue, args.delta_proportion,
                                     args.num_threads)
        diff_transcripts_path = os.path.join(args.out_dir,
                                             'differential_transcripts.tsv')

    filtered_diff_transcripts_path = os.path.join(
        args.out_dir, 'differential_transcripts_filtered.tsv')

    if args.plot_top_n == 0:
        return

    genes_to_process = get_top_n_genes(filtered_diff_transcripts_path,
                                       args.plot_top_n)
    with tempfile.TemporaryDirectory(suffix='_tmp',
                                     prefix='rmats_long_',
                                     dir=args.out_dir) as temp_dir:
        other_gene_abundance_totals = get_other_gene_abundance_totals_by_gene(
            genes_to_process, args.abundance)
        sorted_paths = sort_files_by_genes(genes_to_process, args.abundance,
                                           args.updated_gtf, args.gencode_gtf,
                                           diff_transcripts_path, temp_dir)
        process_genes(genes_to_process, temp_dir, sorted_paths,
                      other_gene_abundance_totals, args.out_dir, args.group_1,
                      args.group_1_name, args.group_2, args.group_2_name,
                      args.plot_file_type)


def main():
    args = parse_args()
    rmats_long(args)
    print('finished')


if __name__ == '__main__':
    main()
