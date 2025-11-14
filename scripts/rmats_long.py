import argparse
import datetime
import math
import multiprocessing
import os
import os.path
import shutil
import tempfile

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(description=(
        'Identify significant splicing changes and produce plots'))
    gene_group_str = 'Gene isoforms'
    asm_group_str = 'ASM isoforms'
    use_either_str = 'Use either {} or {}'
    gene_group = parser.add_argument_group(
        gene_group_str, use_either_str.format(gene_group_str, asm_group_str))
    asm_group = parser.add_argument_group(
        asm_group_str, use_either_str.format(asm_group_str, gene_group_str))
    required_group = parser.add_argument_group('Required')
    optional_group = parser.add_argument_group('Optional')
    gene_group.add_argument(
        '--abundance', help='The path to the abundance.esp file from ESPRESSO')
    gene_group.add_argument(
        '--updated-gtf', help='The path to the updated.gtf file from ESPRESSO')
    gene_group.add_argument('--use-drimseq',
                            action='store_true',
                            help='Use DRIMSeq instead of run_stat_model.py')
    asm_group.add_argument(
        '--event-dir',
        help='The output directory from detect_splicing_events.py')
    asm_group.add_argument('--asm-counts',
                           help='The output file from count_reads_for_asms.py')
    asm_group.add_argument(
        '--align-dir',
        help=('The output directory from'
              ' organize_alignment_info_by_gene_and_chr.py'))
    asm_group.add_argument(
        '--gtf-dir',
        help='The output directory from organize_gene_info_by_chr.py')
    required_group.add_argument(
        '--group-1',
        required=True,
        help=('The path to a file listing the sample names for group 1. The'
              ' file should have a single line with the sample names as a'
              ' comma separated list. The sample names should match with the'
              ' ESPRESSO abundance column names or --asm-counts names.'))
    required_group.add_argument(
        '--group-2',
        required=True,
        help='The path to a file listing the sample names for group 2')
    required_group.add_argument('--out-dir',
                                required=True,
                                help='The path to use as the output directory')
    optional_group.add_argument(
        '--gencode-gtf',
        help=('The path to a gencode annotation.gtf file. Will be used to'
              ' identify the Ensembl canonical isoform and the gene name'))
    optional_group.add_argument(
        '--group-1-name',
        default='group 1',
        help='A name for group 1 (default %(default)s)')
    optional_group.add_argument(
        '--group-2-name',
        default='group 2',
        help='A name for group 2 (default %(default)s)')
    optional_group.add_argument(
        '--num-threads',
        type=int,
        default=1,
        help='The number of threads to use (default %(default)s)')
    optional_group.add_argument(
        '--process-top-n',
        type=int,
        help=('Generate plots and classify isoform differences for the'
              ' top "n" significant genes. By default all significant'
              ' genes are processed'))
    optional_group.add_argument(
        '--process-selected',
        help=('A comma separated list of gene IDs or ASM IDs to'
              ' generate plots and classify isoform differences for'))
    optional_group.add_argument(
        '--plot-file-type',
        choices=['.pdf', '.png', 'all'],
        default='.png',
        help='The file type for output plots (default %(default)s))')
    optional_group.add_argument(
        '--intron-scaling',
        type=int,
        default=1,
        help=('The factor to use to reduce intron length in the plot.'
              ' A value of 2 would reduce introns to 1/2 of the'
              ' original plot length (default %(default)s)'))
    optional_group.add_argument(
        '--max-transcripts',
        type=int,
        default=7,
        help=('How many transcripts to plot individually.'
              ' The remaining transcripts in the gene will be grouped together'
              ' (default %(default)s)'))
    optional_group.add_argument(
        '--diff-transcripts',
        help=('The path to the differential transcript results. If given then'
              ' skip the differential isoform calculation.'))
    optional_group.add_argument(
        '--adj-pvalue',
        type=float,
        default=0.05,
        help='The cutoff for adjusted p-value (default %(default)s)')
    optional_group.add_argument(
        '--use-unadjusted-pvalue',
        action='store_true',
        help='Use pvalue instead of adj_pvalue for the cutoff')
    optional_group.add_argument(
        '--delta-proportion',
        type=float,
        default=0.05,
        help='The cutoff for delta isoform proportion (default %(default)s)')
    optional_group.add_argument(
        '--compare-all-within-gene',
        action='store_true',
        help=('Compare the most significant isoform to all other isoforms in'
              ' the gene. By default, the most significant isoform is only'
              ' compared to the most significant isoform with a delta'
              ' proportion in the opposite direction.'))
    optional_group.add_argument(
        '--covar-tsv',
        help=('A .tsv with 1 line per sample. The first line has the column'
              ' names. The first column is sample_id. Each additional column'
              ' is a covariate.'))
    optional_group.add_argument(
        '--min-cpm-per-asm',
        type=float,
        default=0,
        help=('Only consider ASMs where at least 1 sample has at least'
              ' this CPM of reads assigned to the ASM.'
              ' (default %(default)s)'))
    optional_group.add_argument('--no-splice-graph-plot',
                                action='store_true',
                                help='Do not run plot_splice_graph.py')
    optional_group.add_argument(
        '--min-isoform-reads',
        type=int,
        default=1,
        help=('Only consider isoforms with at least this many reads'
              ' (default %(default)s)'))
    optional_group.add_argument(
        '--limit-asm-to-top-n-isoforms',
        type=int,
        default=50,
        help=('Only consider the top N isoforms with the highest'
              ' total proportion across samples for each ASM'
              ' (default %(default)s)'))
    optional_group.add_argument(
        '--average-reads-per-group',
        type=float,
        default=10,
        help=('For each sample group require the average read count to be at'
              ' least this value in order to be significant'
              ' (default %(default)s)'))
    optional_group.add_argument(
        '--average-cpm-per-group',
        type=float,
        help=('For each sample group require the average CPM to be at least'
              ' this value in order to be significant'))
    optional_group.add_argument(
        '--min-cpm-per-group',
        type=float,
        help=('For each sample group require the min CPM to be at least'
              ' this value in order to be significant'))
    optional_group.add_argument(
        '--asm-proportion-of-gene',
        type=float,
        default=0.05,
        help=('Require the ASM CPM to be at least this proportion of the'
              ' gene CPM in at least 1 sample in order to be significant'
              ' (default %(default)s)'))

    script_command_line = rmats_long_utils.get_script_command_line()
    args = parser.parse_args()
    gene_args = [args.abundance, args.updated_gtf]
    asm_args = [args.event_dir, args.asm_counts, args.align_dir, args.gtf_dir]
    any_gene_arg = any(x is not None for x in gene_args)
    any_asm_arg = any(x is not None for x in asm_args)
    if not (any_gene_arg or any_asm_arg):
        parser.error('Either arguments for {} or {} are required'.format(
            gene_group_str, asm_group_str))

    if any_gene_arg and any_asm_arg:
        print('Warning: either the arguments for {} or {}'
              ' should be used'.format(gene_group_str, asm_group_str))

    if any_gene_arg and (args.abundance is None):
        parser.error('--abundance is required when using Gene arguments')

    if any_asm_arg:
        if args.align_dir is None:
            parser.error('--align-dir is required when using ASM arguments')
        if args.event_dir is None:
            parser.error('--event-dir is required when using ASM arguments')

    if args.process_selected:
        args.process_selected = args.process_selected.split(',')

    return {'script_command_line': script_command_line, 'args': args}


def process_genes(genes, temp_dir, sorted_paths, other_gene_abundance_totals,
                  out_dir, group_1, group_1_name, group_2, group_2_name,
                  plot_file_type, intron_scaling, max_transcripts,
                  compare_all_within_gene, num_threads, summary):
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
                            group_2_name, plot_file_type, intron_scaling,
                            max_transcripts, compare_all_within_gene,
                            num_threads, summary)
                else:
                    process_genes_with_handles(
                        genes, temp_dir, abundance_handle, updated_gtf_handle,
                        diff_transcripts_handle, None,
                        other_gene_abundance_totals, out_dir, group_1,
                        group_1_name, group_2, group_2_name, plot_file_type,
                        intron_scaling, max_transcripts,
                        compare_all_within_gene, num_threads, summary)


def process_genes_with_handles(genes, temp_dir, abundance_handle,
                               updated_gtf_handle, diff_transcripts_handle,
                               gencode_gtf_handle, other_gene_abundance_totals,
                               out_dir, group_1, group_1_name, group_2,
                               group_2_name, plot_file_type, intron_scaling,
                               max_transcripts, compare_all_within_gene,
                               num_threads, summary):
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
    thread_details = create_gene_threads(num_threads)
    try:
        for gene in genes:
            temp_files_for_gene = write_temp_files_for_gene(
                gene, temp_dir, abundance_info, updated_gtf_info,
                diff_transcripts_info, gencode_gtf_info,
                other_gene_abundance_totals)
            process_gene(gene, temp_files_for_gene['abundance'],
                         temp_files_for_gene['updated_gtf'],
                         temp_files_for_gene['gencode_gtf'],
                         temp_files_for_gene['diff_transcripts'], out_dir,
                         group_1, group_1_name, group_2, group_2_name,
                         plot_file_type, intron_scaling, max_transcripts,
                         compare_all_within_gene, summary, thread_details)

        process_genes_on_main_thread(thread_details, summary)
    finally:
        cleanup_gene_threads(thread_details)

    combine_summary_values(thread_details, summary)


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
    with tempfile.NamedTemporaryFile('wt',
                                     prefix='gene_',
                                     suffix='_abundance.esp',
                                     dir=temp_dir,
                                     delete=False) as handle:
        temp_abundance_path = handle.name

    write_temp_file_for_gene(gene, temp_abundance_path, abundance_info)
    add_other_gene_totals_to_abundance(other_gene_abundance_totals[gene],
                                       temp_abundance_path)

    with tempfile.NamedTemporaryFile('wt',
                                     prefix='gene_',
                                     suffix='_updated.gtf',
                                     dir=temp_dir,
                                     delete=False) as handle:
        temp_updated_gtf_path = handle.name

    write_temp_file_for_gene(gene, temp_updated_gtf_path, updated_gtf_info)

    with tempfile.NamedTemporaryFile('wt',
                                     prefix='gene_',
                                     suffix='_gencode.gtf',
                                     dir=temp_dir,
                                     delete=False) as handle:
        temp_gencode_gtf_path = handle.name

    write_temp_file_for_gene(gene, temp_gencode_gtf_path, gencode_gtf_info)

    with tempfile.NamedTemporaryFile('wt',
                                     prefix='gene_',
                                     suffix='_diff_transcripts.tsv',
                                     dir=temp_dir,
                                     delete=False) as handle:
        temp_diff_transcripts_path = handle.name

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
                 group_2_name, plot_file_type, intron_scaling, max_transcripts,
                 compare_all_within_gene, summary, thread_details):
    gene_out_dir = os.path.join(out_dir, 'results_by_gene', gene_id)
    rmats_long_utils.create_output_dir(gene_out_dir, check_empty=False)
    gene_name = None
    if thread_details['threads']:
        thread_args = {
            'gene_name': gene_name,
            'gene_out_dir': gene_out_dir,
            'gene_id': gene_id,
            'diff_transcripts': diff_transcripts,
            'abundance': abundance,
            'updated_gtf': updated_gtf,
            'gencode_gtf': gencode_gtf,
            'group_1': group_1,
            'group_1_name': group_1_name,
            'group_2': group_2,
            'group_2_name': group_2_name,
            'plot_file_type': plot_file_type,
            'intron_scaling': intron_scaling,
            'max_transcripts': max_transcripts,
            'compare_all_within_gene': compare_all_within_gene
        }
        while True:
            was_put = rmats_long_utils.try_put_to_queue_with_short_wait(
                thread_details['in_queue'], thread_args)
            if was_put:
                break

            # The queue could be full due to thread errors
            rmats_long_utils.raise_exception_if_thread_exited_early(
                thread_details['threads'])
    else:
        process_gene_thread(gene_name, gene_out_dir, gene_id, diff_transcripts,
                            abundance, updated_gtf, gencode_gtf, group_1,
                            group_1_name, group_2, group_2_name,
                            plot_file_type, intron_scaling, max_transcripts,
                            compare_all_within_gene, summary)


def count_significant_isoforms(out_dir, diff_transcripts, adj_pvalue,
                               use_unadjusted_pvalue, delta_proportion,
                               average_reads_per_group, average_cpm_per_group,
                               min_cpm_per_group, asm_proportion_of_gene,
                               group_1_path, group_2_path, is_asm):
    python_executable = rmats_long_utils.get_python_executable()
    script_dir = rmats_long_utils.get_script_dir()
    count_script = os.path.join(script_dir, 'count_significant_isoforms.py')
    if is_asm:
        out_tsv = os.path.join(out_dir, 'differential_isoforms_filtered.tsv')
        diff_asms = os.path.join(out_dir, 'differential_asms.tsv')
    else:
        out_tsv = os.path.join(out_dir,
                               'differential_transcripts_filtered.tsv')
        diff_asms = os.path.join(out_dir, 'differential_genes.tsv')

    command = [
        python_executable, count_script, '--diff-transcripts',
        diff_transcripts, '--diff-asms', diff_asms, '--out-tsv', out_tsv,
        '--adj-pvalue',
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


def detect_differential_isoforms(
        abundance, asm_counts, out_dir, group_1, group_2, adj_pvalue,
        use_unadjusted_pvalue, delta_proportion, num_threads, covar_tsv,
        use_drimseq, min_isoform_reads, min_cpm_per_asm,
        limit_asm_to_top_n_isoforms, average_reads_per_group,
        average_cpm_per_group, min_cpm_per_group, asm_proportion_of_gene,
        gene_cpm_tsv, sample_read_total_path):
    python_executable = rmats_long_utils.get_python_executable()
    script_dir = rmats_long_utils.get_script_dir()
    detect_script = os.path.join(script_dir, 'detect_differential_isoforms.py')
    command = [
        python_executable, detect_script, '--out-dir', out_dir, '--group-1',
        group_1, '--group-2', group_2, '--adj-pvalue',
        str(adj_pvalue), '--delta-proportion',
        str(delta_proportion), '--num-threads',
        str(num_threads), '--min-isoform-reads',
        str(min_isoform_reads), '--min-cpm-per-asm',
        str(min_cpm_per_asm), '--sample-read-total-tsv',
        sample_read_total_path, '--limit-asm-to-top-n-isoforms',
        str(limit_asm_to_top_n_isoforms), '--average-reads-per-group',
        str(average_reads_per_group), '--gene-cpm-tsv', gene_cpm_tsv,
        '--asm-proportion-of-gene',
        str(asm_proportion_of_gene)
    ]
    if abundance:
        command.extend(['--abundance', abundance])
    else:
        command.extend(['--asm-counts', asm_counts])

    if use_unadjusted_pvalue:
        command.append('--use-unadjusted-pvalue')

    if covar_tsv:
        command.extend(['--covar-tsv', covar_tsv])

    if use_drimseq:
        command.append('--use-drimseq')

    if average_cpm_per_group is not None:
        command.extend(['--average-cpm-per-group', str(average_cpm_per_group)])

    if min_cpm_per_group is not None:
        command.extend(['--min-cpm-per-group', str(min_cpm_per_group)])

    rmats_long_utils.run_command(command)


def visualize_asm_in_gene(event_dir, chr_id, asm_id_str, gencode_gtf, asm_gtf,
                          diff_transcripts, gene_out_dir, plot_file_type,
                          intron_scaling, max_transcripts, gene_id, gene_name,
                          asm_strand, asm_start, asm_end):
    if gencode_gtf is None:
        return

    start_coord = None
    if asm_start != 'source':
        start_coord = asm_start

    end_coord = None
    if asm_end != 'sink':
        end_coord = asm_end

    if asm_strand == '-':
        start_coord, end_coord = end_coord, start_coord

    graph_path = rmats_long_utils.get_graph_file_path(event_dir, chr_id)
    with tempfile.TemporaryDirectory(suffix='_tmp',
                                     prefix='visualize_asm_in_gene_',
                                     dir=gene_out_dir) as temp_dir:
        python_executable = rmats_long_utils.get_python_executable()
        script_dir = rmats_long_utils.get_script_dir()
        visualize_script = os.path.join(script_dir, 'visualize_isoforms.py')
        command = [
            python_executable, visualize_script, '--gene-id', gene_id,
            '--gene-name', gene_name, '--asm-id', asm_id_str, '--asm-gtf',
            asm_gtf, '--gencode-gtf', gencode_gtf, '--updated-gtf',
            gencode_gtf, '--diff-transcripts', diff_transcripts, '--out-dir',
            temp_dir, '--plot-file-type', plot_file_type, '--intron-scaling',
            str(intron_scaling), '--max-transcripts',
            str(max_transcripts), '--graph-file', graph_path
        ]
        if start_coord is not None:
            command.extend(['--start-coord', str(start_coord)])

        if end_coord is not None:
            command.extend(['--end-coord', str(end_coord)])

        rmats_long_utils.run_command(command)

        files = os.listdir(temp_dir)
        rename_pairs = [('_structure.', '{}_in_gene{}'),
                        ('_structure_sj_counts.', '{}_in_gene_sj_counts{}')]
        for name in files:
            for rename_pair in rename_pairs:
                expected_partial_name = rename_pair[0]
                new_name_template = rename_pair[1]
                if expected_partial_name not in name:
                    continue

                tmp_path = os.path.join(temp_dir, name)
                _, suffix = os.path.splitext(name)
                new_name = new_name_template.format(asm_id_str, suffix)
                new_path = os.path.join(gene_out_dir, new_name)
                shutil.move(tmp_path, new_path)


def visualize_isoforms(event_dir, chr_id, gene_id, abundance, updated_gtf,
                       gencode_gtf, diff_transcripts, out_dir, group_1,
                       group_1_name, group_2, group_2_name, plot_file_type,
                       intron_scaling, max_transcripts, asm_id_str, gene_name,
                       is_asm):
    python_executable = rmats_long_utils.get_python_executable()
    script_dir = rmats_long_utils.get_script_dir()
    visualize_script = os.path.join(script_dir, 'visualize_isoforms.py')
    command = [
        python_executable, visualize_script, '--gene-id', gene_id,
        '--abundance', abundance, '--updated-gtf', updated_gtf,
        '--diff-transcripts', diff_transcripts, '--out-dir', out_dir,
        '--group-1', group_1, '--group-2', group_2, '--group-1-name',
        group_1_name, '--group-2-name', group_2_name, '--plot-file-type',
        plot_file_type, '--intron-scaling',
        str(intron_scaling), '--max-transcripts',
        str(max_transcripts)
    ]
    if gencode_gtf:
        command.extend(['--gencode-gtf', gencode_gtf])
    if gene_name:
        command.extend(['--gene-name', gene_name])
    if is_asm:
        command.append('--is-asm')
    if event_dir:
        graph_path = rmats_long_utils.get_graph_file_path(event_dir, chr_id)
        command.extend(['--graph-file', graph_path])
    if asm_id_str:
        command.extend(['--asm-id', asm_id_str])

    rmats_long_utils.run_command(command)


def classify_isoforms(gene, out_dir, updated_gtf, gencode_gtf,
                      main_transcript_id, second_transcript_id, summary):
    python_executable = rmats_long_utils.get_python_executable()
    script_dir = rmats_long_utils.get_script_dir()
    classify_script = os.path.join(script_dir,
                                   'classify_isoform_differences.py')
    if second_transcript_id is None:
        # comparing to all other isoforms in the gene
        out_tsv = os.path.join(
            out_dir, '{}_isoform_differences_from_{}.tsv'.format(
                gene, main_transcript_id))
        command = [
            python_executable, classify_script, '--updated-gtf', updated_gtf,
            '--out-tsv', out_tsv, '--main-transcript-id', main_transcript_id
        ]
    else:
        out_tsv = os.path.join(
            out_dir, '{}_isoform_differences_{}_to_{}.tsv'.format(
                gene, main_transcript_id, second_transcript_id))
        command = [
            python_executable, classify_script, '--updated-gtf', updated_gtf,
            '--out-tsv', out_tsv, '--main-transcript-id', main_transcript_id,
            '--second-transcript-id', second_transcript_id
        ]

    if gencode_gtf:
        command.extend(['--gencode-gtf', gencode_gtf])

    rmats_long_utils.run_command(command)
    summarize_classification(summary, out_tsv)


def get_gene_id_order_for_process_selected(ids_to_process,
                                           diff_transcripts_path):
    gene_id_order = list()
    gene_ids_seen = set()
    gene_by_id = dict()
    with open(diff_transcripts_path, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            gene_id = row['gene_id']
            if gene_id not in gene_ids_seen:
                gene_id_order.append(gene_id)
                gene_ids_seen.add(gene_id)

            if gene_id in ids_to_process:
                gene_by_id[gene_id] = gene_id
                continue

            asm_id = row.get('asm_id')
            if asm_id in ids_to_process:
                gene_by_id[asm_id] = gene_id

    for id_to_process in ids_to_process:
        gene_id = gene_by_id.get(id_to_process)
        if not gene_id:
            raise Exception('{} not found in {}'.format(
                id_to_process, diff_transcripts_path))

    selected_genes = set(gene_by_id.values())
    final_gene_id_order = list()
    for gene_id in gene_id_order:
        if gene_id in selected_genes:
            final_gene_id_order.append(gene_id)

    return final_gene_id_order


def get_top_n_ids_to_process(diff_transcripts, process_top_n,
                             use_unadjusted_pvalue, is_asm):
    gene_id_order = list()
    gene_ids_seen = set()
    ids_to_sort = dict()
    with open(diff_transcripts, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            if use_unadjusted_pvalue:
                pvalue = rmats_long_utils.parse_float(row['pvalue'])
            else:
                pvalue = rmats_long_utils.parse_float(row['adj_pvalue'])

            if math.isnan(pvalue):
                pvalue = 2  # to be last in sort

            abs_delta = abs(
                rmats_long_utils.parse_float(row['delta_isoform_proportion']))
            if math.isnan(abs_delta):
                abs_delta = -1  # to be last in sort

            gene_id = row['gene_id']
            if gene_id not in gene_ids_seen:
                gene_id_order.append(gene_id)
                gene_ids_seen.add(gene_id)

            if is_asm:
                row_id = row['asm_id']
            else:
                row_id = gene_id

            sort_tuple = (pvalue, -abs_delta, row_id, gene_id)
            existing = ids_to_sort.get(row_id)
            if (not existing) or (sort_tuple < existing):
                ids_to_sort[row_id] = sort_tuple

    sorted_ids = sorted(ids_to_sort.keys(),
                        key=lambda row_id: ids_to_sort[row_id])
    if (process_top_n is not None) and process_top_n >= 0:
        sorted_ids = sorted_ids[:process_top_n]

    selected_genes = set()
    for row_id in sorted_ids:
        sort_tuple = ids_to_sort[row_id]
        gene_id = sort_tuple[3]
        selected_genes.add(gene_id)

    final_gene_id_order = list()
    for gene_id in gene_id_order:
        if gene_id in selected_genes:
            final_gene_id_order.append(gene_id)

    return {'ids': sorted_ids, 'gene_order': final_gene_id_order}


def write_filtered_gtf_for_gene_from_handle(gene_id, gene_gtf_path, in_handle):
    # in_handle could have extra genes.
    # Skip other genes until gene_id is found.
    # After writing lines for gene_id and finding a new gene,
    # seek back to the start of the line
    found_gene = False
    with open(gene_gtf_path, 'wt') as out_handle:
        while True:
            offset = in_handle.tell()
            line = in_handle.readline()
            if not line:
                return

            gene, line = parse_sorted_gene_line(line)
            if gene is None:
                continue

            if gene != gene_id:
                if found_gene:
                    in_handle.seek(offset)
                    return

                continue

            found_gene = True
            out_handle.write(line)


def sort_gtf_by_gene_id(orig_gtf_path, gene_id_order, temp_dir):
    if not orig_gtf_path:
        return None

    sorted_path = os.path.join(temp_dir, 'sorted_by_gene_id.gtf')
    sort_file_by_genes(gene_id_order, orig_gtf_path, sorted_path,
                       get_genes_from_line_gtf)
    return sorted_path


def parse_sorted_gene_line(string):
    gene, line = string.split('\t', maxsplit=1)
    if gene == 'None':
        return None, line

    without_quotes = gene.strip('"')
    return without_quotes, line


def sort_file_by_genes(genes, orig_path, sorted_path, get_genes_from_line):
    # header and comment lines will be put under None in lines_by_gene
    lines_by_gene = {None: list()}
    for gene in genes:
        lines_by_gene[gene] = list()

    with open(orig_path, 'rt') as in_handle:
        for line_i, line in enumerate(in_handle):
            genes_for_line = get_genes_from_line(line_i, line)
            for gene_for_line in genes_for_line:
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


def get_genes_from_line_gtf(line_i, line):
    parsed = rmats_long_utils.parse_gtf_line(line)
    if parsed is None:
        return [None]

    gene_id_str = parsed['attributes'].get('gene_id')
    gene_ids = gene_id_str.split(',')
    return gene_ids


def get_genes_from_line_abundance(line_i, line):
    if line_i == 0:
        return [None]

    columns = rmats_long_utils.read_tsv_line(line)
    gene_id_str = columns[2]
    gene_ids = gene_id_str.split(',')
    return gene_ids


def get_genes_from_line_diff_transcripts(line_i, line):
    if line_i == 0:
        return [None]

    columns = rmats_long_utils.read_tsv_line(line)
    gene_id = columns[0]
    return [gene_id]


def write_sample_read_total_from_abundance(total_by_sample, read_total_path):
    with open(read_total_path, 'wt') as read_handle:
        headers = ['sample', 'total']
        rmats_long_utils.write_tsv_line(read_handle, headers)
        for sample, total in total_by_sample.items():
            # require the total read count to be an integer
            rounded_total = math.ceil(total)
            columns = [sample, str(rounded_total)]
            rmats_long_utils.write_tsv_line(read_handle, columns)


def write_gene_cpm_from_abundance(total_by_gene_by_sample, total_by_sample,
                                  gene_cpm_path):
    genes = sorted(total_by_gene_by_sample.keys())
    samples = sorted(total_by_sample.keys())
    headers = ['gene_id']
    headers.extend(samples)
    with open(gene_cpm_path, 'wt') as cpm_handle:
        rmats_long_utils.write_tsv_line(cpm_handle, headers)
        for gene in genes:
            columns = [gene]
            gene_total_by_sample = total_by_gene_by_sample[gene]
            for sample in samples:
                gene_total = gene_total_by_sample.get(sample, 0)
                sample_total = total_by_sample[sample]
                if sample_total == 0:
                    cpm = 0
                else:
                    cpm = gene_total * (1e6 / sample_total)

                formatted = rmats_long_utils.format_float(cpm)
                columns.append(formatted)

            rmats_long_utils.write_tsv_line(cpm_handle, columns)


def write_sample_read_total_and_gene_cpm_from_abundance(abundance, out_dir):
    read_total_path = os.path.join(out_dir, 'sample_read_totals.tsv')
    gene_cpm_path = os.path.join(out_dir, 'sample_gene_cpm.tsv')
    parsed = rmats_long_utils.parse_abundance_file(abundance)
    total_by_sample = parsed['total_by_sample']
    total_by_gene_by_sample = parsed['total_by_gene_by_sample']
    write_sample_read_total_from_abundance(total_by_sample, read_total_path)
    write_gene_cpm_from_abundance(total_by_gene_by_sample, total_by_sample,
                                  gene_cpm_path)
    return {'read_total': read_total_path, 'gene_cpm': gene_cpm_path}


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
                       get_genes_from_line_abundance)

    sorted_updated_gtf = os.path.join(temp_dir, 'updated.gtf')
    sorted_paths['updated_gtf'] = sorted_updated_gtf
    sort_file_by_genes(genes, updated_gtf, sorted_updated_gtf,
                       get_genes_from_line_gtf)

    sorted_diff_transcripts = os.path.join(temp_dir, 'diff_transcripts.tsv')
    sorted_paths['diff_transcripts'] = sorted_diff_transcripts
    sort_file_by_genes(genes, diff_transcripts, sorted_diff_transcripts,
                       get_genes_from_line_diff_transcripts)
    if gencode_gtf:
        sorted_gencode_gtf = os.path.join(temp_dir, 'gencode.gtf')
        sorted_paths['gencode_gtf'] = sorted_gencode_gtf
        sort_file_by_genes(genes, gencode_gtf, sorted_gencode_gtf,
                           get_genes_from_line_gtf)

    return sorted_paths


# Can be run with duplicate_asms to overwrite previous summary values
def summarize_differential_transcripts(summary,
                                       filtered_diff_transcripts_path,
                                       adj_pvalue,
                                       delta_proportion,
                                       is_asm,
                                       duplicate_asms=None):
    if duplicate_asms is None:
        values_for_dups = dict()
    else:
        values_for_dups = dict()
        # All asm_ids in the same duplicate set will have the same list.
        # The list will keep the priority order of IDs in the duplicate set.
        # For IDs that are part of a duplicate set, update the values in the
        # list when the row is read. After all rows are read, the values
        # for the highest priority ID in each set will be used.
        for dup_ids in duplicate_asms.values():
            values = list()
            for asm_id in dup_ids:
                values_for_dups[asm_id] = values
                values.append({'asm': asm_id, 'isoforms': list()})

    summary['pvalue_threshold'] = adj_pvalue
    summary['delta_proportion_threshold'] = delta_proportion
    sig_asms = set()
    sig_genes = set()
    sig_isoforms = set()
    with open(filtered_diff_transcripts_path, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            if is_asm:
                asm_id = row['asm_id']
                dup_values = values_for_dups.get(asm_id)
                if dup_values:
                    for asm_dict in dup_values:
                        if asm_dict['asm'] == asm_id:
                            asm_dict['gene'] = row['gene_id']
                            asm_dict['isoforms'].append(row['isoform_id'])
                            break

                    continue

                sig_asms.add(asm_id)
                sig_genes.add(row['gene_id'])
                sig_isoforms.add(row['isoform_id'])
            else:
                sig_genes.add(row['gene_id'])
                sig_isoforms.add(row['feature_id'])

    for dup_values in values_for_dups.values():
        for asm_dict in dup_values:
            if 'gene' in asm_dict:
                sig_asms.add(asm_dict['asm'])
                sig_genes.add(asm_dict['gene'])
                for isoform in asm_dict['isoforms']:
                    sig_isoforms.add(isoform)

                break

        # Clear the list so that it is only processed once
        dup_values.clear()

    if is_asm:
        summary['significant_asms'] = len(sig_asms)

    summary['significant_genes'] = len(sig_genes)
    summary['significant_isoforms'] = len(sig_isoforms)


def initialize_summary_classifications(summary):
    classifications = summary.get('classifications')
    if not classifications:
        classifications = {
            'total': 0,
            'SE': 0,
            'A5SS': 0,
            'A3SS': 0,
            'MXE': 0,
            'RI': 0,
            'AFE': 0,
            'ALE': 0,
            'COMPLEX': 0,
            'combo': 0,
            'ENDS': 0
        }
        summary['classifications'] = classifications

    return classifications


def summarize_classification(summary, out_tsv):
    classifications = initialize_summary_classifications(summary)
    by_pair = dict()
    with open(out_tsv, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            t_1 = row['transcript1']
            t_2 = row['transcript2']
            pair = (t_1, t_2)
            event = row['event']
            if pair in by_pair:
                by_pair[pair] = 'combo'
            else:
                by_pair[pair] = event

    if not by_pair:
        return

    total = len(by_pair)
    classifications['total'] += total
    for classification in by_pair.values():
        classifications[classification] += 1


def write_summary(summary, summary_path, use_unadjusted_pvalue, command_line,
                  commit_id):
    is_asm = summary.get('significant_asms') is not None
    with open(summary_path, 'wt') as handle:
        handle.write('## {}\n'.format(command_line))
        handle.write('## source code commit: {}\n'.format(commit_id))
        handle.write('## significant differential isoform usage\n')
        handle.write('total significant isoforms: {}\n'.format(
            summary['significant_isoforms']))
        handle.write('total genes with significant isoforms: {}\n'.format(
            summary['significant_genes']))
        if is_asm:
            handle.write('total ASMs with significant isoforms: {}\n'.format(
                summary['significant_asms']))

        pvalue_line = 'pvalue threshold: {}\n'.format(
            summary['pvalue_threshold'])
        if not use_unadjusted_pvalue:
            pvalue_line = 'adjusted {}'.format(pvalue_line)

        handle.write(pvalue_line)

        handle.write('delta isoform proportion threshold: {}\n'.format(
            summary['delta_proportion_threshold']))
        classifications = summary.get('classifications')
        if classifications:
            handle.write(
                '## alternative splicing classifications between isoform pairs\n'
            )
            handle.write('total classified isoform pairs: {}\n'.format(
                classifications['total']))
            handle.write('exon skipping: {}\n'.format(classifications['SE']))
            handle.write("alternative 5'-splice site: {}\n".format(
                classifications['A5SS']))
            handle.write("alternative 3'-splice site: {}\n".format(
                classifications['A3SS']))
            handle.write('mutually exclusive exons: {}\n'.format(
                classifications['MXE']))
            handle.write('intron retention: {}\n'.format(
                classifications['RI']))
            handle.write('alternative first exon: {}\n'.format(
                classifications['AFE']))
            handle.write('alternative last exon: {}\n'.format(
                classifications['ALE']))
            handle.write('complex: {}\n'.format(classifications['COMPLEX']))
            handle.write('combinatorial: {}\n'.format(
                classifications['combo']))
            handle.write('alternative endpoints: {}\n'.format(
                classifications['ENDS']))

        if is_asm:
            handle.write('## Number of isoforms per ASM\n')
            count_by_num_isoforms = summary.get('asms_by_num_isoforms', dict())
            min_num = 2
            max_num = 2
            for num_isoforms in count_by_num_isoforms.keys():
                min_num = min(min_num, num_isoforms)
                max_num = max(max_num, num_isoforms)

            for num_isoforms in range(2, max_num + 1):
                count = count_by_num_isoforms.get(num_isoforms, 0)
                handle.write('total ASMs with {} isoforms: {}\n'.format(
                    num_isoforms, count))

    return summary_path


def visualize_summary(summary_path, summary_dir):
    script_dir = rmats_long_utils.get_script_dir()
    visualize_script = os.path.join(script_dir, 'visualize_summary.R')
    command = ['Rscript', visualize_script, summary_path, summary_dir]
    rmats_long_utils.run_command(command)


def parse_asm_id(asm_id):
    parts = asm_id.split('_')
    if len(parts) != 2:
        raise Exception('could not parse asm_id: {}'.format(asm_id))

    chr_id = parts[0]
    event_i = int(parts[1])
    return {'chr': chr_id, 'event_i': event_i}


def sort_asm_ids_by_chr(asm_ids):
    by_chr = dict()
    for asm in asm_ids:
        parsed = parse_asm_id(asm)
        chr_id = parsed['chr']
        for_chr = rmats_long_utils.try_get_or_set_default(
            by_chr, chr_id, list())
        for_chr.append(parsed)

    for events in by_chr.values():
        events.sort(key=lambda x: x['event_i'])

    return by_chr


def load_diff_asm_isoforms_by_chr(path):
    by_chr = dict()
    with open(path, 'rt') as handle:
        rows = rmats_long_utils.row_iterator_for_tsv_with_header(handle)
        for row in rows:
            parsed_event = parse_asm_id(row['asm_id'])
            chr_id = parsed_event['chr']
            for_chr = rmats_long_utils.try_get_or_set_default(
                by_chr, chr_id, list())
            for_chr.append(row)

    return by_chr


def load_gene_id_to_name(chr_id, gtf_dir):
    gene_id_to_name = dict()
    if gtf_dir is None:
        return gene_id_to_name

    path = rmats_long_utils.get_chr_file_path(gtf_dir, chr_id)
    with open(path, 'rt') as handle:
        for line in handle:
            columns = rmats_long_utils.read_tsv_line(line)
            # Gene lines have 4 columns
            if len(columns) != 4:
                continue

            gene_id = columns[2]
            gene_name = columns[3]
            if not gene_name:
                continue

            gene_id_to_name[gene_id] = gene_name

    return gene_id_to_name


def write_abundance_for_asm_id(asm_id, asm_id_str, temp_dir, samples,
                               read_count_rows, isoform_ids, total_by_sample):
    headers = ['transcript_ID', 'transcript_name', 'gene_ID']
    headers.extend(samples)
    asm_chr = asm_id['chr']
    asm_i = asm_id['event_i']
    out_columns = [None, 'NA', asm_id_str]
    remaining_for_samples = list()
    for sample in samples:
        out_columns.append('0')
        sample_total = total_by_sample[sample]
        remaining_for_samples.append(sample_total)

    counts_by_isoform_by_sample = dict()
    # Initialize all isoform_ids for this event.
    # That way any isoforms not found in read_count_rows will still
    # have a line in the abundance file.
    for isoform_id in isoform_ids:
        counts_by_isoform_by_sample[isoform_id] = dict()

    for row in read_count_rows:
        isoform = row['isoform_id']
        sample = row['sample_id']
        count = float(row['count'])
        counts_by_sample = rmats_long_utils.try_get_or_set_default(
            counts_by_isoform_by_sample, isoform, dict())
        counts_by_sample[sample] = count

    prefix = '{}_{}'.format(asm_chr, asm_i)
    with tempfile.NamedTemporaryFile('wt',
                                     prefix=prefix,
                                     suffix='.tsv',
                                     dir=temp_dir,
                                     delete=False) as handle:
        abundance = handle.name
        rmats_long_utils.write_tsv_line(handle, headers)
        for isoform, counts_by_sample in counts_by_isoform_by_sample.items():
            out_columns[0] = isoform
            for sample_i, sample in enumerate(samples):
                value = counts_by_sample.get(sample, 0)
                formatted = rmats_long_utils.format_float(value, digits=2)
                out_columns[sample_i + 3] = formatted
                remaining_for_samples[sample_i] -= value

            rmats_long_utils.write_tsv_line(handle, out_columns)

        other_abundance_columns = ['NA', 'NA', 'NA']
        for remaining in remaining_for_samples:
            other_abundance_columns.append(str(remaining))

        rmats_long_utils.write_tsv_line(handle, other_abundance_columns)

    return abundance


def write_gtf_transcript_or_exon_line(feature_name, gene_id, transcript_id,
                                      chr_name, start, end, strand, handle):
    columns = list()
    columns.append(chr_name)
    columns.append('.')  # source
    columns.append(feature_name)
    columns.append(str(start))
    columns.append(str(end))
    columns.append('.')  # score
    columns.append(strand)
    columns.append('.')  # frame
    attributes_str = 'gene_id "{}";transcript_id "{}";'.format(
        gene_id, transcript_id)
    columns.append(attributes_str)
    rmats_long_utils.write_tsv_line(handle, columns)


def write_gtf_exon_line(gene_id, transcript_id, chr_name, start, end, strand,
                        handle):
    write_gtf_transcript_or_exon_line('exon', gene_id, transcript_id, chr_name,
                                      start, end, strand, handle)


def write_gtf_transcript_line(gene_id, transcript_id, chr_name, start, end,
                              strand, handle):
    write_gtf_transcript_or_exon_line('transcript', gene_id, transcript_id,
                                      chr_name, start, end, strand, handle)


def write_gtf_for_asm_id(asm_id, chr_name, temp_dir, event_reader,
                         num_isoforms_by_asm, summary):
    asm_chr = asm_id['chr']
    asm_i = asm_id['event_i']
    gene_id = None
    asm_strand = None
    asm_start = None
    asm_end = None
    asm_id_str = '{}_{}'.format(asm_chr, asm_i)
    isoform_ids = list()
    with tempfile.NamedTemporaryFile('wt',
                                     prefix=asm_id_str,
                                     suffix='.gtf',
                                     dir=temp_dir,
                                     delete=False) as handle:
        gtf = handle.name
        while True:
            row = event_reader.get_next_row()
            if row is None:
                break

            row_id = row['asm_id']
            parsed_event = parse_asm_id(row_id)
            if parsed_event['chr'] != asm_chr:
                raise Exception(
                    'Unexpected chr found in event_reader {} {}'.format(
                        asm_id, row))

            if parsed_event['event_i'] > asm_i:
                break

            event_reader.mark_row_as_used()
            if parsed_event['event_i'] < asm_i:
                continue

            gene_id = event_reader.gene_id
            asm_strand = row['strand']
            asm_start = row['start']
            asm_end = row['end']
            isoforms_str = row['isoforms']
            isoforms = rmats_long_utils.parse_isoforms_str(isoforms_str)
            isoform_ids_str = row['isoform_ids']
            isoform_ids = rmats_long_utils.parse_isoform_ids_str(
                isoform_ids_str)
            by_num_isoforms = rmats_long_utils.try_get_or_set_default(
                summary, 'asms_by_num_isoforms', dict())
            num_isoforms = len(isoforms)
            num_isoforms_by_asm[asm_id_str] = num_isoforms
            old_count = by_num_isoforms.get(num_isoforms, 0)
            by_num_isoforms[num_isoforms] = old_count + 1
            for isoform_i, isoform in enumerate(isoforms):
                isoform_id = isoform_ids[isoform_i]
                coord_i = 0
                last_start_i = len(isoform) - 2
                transcript_start = isoform[coord_i]
                transcript_end = isoform[last_start_i + 1]
                write_gtf_transcript_line(asm_id_str, isoform_id, chr_name,
                                          transcript_start, transcript_end,
                                          asm_strand, handle)
                while coord_i <= last_start_i:
                    exon_start = isoform[coord_i]
                    exon_end = isoform[coord_i + 1]
                    coord_i += 2
                    write_gtf_exon_line(asm_id_str, isoform_id, chr_name,
                                        exon_start, exon_end, asm_strand,
                                        handle)

    return {
        'gtf': gtf,
        'gene_id': gene_id,
        'isoform_ids': isoform_ids,
        'asm_strand': asm_strand,
        'asm_start': asm_start,
        'asm_end': asm_end
    }


def write_diff_transcripts_for_asm_id(asm_id, temp_dir, diff_info):
    headers = None
    with tempfile.NamedTemporaryFile('wt',
                                     prefix='diff',
                                     suffix='.tsv',
                                     dir=temp_dir,
                                     delete=False) as handle:
        tsv_path = handle.name
        for info in diff_info:
            parsed_asm_id = parse_asm_id(info['asm_id'])
            if parsed_asm_id != asm_id:
                continue

            if headers is None:
                headers = list()
                renamed_headers = list()
                for key in info.keys():
                    headers.append(key)
                    if key == 'asm_id':
                        renamed_headers.append('gene_id')
                    elif key == 'gene_id':
                        renamed_headers.append('orig_gene_id')
                    elif key == 'isoform_id':
                        renamed_headers.append('feature_id')
                    else:
                        renamed_headers.append(key)

                rmats_long_utils.write_tsv_line(handle, renamed_headers)

            columns = list()
            for header in headers:
                columns.append(info[header])

            rmats_long_utils.write_tsv_line(handle, columns)

    return tsv_path


def plot_splice_graph(event_dir, chr_name, gene_id, asm_id, out_path):
    python_executable = rmats_long_utils.get_python_executable()
    script_dir = rmats_long_utils.get_script_dir()
    plot_script = os.path.join(script_dir, 'plot_splice_graph.py')
    command = [
        python_executable, plot_script, '--event-dir', event_dir, '--chr',
        chr_name, '--gene-id', gene_id, '--asm-id', asm_id, '--out-file',
        out_path
    ]
    rmats_long_utils.run_command(command)


def process_asm_event_with_count_rows(
        asm_id, temp_dir, id_to_chr, read_count_rows, samples, diff_info,
        event_dir, event_reader, gene_id_to_name, seen_gene_ids,
        total_by_sample, gencode_gtf_handle, out_dir, group_1, group_1_name,
        group_2, group_2_name, plot_file_type, intron_scaling, max_transcripts,
        no_splice_graph_plot, num_isoforms_by_asm, summary, thread_details):
    asm_chr = asm_id['chr']
    asm_i = asm_id['event_i']
    asm_id_str = '{}_{}'.format(asm_chr, asm_i)
    chr_id = int(asm_chr)
    chr_name = id_to_chr[chr_id]
    gtf_result = write_gtf_for_asm_id(asm_id, chr_name, temp_dir, event_reader,
                                      num_isoforms_by_asm, summary)
    updated_gtf = gtf_result['gtf']
    gene_id = gtf_result['gene_id']
    isoform_ids = gtf_result['isoform_ids']
    asm_strand = gtf_result['asm_strand']
    asm_start = gtf_result['asm_start']
    asm_end = gtf_result['asm_end']
    gene_name = gene_id_to_name.get(gene_id, gene_id)
    gene_out_dir = os.path.join(out_dir, 'results_by_gene', gene_id)
    rmats_long_utils.create_output_dir(gene_out_dir, check_empty=False)

    gene_gtf_path = None
    if gencode_gtf_handle is not None:
        gene_gtf_path = os.path.join(gene_out_dir, 'filtered.gtf')
        if gene_id not in seen_gene_ids:
            write_filtered_gtf_for_gene_from_handle(gene_id, gene_gtf_path,
                                                    gencode_gtf_handle)

    seen_gene_ids.add(gene_id)

    if thread_details['threads']:
        thread_args = {
            'asm_id': asm_id,
            'asm_id_str': asm_id_str,
            'asm_strand': asm_strand,
            'asm_start': asm_start,
            'asm_end': asm_end,
            'temp_dir': temp_dir,
            'read_count_rows': read_count_rows,
            'samples': samples,
            'diff_info': diff_info,
            'event_dir': event_dir,
            'total_by_sample': total_by_sample,
            'gene_out_dir': gene_out_dir,
            'group_1': group_1,
            'group_1_name': group_1_name,
            'group_2': group_2,
            'group_2_name': group_2_name,
            'plot_file_type': plot_file_type,
            'intron_scaling': intron_scaling,
            'max_transcripts': max_transcripts,
            'no_splice_graph_plot': no_splice_graph_plot,
            'gene_id': gene_id,
            'isoform_ids': isoform_ids,
            'gene_name': gene_name,
            'chr_name': chr_name,
            'chr_id': chr_id,
            'updated_gtf': updated_gtf,
            'gencode_gtf': gene_gtf_path,
        }
        while True:
            was_put = rmats_long_utils.try_put_to_queue_with_short_wait(
                thread_details['in_queue'], thread_args)
            if was_put:
                break

            # The queue could be full due to thread errors
            rmats_long_utils.raise_exception_if_thread_exited_early(
                thread_details['threads'])
    else:
        process_asm_event_thread(asm_id, asm_id_str, asm_strand, asm_start,
                                 asm_end, temp_dir, read_count_rows, samples,
                                 diff_info, event_dir, total_by_sample,
                                 gene_out_dir, group_1, group_1_name, group_2,
                                 group_2_name, plot_file_type, intron_scaling,
                                 max_transcripts, no_splice_graph_plot,
                                 gene_id, isoform_ids, gene_name, chr_name,
                                 chr_id, updated_gtf, gene_gtf_path, summary)


class EventReader:
    def __init__(self, event_dir):
        self.event_dir = event_dir
        self.next_line = None
        self.handle = None
        self.headers = list()
        self.gene_i = None
        self.gene_id = None

    def close(self):
        if self.handle is not None:
            self.handle.close()
            self.handle = None
            self.next_line = None
            self.headers = list()
            self.gene_i = None
            self.gene_id = None

    def initialize_for_chr(self, chr_id):
        self.close()
        path = rmats_long_utils.get_chr_file_path(self.event_dir, chr_id)
        self.handle = open(path, 'rt')
        first_line = self.handle.readline()
        if not first_line:
            return

        columns = rmats_long_utils.read_tsv_line(first_line)
        self.headers = columns

    def get_next_row(self):
        if self.next_line is None:
            self.next_line = self.handle.readline()

        if not self.next_line:
            return None

        columns = rmats_long_utils.read_tsv_line(self.next_line)
        row = dict(zip(self.headers, columns))
        gene_i = row['gene_i']
        gene_id = row['gene_id']
        if (self.gene_i is None) or (self.gene_i != gene_i):
            self.gene_i = gene_i
            self.gene_id = gene_id

        return row

    def mark_row_as_used(self):
        self.next_line = None


def combine_summary_values(thread_details, summary):
    threads = thread_details['threads']
    if not threads:
        return

    remaining_summaries = len(threads)
    classifications = initialize_summary_classifications(summary)
    out_queue = thread_details['out_queue']
    while remaining_summaries != 0:
        other = rmats_long_utils.try_get_from_queue_with_short_wait(out_queue)
        if other is None:
            # The worker threads have already been joined.
            # The summaries should already be on the out_queue.
            raise Exception(
                'Missing {} summary file(s)'.format(remaining_summaries))

        remaining_summaries -= 1
        other_dict = other.get('classifications')
        if not other_dict:
            continue

        for k, v in other_dict.items():
            old_val = classifications.get(k)
            classifications[k] = old_val + v


def rm_gene_temp_files(abundance, updated_gtf, gencode_gtf, diff_transcripts):
    paths = [abundance, updated_gtf, gencode_gtf, diff_transcripts]
    for path in paths:
        if path:
            os.remove(path)


def process_gene_thread(gene_name, gene_out_dir, gene_id, diff_transcripts,
                        abundance, updated_gtf, gencode_gtf, group_1,
                        group_1_name, group_2, group_2_name, plot_file_type,
                        intron_scaling, max_transcripts,
                        compare_all_within_gene, summary):
    is_asm = False
    no_event_dir = None
    no_chr_id = None
    no_asm_id = None
    visualize_isoforms(no_event_dir, no_chr_id, gene_id, abundance,
                       updated_gtf, gencode_gtf, diff_transcripts,
                       gene_out_dir, group_1, group_1_name, group_2,
                       group_2_name, plot_file_type, intron_scaling,
                       max_transcripts, no_asm_id, gene_name, is_asm)

    significant_transcripts = rmats_long_utils.select_significant_transcripts(
        gene_id, diff_transcripts)
    most_significant = significant_transcripts[0]
    if compare_all_within_gene:
        classify_isoforms(gene_id, gene_out_dir, updated_gtf, gencode_gtf,
                          most_significant, None, summary)
        rm_gene_temp_files(abundance, updated_gtf, gencode_gtf,
                           diff_transcripts)
        return

    if len(significant_transcripts) != 2:
        rm_gene_temp_files(abundance, updated_gtf, gencode_gtf,
                           diff_transcripts)
        return

    second_transcript = significant_transcripts[1]
    classify_isoforms(gene_id, gene_out_dir, updated_gtf, gencode_gtf,
                      most_significant, second_transcript, summary)
    rm_gene_temp_files(abundance, updated_gtf, gencode_gtf, diff_transcripts)


def run_gene_thread(in_queue, out_queue, signal_queue):
    summary = dict()
    while True:
        signal = rmats_long_utils.try_get_from_queue_without_wait(signal_queue)
        if signal is not None:
            out_queue.put(summary)
            return

        args = rmats_long_utils.try_get_from_queue_with_short_wait(in_queue)
        if args is None:
            continue

        # After the main thread has put all the actual gene args on the
        # in_queue, it will put an int for each thread (including the main thread).
        # Those values allow checking if the in_queue has been fully processed.
        if isinstance(args, int):
            break

        process_gene_thread(args['gene_name'], args['gene_out_dir'],
                            args['gene_id'], args['diff_transcripts'],
                            args['abundance'], args['updated_gtf'],
                            args['gencode_gtf'], args['group_1'],
                            args['group_1_name'], args['group_2'],
                            args['group_2_name'], args['plot_file_type'],
                            args['intron_scaling'], args['max_transcripts'],
                            args['compare_all_within_gene'], summary)

    signal_queue.get()
    out_queue.put(summary)


def try_enqueue_thread_ints(remaining_thread_i, in_queue):
    while remaining_thread_i >= 0:
        was_put = rmats_long_utils.try_put_to_queue_with_short_wait(
            in_queue, remaining_thread_i)
        if not was_put:
            break

        remaining_thread_i -= 1

    return remaining_thread_i


def process_genes_on_main_thread(thread_details, summary):
    threads = thread_details['threads']
    if not threads:
        return

    in_queue = thread_details['in_queue']
    remaining_thread_i = len(threads)
    # Loop until one of the ints is read back from in_queue
    while True:
        remaining_thread_i = try_enqueue_thread_ints(remaining_thread_i,
                                                     in_queue)
        args = rmats_long_utils.try_get_from_queue_with_short_wait(in_queue)
        if args is None:
            continue

        if isinstance(args, int):
            rmats_long_utils.drain_queue(in_queue)
            return

        process_gene_thread(args['gene_name'], args['gene_out_dir'],
                            args['gene_id'], args['diff_transcripts'],
                            args['abundance'], args['updated_gtf'],
                            args['gencode_gtf'], args['group_1'],
                            args['group_1_name'], args['group_2'],
                            args['group_2_name'], args['plot_file_type'],
                            args['intron_scaling'], args['max_transcripts'],
                            args['compare_all_within_gene'], summary)


def create_gene_threads(num_threads):
    threads = list()
    details = {
        'threads': threads,
        'in_queue': None,
        'out_queue': None,
        'signal_queue': None
    }
    if num_threads == 1:
        return details

    # Every thread should be able to have something on the queue.
    # queue_size_mult allows some extra room
    queue_size_mult = 10
    num_workers = num_threads - 1
    signal_queue = multiprocessing.Queue(num_workers)
    queue_size = num_workers * queue_size_mult
    in_queue = multiprocessing.Queue(queue_size)
    out_queue = multiprocessing.Queue(queue_size)
    details['in_queue'] = in_queue
    details['out_queue'] = out_queue
    details['signal_queue'] = signal_queue
    for _ in range(num_workers):
        thread = multiprocessing.Process(target=run_gene_thread,
                                         args=(in_queue, out_queue,
                                               signal_queue))
        threads.append(thread)
        thread.start()

    return details


def cleanup_gene_threads(details):
    threads = details['threads']
    if not threads:
        return

    in_queue = details['in_queue']
    signal_queue = details['signal_queue']
    for _ in threads:
        signal_queue.put(True)

    rmats_long_utils.drain_queue(in_queue)

    for thread in threads:
        thread.join()
        if thread.exitcode != 0:
            raise Exception('thread exited with value: {}'.format(
                thread.exitcode))


def process_asm_event_thread(asm_id, asm_id_str, asm_strand, asm_start,
                             asm_end, temp_dir, read_count_rows, samples,
                             diff_info, event_dir, total_by_sample,
                             gene_out_dir, group_1, group_1_name, group_2,
                             group_2_name, plot_file_type, intron_scaling,
                             max_transcripts, no_splice_graph_plot, gene_id,
                             isoform_ids, gene_name, chr_name, chr_id,
                             updated_gtf, gencode_gtf, summary):
    abundance = write_abundance_for_asm_id(asm_id, asm_id_str, temp_dir,
                                           samples, read_count_rows,
                                           isoform_ids, total_by_sample)
    diff_transcripts = write_diff_transcripts_for_asm_id(
        asm_id, temp_dir, diff_info)
    no_gencode_gtf = None
    is_asm = True
    visualize_isoforms(event_dir, chr_id, gene_id, abundance, updated_gtf,
                       no_gencode_gtf, diff_transcripts, gene_out_dir, group_1,
                       group_1_name, group_2, group_2_name, plot_file_type,
                       intron_scaling, max_transcripts, asm_id_str, gene_name,
                       is_asm)

    visualize_asm_in_gene(event_dir, chr_id, asm_id_str, gencode_gtf,
                          updated_gtf, diff_transcripts, gene_out_dir,
                          plot_file_type, intron_scaling, max_transcripts,
                          gene_id, gene_name, asm_strand, asm_start, asm_end)

    significant_transcripts = rmats_long_utils.select_significant_transcripts(
        asm_id_str, diff_transcripts)
    if len(significant_transcripts) == 2:
        most_significant = significant_transcripts[0]
        second_transcript = significant_transcripts[1]
        classify_isoforms(asm_id_str, gene_out_dir, updated_gtf,
                          no_gencode_gtf, most_significant, second_transcript,
                          summary)

    plot_file_types = [plot_file_type]
    if plot_file_type == 'all':
        plot_file_types = ['.pdf', '.png']

    for file_type in plot_file_types:
        graph_out = os.path.join(gene_out_dir,
                                 '{}_graph{}'.format(asm_id_str, file_type))
        if not no_splice_graph_plot:
            plot_splice_graph(event_dir, chr_name, gene_id, asm_id_str,
                              graph_out)

    os.remove(abundance)
    os.remove(updated_gtf)


def run_asm_thread(in_queue, out_queue, signal_queue):
    summary = dict()
    while True:
        signal = rmats_long_utils.try_get_from_queue_without_wait(signal_queue)
        if signal is not None:
            out_queue.put(summary)
            return

        args = rmats_long_utils.try_get_from_queue_with_short_wait(in_queue)
        if args is None:
            continue

        # Similar to run_gene_thread, expect an int to indicate that
        # the queue has been fully processed.
        if isinstance(args, int):
            break

        process_asm_event_thread(
            args['asm_id'], args['asm_id_str'], args['asm_strand'],
            args['asm_start'], args['asm_end'], args['temp_dir'],
            args['read_count_rows'], args['samples'], args['diff_info'],
            args['event_dir'], args['total_by_sample'], args['gene_out_dir'],
            args['group_1'], args['group_1_name'], args['group_2'],
            args['group_2_name'], args['plot_file_type'],
            args['intron_scaling'], args['max_transcripts'],
            args['no_splice_graph_plot'], args['gene_id'], args['isoform_ids'],
            args['gene_name'], args['chr_name'], args['chr_id'],
            args['updated_gtf'], args['gencode_gtf'], summary)

    signal_queue.get()
    out_queue.put(summary)


def process_asms_on_main_thread(thread_details, summary):
    threads = thread_details['threads']
    if not threads:
        return

    in_queue = thread_details['in_queue']
    remaining_thread_i = len(threads)
    # Loop until one of the ints is read back from in_queue
    while True:
        remaining_thread_i = try_enqueue_thread_ints(remaining_thread_i,
                                                     in_queue)
        args = rmats_long_utils.try_get_from_queue_with_short_wait(in_queue)
        if args is None:
            continue

        if isinstance(args, int):
            rmats_long_utils.drain_queue(in_queue)
            return

        process_asm_event_thread(
            args['asm_id'], args['asm_id_str'], args['asm_strand'],
            args['asm_start'], args['asm_end'], args['temp_dir'],
            args['read_count_rows'], args['samples'], args['diff_info'],
            args['event_dir'], args['total_by_sample'], args['gene_out_dir'],
            args['group_1'], args['group_1_name'], args['group_2'],
            args['group_2_name'], args['plot_file_type'],
            args['intron_scaling'], args['max_transcripts'],
            args['no_splice_graph_plot'], args['gene_id'], args['isoform_ids'],
            args['gene_name'], args['chr_name'], args['chr_id'],
            args['updated_gtf'], args['gencode_gtf'], summary)


def create_asm_threads(num_threads):
    threads = list()
    details = {
        'threads': threads,
        'in_queue': None,
        'out_queue': None,
        'signal_queue': None
    }
    if num_threads == 1:
        return details

    # Every thread should be able to have something on the queue.
    # queue_size_mult allows some extra room
    queue_size_mult = 10
    num_workers = num_threads - 1
    signal_queue = multiprocessing.Queue(num_workers)
    queue_size = num_workers * queue_size_mult
    in_queue = multiprocessing.Queue(queue_size)
    out_queue = multiprocessing.Queue(queue_size)
    details['in_queue'] = in_queue
    details['out_queue'] = out_queue
    details['signal_queue'] = signal_queue
    for _ in range(num_workers):
        thread = multiprocessing.Process(target=run_asm_thread,
                                         args=(in_queue, out_queue,
                                               signal_queue))
        threads.append(thread)
        thread.start()

    return details


def cleanup_asm_threads(details):
    threads = details['threads']
    if not threads:
        return

    in_queue = details['in_queue']
    signal_queue = details['signal_queue']
    for _ in threads:
        signal_queue.put(True)

    rmats_long_utils.drain_queue(in_queue)

    for thread in threads:
        thread.join()
        if thread.exitcode != 0:
            raise Exception('thread exited with value: {}'.format(
                thread.exitcode))


def process_top_asm_events(asm_ids, temp_dir, id_to_chr, out_dir, group_1,
                           group_1_name, group_2, group_2_name, plot_file_type,
                           intron_scaling, max_transcripts,
                           no_splice_graph_plot, diff_transcripts_path,
                           event_dir, gtf_dir, gencode_gtf, num_threads,
                           sample_read_total_path, num_isoforms_by_asm,
                           summary):
    total_by_sample = rmats_long_utils.parse_sample_totals_file(
        sample_read_total_path)
    event_reader = EventReader(event_dir)
    thread_details = create_asm_threads(num_threads)
    try:
        if gencode_gtf:
            with open(gencode_gtf, 'rt') as gencode_gtf_handle:
                process_top_asm_events_with_event_reader(
                    asm_ids, temp_dir, id_to_chr, out_dir, group_1,
                    group_1_name, group_2, group_2_name, plot_file_type,
                    intron_scaling, max_transcripts, no_splice_graph_plot,
                    diff_transcripts_path, event_dir, event_reader, gtf_dir,
                    gencode_gtf_handle, total_by_sample, num_isoforms_by_asm,
                    summary, thread_details)
        else:
            process_top_asm_events_with_event_reader(
                asm_ids, temp_dir, id_to_chr, out_dir, group_1, group_1_name,
                group_2, group_2_name, plot_file_type, intron_scaling,
                max_transcripts, no_splice_graph_plot, diff_transcripts_path,
                event_dir, event_reader, gtf_dir, None, total_by_sample,
                num_isoforms_by_asm, summary, thread_details)
    finally:
        event_reader.close()
        cleanup_asm_threads(thread_details)

    combine_summary_values(thread_details, summary)


def process_top_asm_events_with_event_reader(
        asm_ids, temp_dir, id_to_chr, out_dir, group_1, group_1_name, group_2,
        group_2_name, plot_file_type, intron_scaling, max_transcripts,
        no_splice_graph_plot, diff_transcripts_path, event_dir, event_reader,
        gtf_dir, gencode_gtf_handle, total_by_sample, num_isoforms_by_asm,
        summary, thread_details):
    asms_by_chr = sort_asm_ids_by_chr(asm_ids)
    diff_info_by_chr = load_diff_asm_isoforms_by_chr(diff_transcripts_path)
    current_chr = None
    chr_asms = None
    chr_diff_info = None
    gene_id_to_name = None
    seen_gene_ids = set()
    current_asm_id = None
    read_count_rows = list()
    samples = sorted(total_by_sample.keys())
    count_path = os.path.join(out_dir, 'count.tsv')
    with open(count_path, 'rt') as counts_handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(
                counts_handle):
            parsed_asm_id = parse_asm_id(row['asm_id'])
            row_chr = parsed_asm_id['chr']
            if (current_chr is None) or (row_chr != current_chr):
                if read_count_rows:
                    process_asm_event_with_count_rows(
                        current_asm_id, temp_dir, id_to_chr, read_count_rows,
                        samples, chr_diff_info, event_dir, event_reader,
                        gene_id_to_name, seen_gene_ids, total_by_sample,
                        gencode_gtf_handle, out_dir, group_1, group_1_name,
                        group_2, group_2_name, plot_file_type, intron_scaling,
                        max_transcripts, no_splice_graph_plot,
                        num_isoforms_by_asm, summary, thread_details)
                    read_count_rows = list()
                    current_asm_id = None

                current_chr = row_chr
                chr_asms = asms_by_chr.get(row_chr)
                chr_diff_info = diff_info_by_chr.get(row_chr)
                event_reader.initialize_for_chr(current_chr)
                gene_id_to_name = load_gene_id_to_name(current_chr, gtf_dir)
                seen_gene_ids = set()

            if not chr_asms:
                continue

            if (((current_asm_id is not None)
                 and (parsed_asm_id != current_asm_id))):
                if read_count_rows:
                    process_asm_event_with_count_rows(
                        current_asm_id, temp_dir, id_to_chr, read_count_rows,
                        samples, chr_diff_info, event_dir, event_reader,
                        gene_id_to_name, seen_gene_ids, total_by_sample,
                        gencode_gtf_handle, out_dir, group_1, group_1_name,
                        group_2, group_2_name, plot_file_type, intron_scaling,
                        max_transcripts, no_splice_graph_plot,
                        num_isoforms_by_asm, summary, thread_details)
                    read_count_rows = list()
                    current_asm_id = None

            if parsed_asm_id not in chr_asms:
                continue

            current_asm_id = parsed_asm_id
            read_count_rows.append(row)

    if read_count_rows:
        process_asm_event_with_count_rows(
            current_asm_id, temp_dir, id_to_chr, read_count_rows, samples,
            chr_diff_info, event_dir, event_reader, gene_id_to_name,
            seen_gene_ids, total_by_sample, gencode_gtf_handle, out_dir,
            group_1, group_1_name, group_2, group_2_name, plot_file_type,
            intron_scaling, max_transcripts, no_splice_graph_plot,
            num_isoforms_by_asm, summary, thread_details)

    process_asms_on_main_thread(thread_details, summary)


def get_gene_dirs_from_out_dir(out_dir):
    gene_dirs = list()
    by_gene_dir = os.path.join(out_dir, 'results_by_gene')
    if not os.path.exists(by_gene_dir):
        return gene_dirs

    names = os.listdir(by_gene_dir)
    for name in names:
        path = os.path.join(by_gene_dir, name)
        if os.path.isdir(path):
            gene_dirs.append(path)

    return gene_dirs


def parse_isoform_diff_file(name, path):
    expected_str = '_isoform_differences_'
    parts = name.split(expected_str)
    if len(parts) == 1:
        return None

    asm_id = parts[0]
    events = dict()
    with open(path, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            event_type = row['event']
            coords = row['coordinates']
            coords_for_type = events.get(event_type)
            if coords_for_type is None:
                coords_for_type = list()
                events[event_type] = coords_for_type

            coords_for_type.append(coords)

    event_types = list(events.keys())
    for event_type in event_types:
        list_value = events[event_type]
        events[event_type] = tuple(sorted(list_value))

    return {'asm_id': asm_id, 'events': events}


def check_for_duplicate_asms_in_isoform_diff_files(gene_dir,
                                                   ids_to_process_set,
                                                   summary):
    classifications = initialize_summary_classifications(summary)
    any_duplicate = False
    seen_events = dict()
    names = os.listdir(gene_dir)
    names.sort()
    for name in names:
        path = os.path.join(gene_dir, name)
        parsed_diff = parse_isoform_diff_file(name, path)
        if parsed_diff is None:
            continue

        asm_id = parsed_diff['asm_id']
        events = parsed_diff['events']
        event_key = tuple(sorted(events.items()))
        details = seen_events.get(event_key)
        if details:
            any_duplicate = True
            details['asms'].append(asm_id)
        else:
            seen_events[event_key] = {'asms': [asm_id], 'events': events}

    for event_key, details in seen_events.items():
        processed_asms = [
            x for x in details['asms'] if x in ids_to_process_set
        ]
        num_processed_dups = len(processed_asms) - 1
        if num_processed_dups < 1:
            continue

        events = details['events']
        if len(events) == 1:
            only_event_type = event_key[0][0]
            events_for_type = event_key[0][1]
            if len(events_for_type) == 1:
                event_type = only_event_type
            else:
                event_type = 'combo'
        else:
            event_type = 'combo'

        old_count = classifications.get(event_type)
        if old_count is not None:
            classifications[event_type] = old_count - num_processed_dups
            classifications['total'] -= num_processed_dups

    return {'any_dup': any_duplicate, 'seen_events': seen_events}


def create_or_update_num_isoforms_file(gene_dir, seen_events,
                                       num_isoforms_by_asm):
    num_isoforms_file = os.path.join(gene_dir, 'isoform_count_by_asm.tsv')
    gene_asm_ids = list()
    for details in seen_events.values():
        asm_ids = details['asms']
        gene_asm_ids.extend(asm_ids)

    gene_asm_id_set = set(gene_asm_ids)
    if os.path.exists(num_isoforms_file):
        with open(num_isoforms_file, 'rt') as in_handle:
            for row in rmats_long_utils.row_iterator_for_tsv_with_header(
                    in_handle):
                asm_id = row['asm_id']
                if asm_id not in gene_asm_id_set:
                    continue

                if asm_id in num_isoforms_by_asm:
                    continue

                num_isoforms_by_asm[asm_id] = int(row['isoform_count'])

    gene_asm_ids.sort()
    with open(num_isoforms_file, 'wt') as out_handle:
        headers = ['asm_id', 'isoform_count']
        rmats_long_utils.write_tsv_line(out_handle, headers)
        for asm_id in gene_asm_ids:
            isoform_count = num_isoforms_by_asm.get(asm_id)
            if isoform_count is None:
                continue

            rmats_long_utils.write_tsv_line(
                out_handle, [asm_id, str(isoform_count)])


def load_duplicate_asms_from_tsv(path, duplicate_asms):
    if not os.path.exists(path):
        return

    with open(path, 'rt') as handle:
        for line in handle:
            asm_ids = rmats_long_utils.read_tsv_line(line)
            for asm_id in asm_ids:
                duplicate_asms[asm_id] = asm_ids


def check_for_duplicate_asms_in_gene(gene_dir, gene_was_processed,
                                     num_isoforms_by_asm, asm_pvalues,
                                     duplicate_asms, ids_to_process_set,
                                     summary):
    dup_file = os.path.join(gene_dir, 'duplicate_asms.tsv')
    if not gene_was_processed:
        load_duplicate_asms_from_tsv(dup_file, duplicate_asms)
        return

    isoform_file_results = check_for_duplicate_asms_in_isoform_diff_files(
        gene_dir, ids_to_process_set, summary)
    seen_events = isoform_file_results['seen_events']
    create_or_update_num_isoforms_file(gene_dir, seen_events,
                                       num_isoforms_by_asm)
    if not isoform_file_results['any_dup']:
        return

    with open(dup_file, 'wt') as handle:
        sorted_event_keys = sorted(seen_events.keys())
        for event_key in sorted_event_keys:
            seen_event_details = seen_events[event_key]
            asms_for_event = seen_event_details['asms']
            asm_sort_keys = list()
            for asm_id in asms_for_event:
                num_isoforms = num_isoforms_by_asm[asm_id]
                pvalue = asm_pvalues[asm_id]
                if math.isnan(pvalue):
                    pvalue = 2  # sort NaN last

                asm_sort_keys.append((num_isoforms, pvalue, asm_id))

            if len(asm_sort_keys) <= 1:
                continue

            asm_sort_keys.sort()
            asm_ids = list()
            for sort_key in asm_sort_keys:
                asm_id = sort_key[2]
                asm_ids.append(asm_id)
                duplicate_asms[asm_id] = asm_ids

            rmats_long_utils.write_tsv_line(handle, asm_ids)


def check_for_duplicate_asms(out_dir, num_isoforms_by_asm, asm_pvalues,
                             gene_id_order, ids_to_process_set, summary):
    genes_to_process_set = set(gene_id_order)
    duplicate_asms = dict()
    gene_dirs = get_gene_dirs_from_out_dir(out_dir)
    for gene_dir in gene_dirs:
        gene_id = os.path.basename(gene_dir)
        gene_was_processed = gene_id in genes_to_process_set
        check_for_duplicate_asms_in_gene(gene_dir, gene_was_processed,
                                         num_isoforms_by_asm, asm_pvalues,
                                         duplicate_asms, ids_to_process_set,
                                         summary)

    return duplicate_asms


def update_summary_asm_num_isoforms_counts(num_isoforms_by_asm, duplicate_asms,
                                           ids_to_process_set, summary):
    by_num_isoforms = summary.get('asms_by_num_isoforms')
    if not by_num_isoforms:
        return

    for asm_id in ids_to_process_set:
        duplicates = duplicate_asms.get(asm_id)
        if not duplicates:
            continue

        processed_dups = [x for x in duplicates if x in ids_to_process_set]
        if asm_id == processed_dups[0]:
            continue

        num_isoforms = num_isoforms_by_asm[asm_id]
        old_count = by_num_isoforms.get(num_isoforms)
        if old_count is None:
            continue

        if old_count <= 1:
            del by_num_isoforms[num_isoforms]
            continue

        by_num_isoforms[num_isoforms] = old_count - 1


def read_asm_pvalues(out_dir):
    asm_pvalues = dict()
    diff_asms_path = os.path.join(out_dir, 'differential_asms.tsv')
    with open(diff_asms_path, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            asm_id = row['asm_id']
            pvalue_string = row['pvalue']
            pvalue = rmats_long_utils.parse_float(pvalue_string)
            asm_pvalues[asm_id] = pvalue

    return asm_pvalues


def append_is_strict_column_with_handles(in_handle, out_handle, event_reader):
    strict_header = 'is_strict'
    current_chr = None
    for line_i, line in enumerate(in_handle):
        columns = rmats_long_utils.read_tsv_line(line)
        if line_i == 0:
            orig_headers = columns
            if strict_header in orig_headers:
                return False

            new_headers = orig_headers + [strict_header]
            rmats_long_utils.write_tsv_line(out_handle, new_headers)
            continue

        diff_row = dict(zip(orig_headers, columns))
        asm_id = diff_row['asm_id']
        asm_id_details = parse_asm_id(asm_id)
        chr_id = asm_id_details['chr']
        if chr_id != current_chr:
            current_chr = chr_id
            event_reader.initialize_for_chr(current_chr)

        event_row = event_reader.get_next_row()
        is_strict = None
        while event_row:
            event_reader.mark_row_as_used()
            event_asm_id = event_row['asm_id']
            if event_asm_id == asm_id:
                is_strict = event_row[strict_header]
                break

            event_row = event_reader.get_next_row()

        out_columns = list()
        for header in orig_headers:
            out_columns.append(diff_row[header])

        out_columns.append(is_strict)
        rmats_long_utils.write_tsv_line(out_handle, out_columns)

    return True


def append_is_strict_column(out_dir, event_dir):
    diff_asms_path = os.path.join(out_dir, 'differential_asms.tsv')
    with tempfile.TemporaryDirectory(suffix='_tmp',
                                     prefix='append_is_strict_column_',
                                     dir=out_dir) as temp_dir:
        tmp_path = os.path.join(temp_dir, 'differential_asms.tsv')
        with open(diff_asms_path, 'rt') as in_handle:
            with open(tmp_path, 'wt') as out_handle:
                event_reader = EventReader(event_dir)
                try:
                    # May have already been added in a previous run
                    was_added = append_is_strict_column_with_handles(
                        in_handle, out_handle, event_reader)
                finally:
                    event_reader.close()

        if was_added:
            shutil.move(tmp_path, diff_asms_path)


def rmats_long(args):
    summary = dict()
    rmats_long_utils.create_output_dir(args.out_dir)
    is_asm = (args.abundance is None) and (args.align_dir is not None)
    if is_asm:
        sample_read_total_path = os.path.join(args.align_dir,
                                              'sample_read_totals.tsv')
        gene_cpm_tsv_path = os.path.join(args.align_dir, 'sample_gene_cpm.tsv')
    else:
        files_from_abun = write_sample_read_total_and_gene_cpm_from_abundance(
            args.abundance, args.out_dir)
        sample_read_total_path = files_from_abun['read_total']
        gene_cpm_tsv_path = files_from_abun['gene_cpm']

    if args.diff_transcripts:
        count_significant_isoforms(
            args.out_dir, args.diff_transcripts, args.adj_pvalue,
            args.use_unadjusted_pvalue, args.delta_proportion,
            args.average_reads_per_group, args.average_cpm_per_group,
            args.min_cpm_per_group, args.asm_proportion_of_gene, args.group_1,
            args.group_2, is_asm)
        diff_transcripts_path = args.diff_transcripts
    else:
        detect_differential_isoforms(
            args.abundance, args.asm_counts, args.out_dir, args.group_1,
            args.group_2, args.adj_pvalue, args.use_unadjusted_pvalue,
            args.delta_proportion, args.num_threads, args.covar_tsv,
            args.use_drimseq, args.min_isoform_reads, args.min_cpm_per_asm,
            args.limit_asm_to_top_n_isoforms, args.average_reads_per_group,
            args.average_cpm_per_group, args.min_cpm_per_group,
            args.asm_proportion_of_gene, gene_cpm_tsv_path,
            sample_read_total_path)
        if is_asm:
            diff_transcripts_path = os.path.join(args.out_dir,
                                                 'differential_isoforms.tsv')
        else:
            diff_transcripts_path = os.path.join(
                args.out_dir, 'differential_transcripts.tsv')

    if is_asm:
        filtered_diff_transcripts_path = os.path.join(
            args.out_dir, 'differential_isoforms_filtered.tsv')
        append_is_strict_column(args.out_dir, args.event_dir)
    else:
        filtered_diff_transcripts_path = os.path.join(
            args.out_dir, 'differential_transcripts_filtered.tsv')

    summarize_differential_transcripts(summary, filtered_diff_transcripts_path,
                                       args.adj_pvalue, args.delta_proportion,
                                       is_asm)

    if args.process_selected:
        ids_to_process = args.process_selected
        gene_id_order = get_gene_id_order_for_process_selected(
            ids_to_process, diff_transcripts_path)
    elif args.process_top_n == 0:
        return summary
    else:
        top_n_result = get_top_n_ids_to_process(filtered_diff_transcripts_path,
                                                args.process_top_n,
                                                args.use_unadjusted_pvalue,
                                                is_asm)
        ids_to_process = top_n_result['ids']
        gene_id_order = top_n_result['gene_order']

    with tempfile.TemporaryDirectory(suffix='_tmp',
                                     prefix='rmats_long_',
                                     dir=args.out_dir) as temp_dir:
        if is_asm:
            chr_id_mapping = rmats_long_utils.parse_chr_mapping(args.event_dir)
            id_to_chr = chr_id_mapping['id_to_chr']
            gencode_gtf = sort_gtf_by_gene_id(args.gencode_gtf, gene_id_order,
                                              temp_dir)
            num_isoforms_by_asm = dict()
            process_top_asm_events(
                ids_to_process, temp_dir, id_to_chr, args.out_dir,
                args.group_1, args.group_1_name, args.group_2,
                args.group_2_name, args.plot_file_type, args.intron_scaling,
                args.max_transcripts, args.no_splice_graph_plot,
                diff_transcripts_path, args.event_dir, args.gtf_dir,
                gencode_gtf, args.num_threads, sample_read_total_path,
                num_isoforms_by_asm, summary)
            asm_pvalues = read_asm_pvalues(args.out_dir)
            ids_to_process_set = set(ids_to_process)
            duplicate_asms = check_for_duplicate_asms(
                args.out_dir, num_isoforms_by_asm, asm_pvalues, gene_id_order,
                ids_to_process_set, summary)
            update_summary_asm_num_isoforms_counts(num_isoforms_by_asm,
                                                   duplicate_asms,
                                                   ids_to_process_set, summary)
            summarize_differential_transcripts(summary,
                                               filtered_diff_transcripts_path,
                                               args.adj_pvalue,
                                               args.delta_proportion,
                                               is_asm,
                                               duplicate_asms=duplicate_asms)
        else:
            other_gene_abundance_totals = get_other_gene_abundance_totals_by_gene(
                ids_to_process, args.abundance)
            sorted_paths = sort_files_by_genes(ids_to_process, args.abundance,
                                               args.updated_gtf,
                                               args.gencode_gtf,
                                               diff_transcripts_path, temp_dir)
            process_genes(ids_to_process, temp_dir, sorted_paths,
                          other_gene_abundance_totals, args.out_dir,
                          args.group_1, args.group_1_name, args.group_2,
                          args.group_2_name, args.plot_file_type,
                          args.intron_scaling, args.max_transcripts,
                          args.compare_all_within_gene, args.num_threads,
                          summary)

    return summary


def get_summary_path_and_dir(out_dir):
    summary_dir = out_dir
    summary_path = os.path.join(summary_dir, 'summary.txt')
    if os.path.exists(summary_path):
        datetime_format = '%Y_%m_%d_%H_%M_%S'
        now = datetime.datetime.now()
        datetime_string = now.strftime(datetime_format)
        summary_dir_name = 'summary_{}'.format(datetime_string)
        summary_dir = os.path.join(out_dir, summary_dir_name)
        rmats_long_utils.create_output_dir(summary_dir, check_empty=False)
        summary_path = os.path.join(summary_dir, 'summary.txt')

    return summary_path, summary_dir


def main():
    args_result = parse_args()
    args = args_result['args']
    command_line = args_result['script_command_line']
    commit_id = rmats_long_utils.get_source_code_commit()
    summary = rmats_long(args)
    summary_path, summary_dir = get_summary_path_and_dir(args.out_dir)
    write_summary(summary, summary_path, args.use_unadjusted_pvalue,
                  command_line, commit_id)
    visualize_summary(summary_path, summary_dir)
    print('\nsummary written to {}'.format(summary_path))


if __name__ == '__main__':
    main()
