import snakemake.utils

snakemake.utils.min_version('8.0.0')

configfile: 'snakemake_config.yaml'

onsuccess:
    print('workflow success')

onerror:
    print('workflow error')

DEFAULT_MEM_MB=4 * 1024  # 4 GB
DEFAULT_TIME_HOURS=12


def result_dir():
    return os.path.join('results', config['run_name'])


def all_input(wildcards):
    inputs = dict()
    quant_full = config.get('quantify_full_length_transcripts')
    quant_asms = config.get('quantify_asms')
    if quant_full:
        inputs['rmats_long_full_len'] = os.path.join(
            result_dir(), 'rmats_long_full_length', 'results_by_gene')

    if quant_asms:
        inputs['rmats_long_asm'] = os.path.join(
            result_dir(), 'rmats_long_asm', 'results_by_gene')
        inputs['asm_gtf'] = os.path.join(result_dir(), 'asm.gtf')

    if quant_full and quant_asms:
        inputs['venn'] = os.path.join(result_dir(), 'significant_genes_venn.png')

    return inputs


localrules: all
rule all:
    input:
        unpack(all_input),


def reference_file_wildcard_constraints():
    reference_files = config.get('reference_files')
    if reference_files:
        file_names = '|'.join([re.escape(file_name)
                               for file_name in reference_files])
        without_gz = '|'.join([re.escape(file_name[:-3])
                               for file_name in reference_files
                               if file_name.endswith('.gz')])
    else:
        no_match = '^$'  # only matches empty string
        file_names = no_match
        without_gz = no_match

    return {'file_names': file_names, 'without_gz': without_gz}


def get_url_for_download_reference_file(wildcards):
    file_name = wildcards.file_name
    return config['reference_files'][file_name]['url']


rule download_reference_file:
    output:
        ref_file=os.path.join('references', '{file_name}'),
    log:
        out=os.path.join('references',
                         'download_reference_file_{file_name}_log.out'),
        err=os.path.join('references',
                         'download_reference_file_{file_name}_log.err'),
    wildcard_constraints:
        file_name=reference_file_wildcard_constraints()['file_names']
    params:
        url=get_url_for_download_reference_file,
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        "curl -L '{params.url}'"
        ' -o {output.ref_file}'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule unzip_reference_file:
    input:
        gz=os.path.join('references', '{file_name}.gz'),
    output:
        un_gz=os.path.join('references', '{file_name}'),
    log:
        out=os.path.join('references',
                         'unzip_reference_file_{file_name}_log.out'),
        err=os.path.join('references',
                         'unzip_reference_file_{file_name}_log.err'),
    wildcard_constraints:
        file_name=reference_file_wildcard_constraints()['without_gz']
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time_hours=DEFAULT_TIME_HOURS,
    shell:
        ' gunzip -c {input.gz}'
        ' 1> {output.un_gz}'
        ' 2> {log.err}'

rule add_novel_transcripts_to_gtf:
    input:
        cleaned=os.path.join(result_dir(), 'espresso_cleaned.gtf'),
        ref=os.path.join('references', config['gtf_name']),
    output:
        combined=os.path.join(result_dir(), 'combined.gtf'),
    log:
        out=os.path.join(result_dir(), 'add_novel_transcripts_to_gtf_log.out'),
        err=os.path.join(result_dir(), 'add_novel_transcripts_to_gtf_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts', 'add_novel_transcripts_to_gtf.py'),
    resources:
        mem_mb=config['add_novel_transcripts_to_gtf_mem_gb'] * 1024,
        time_hours=config['add_novel_transcripts_to_gtf_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' --espresso-gtf {input.cleaned}'
        ' --other-gtf {input.ref}'
        ' --out-gtf {output.combined}'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule clean_espresso_gtf:
    input:
        gtf=config.get('espresso_gtf'),
    output:
        cleaned=os.path.join(result_dir(), 'espresso_cleaned.gtf'),
    log:
        out=os.path.join(result_dir(), 'clean_espresso_gtf_log.out'),
        err=os.path.join(result_dir(), 'clean_espresso_gtf_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts', 'clean_espresso_gtf.py'),
    resources:
        mem_mb=config['clean_espresso_gtf_mem_gb'] * 1024,
        time_hours=config['clean_espresso_gtf_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' --in-gtf {input.gtf}'
        ' --out-gtf {output.cleaned}'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule count_reads_for_asms:
    input:
        event_dir=os.path.join(result_dir(), '{mode}_events'),
        gtf_dir=os.path.join(result_dir(), 'annotation'),
        align_dir=os.path.join(result_dir(), 'alignments'),
    output:
        tsv=os.path.join(result_dir(), '{mode}_counts.tsv'),
    log:
        out=os.path.join(result_dir(), 'count_reads_for_asms_{mode}_log.out'),
        err=os.path.join(result_dir(), 'count_reads_for_asms_{mode}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts', 'count_reads_for_asms.py'),
    threads: config['count_reads_for_asms_threads']
    resources:
        mem_mb=config['count_reads_for_asms_mem_gb'] * 1024,
        time_hours=config['count_reads_for_asms_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' --event-dir {input.event_dir}'
        ' --gtf-dir {input.gtf_dir}'
        ' --align-dir {input.align_dir}'
        ' --out-tsv {output.tsv}'
        ' --num-threads {threads}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def detect_splicing_events_input(wildcards):
    inputs = dict()
    result_path = result_dir()
    inputs['gtf_dir'] = os.path.join(result_path, 'annotation')
    inputs['align_dir'] = os.path.join(result_path, 'alignments')

    return inputs


def min_reads_per_edge_param():
    min_reads = config.get('min_reads_per_edge')
    if min_reads in [None, '']:
        return ''

    return '--min-reads-per-edge {}'.format(min_reads)


def max_nodes_in_event_param():
    max_nodes = config.get('max_nodes_in_event')
    if max_nodes in [None, '']:
        return ''

    return '--max-nodes-in-event {}'.format(max_nodes)


def max_paths_in_event_param():
    max_paths = config.get('max_paths_in_event')
    if max_paths in [None, '']:
        return ''

    return '--max-paths-in-event {}'.format(max_paths)


def full_gene_asm_param(wildcards):
    if wildcards.mode == 'asm':
        return ''

    return '--output-full-gene-asm'


def simplify_gene_isoform_endpoints_param(wildcards):
    if wildcards.mode == 'asm':
        return ''

    simplify = config.get('simplify_gene_isoform_endpoints')
    if not simplify:
        return ''

    return '--simplify-gene-isoform-endpoints'


def filter_gene_isoforms_by_edge_param(wildcards):
    if wildcards.mode == 'asm':
        return ''

    filter_by_edge = config.get('filter_gene_isoforms_by_edge')
    if not filter_by_edge:
        return ''

    return '--filter-gene-isoforms-by-edge'


rule detect_splicing_events:
    input:
        unpack(detect_splicing_events_input),
    output:
        event_dir=directory(os.path.join(result_dir(), '{mode}_events')),
    log:
        out=os.path.join(result_dir(), 'detect_splicing_events_{mode}_log.out'),
        err=os.path.join(result_dir(), 'detect_splicing_events_{mode}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts', 'detect_splicing_events.py'),
        min_reads=min_reads_per_edge_param(),
        max_nodes=max_nodes_in_event_param(),
        max_paths=max_paths_in_event_param(),
        full_gene=full_gene_asm_param,
        simplify_gene_isoform_endpoints=simplify_gene_isoform_endpoints_param,
        filter_gene_isoforms_by_edge=filter_gene_isoforms_by_edge_param,
    threads: config['detect_splicing_events_threads']
    resources:
        mem_mb=config['detect_splicing_events_mem_gb'] * 1024,
        time_hours=config['detect_splicing_events_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' --gtf-dir {input.gtf_dir}'
        ' --align-dir {input.align_dir}'
        ' --out-dir {output.event_dir}'
        ' --num-threads {threads}'
        ' {params.min_reads}'
        ' {params.max_nodes}'
        ' {params.max_paths}'
        ' {params.full_gene}'
        ' {params.simplify_gene_isoform_endpoints}'
        ' {params.filter_gene_isoforms_by_edge}'
        ' 1> {log.out}'
        ' 2> {log.err}'

rule create_gtf_from_asm_definitions:
    input:
        event_dir=os.path.join(result_dir(), 'asm_events'),
    output:
        asm_gtf=os.path.join(result_dir(), 'asm.gtf'),
    log:
        out=os.path.join(result_dir(), 'create_gtf_from_asm_definitions_log.out'),
        err=os.path.join(result_dir(), 'create_gtf_from_asm_definitions_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts', 'create_gtf_from_asm_definitions.py'),
    resources:
        mem_mb=config['create_gtf_from_asm_definitions_mem_gb'] * 1024,
        time_hours=config['create_gtf_from_asm_definitions_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' --event-dir {input.event_dir}'
        ' --out-gtf {output.asm_gtf}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def num_inputs_by_sample():
    by_sample = dict()
    config_keys = ['group_1_samples', 'group_2_samples']
    for config_key in config_keys:
        by_name = config.get(config_key)
        for sample, details in by_name.items():
            num_inputs = len(details)
            by_sample[sample] = num_inputs

    return by_sample


localrules: write_samples_tsv
rule write_samples_tsv:
    output:
        samples_tsv=os.path.join(result_dir(), 'samples.tsv'),
    params:
        num_inputs_by_sample=num_inputs_by_sample(),
        result_dir=result_dir(),
    run:
        with open(output.samples_tsv, 'wt') as handle:
            for sample, num_inputs in params.num_inputs_by_sample.items():
                for i in range(num_inputs):
                    file_name = '{}_{}_simplified.tsv'.format(sample, i)
                    path = os.path.join(params.result_dir, 'simplified',
                                        file_name)
                    tabbed = '\t'.join([sample, path])
                    handle.write('{}\n'.format(tabbed))


def organize_alignment_info_input(wildcards):
    inputs = dict()
    result_path = result_dir()
    inputs['gtf_dir'] = os.path.join(result_path, 'annotation')
    inputs['samples_tsv'] = os.path.join(result_path, 'samples.tsv')
    simplified_tsvs = list()
    num_by_sample = num_inputs_by_sample()
    for sample, num_inputs in num_by_sample.items():
        for i in range(num_inputs):
            file_name = '{}_{}_simplified.tsv'.format(sample, i)
            path = os.path.join(result_path, 'simplified', file_name)
            simplified_tsvs.append(path)

    inputs['simplified_tsvs'] = simplified_tsvs
    return inputs


rule organize_alignment_info:
    input:
        unpack(organize_alignment_info_input),
    output:
        align_dir=directory(os.path.join(result_dir(), 'alignments')),
    log:
        out=os.path.join(result_dir(), 'organize_alignment_info_log.out'),
        err=os.path.join(result_dir(), 'organize_alignment_info_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts',
                            'organize_alignment_info_by_gene_and_chr.py'),
    resources:
        mem_mb=config['organize_alignment_info_mem_gb'] * 1024,
        time_hours=config['organize_alignment_info_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' --gtf-dir {input.gtf_dir}'
        ' --out-dir {output.align_dir}'
        ' --samples-tsv {input.samples_tsv}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def get_combined_or_ref_gtf():
    espresso_gtf = config.get('espresso_gtf')
    if espresso_gtf:
        return os.path.join(result_dir(), 'combined.gtf')

    gtf_name = config['gtf_name']
    return os.path.join('references', gtf_name)


def organize_gene_info_input(wildcards):
    inputs = dict()
    inputs['gtf'] = get_combined_or_ref_gtf()
    return inputs


rule organize_gene_info:
    input:
        unpack(organize_gene_info_input),
    output:
        gtf_dir=directory(os.path.join(result_dir(), 'annotation')),
    log:
        out=os.path.join(result_dir(), 'organize_gene_info_log.out'),
        err=os.path.join(result_dir(), 'organize_gene_info_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts', 'organize_gene_info_by_chr.py'),
    resources:
        mem_mb=config['organize_gene_info_mem_gb'] * 1024,
        time_hours=config['organize_gene_info_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' --gtf {input.gtf}'
        ' --out-dir {output.gtf_dir}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def simplify_alignment_info_input(wildcards):
    inputs = dict()
    sample = wildcards.sample
    sample_i = int(wildcards.sample_i)
    config_keys = ['group_1_samples', 'group_2_samples']
    for config_key in config_keys:
        by_name = config.get(config_key)
        details = by_name.get(sample)
        if not details:
            continue

        sample_file_details = details[sample_i]
        if len(sample_file_details) != 1:
            raise Exception(
                'Each entry under a sample name should be one file: {}'
                .format(sample_file_details))

        sam_path = sample_file_details.get('sam')
        if sam_path:
            inputs['in_file'] = sam_path
            break

        bam_path = sample_file_details.get('bam')
        if bam_path:
            inputs['in_file'] = bam_path
            break

        raise Exception('Expecting a sam or bam file for: {}'.format(
            sample_file_details))

    return inputs


rule simplify_alignment_info:
    input:
        unpack(simplify_alignment_info_input),
    output:
        tsv=os.path.join(result_dir(), 'simplified',
                         '{sample}_{sample_i}_simplified.tsv'),
    log:
        out=os.path.join(result_dir(), 'simplified',
                         'simplify_alignment_info_{sample}_{sample_i}_log.out'),
        err=os.path.join(result_dir(), 'simplified',
                         'simplify_alignment_info_{sample}_{sample_i}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts', 'simplify_alignment_info.py'),
    resources:
        mem_mb=config['simplify_alignment_info_mem_gb'] * 1024,
        time_hours=config['simplify_alignment_info_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' --in-file {input.in_file}'
        ' --out-tsv {output.tsv}'
        ' 1> {log.out}'
        ' 2> {log.err}'


def group_names_param(wildcards):
    config_key = 'group_{}_samples'.format(wildcards.group_num)
    by_name = config.get(config_key)
    if not by_name:
        return list()

    return list(by_name.keys())


localrules: write_group_file
rule write_group_file:
    output:
        group=os.path.join(result_dir(), 'group_{group_num}.txt'),
    wildcard_constraints:
        group_num='1|2'
    params:
        group_names=group_names_param,
    run:
        comma_names = ','.join(params.group_names)
        with open(output.group, 'wt') as handle:
            handle.write('{}\n'.format(comma_names))


def rmats_long_out_dir(wildcards, output):
    return os.path.dirname(output.coeff)


def rmats_long_diff_iso_filt_param(wildcards, output):
    out_dir = rmats_long_out_dir(wildcards, output)
    return os.path.join(out_dir, 'differential_isoforms_filtered.tsv')


def adj_pvalue_param():
    adj_pvalue = config.get('adj_pvalue')
    if adj_pvalue in [None, '']:
        return ''

    return '--adj-pvalue {}'.format(adj_pvalue)


def delta_proportion_param():
    delta = config.get('delta_proportion')
    if delta in [None, '']:
        return ''

    return '--delta-proportion {}'.format(delta)


def min_cpm_per_asm_param():
    min_cpm = config.get('min_cpm_per_asm')
    if min_cpm in [None, '']:
        return ''

    return '--min-cpm-per-asm {}'.format(min_cpm)


def covar_tsv_param():
    covar_tsv = config.get('covar_tsv')
    if not covar_tsv:
        return ''

    return '--covar-tsv {}'.format(covar_tsv)


def min_isoform_reads_param():
    min_reads = config.get('min_isoform_reads')
    if min_reads in [None, '']:
        return ''

    return '--min-isoform-reads {}'.format(min_reads)


def limit_asm_to_top_n_isoforms_param():
    top_n_isoforms = config.get('limit_asm_to_top_n_isoforms')
    if not top_n_isoforms:
        return ''

    return '--limit-asm-to-top-n-isoforms {}'.format(top_n_isoforms)


def average_reads_per_group_param():
    avg_reads = config.get('average_reads_per_group')
    if avg_reads in [None, '']:
        return ''

    return '--average-reads-per-group {}'.format(avg_reads)


def average_cpm_per_group_param():
    avg_cpm = config.get('average_cpm_per_group')
    if not avg_cpm:
        return ''

    return '--average-cpm-per-group {}'.format(avg_cpm)


def min_cpm_per_group_param():
    min_cpm = config.get('min_cpm_per_group')
    if not min_cpm:
        return ''

    return '--min-cpm-per-group {}'.format(min_cpm)


def asm_proportion_of_gene_param():
    asm_prop = config.get('asm_proportion_of_gene')
    if not asm_prop:
        return ''

    return '--asm-proportion-of-gene {}'.format(asm_prop)


def rmats_long_input(wildcards):
    inputs = dict()
    inputs['align_dir'] = os.path.join(result_dir(), 'alignments')
    counts = '{}_counts.tsv'.format(wildcards.mode)
    inputs['asm_counts'] = os.path.join(result_dir(), counts)
    event_dir = '{}_events'.format(wildcards.mode)
    inputs['event_dir'] = os.path.join(result_dir(), event_dir)
    inputs['group_1'] = os.path.join(result_dir(), 'group_1.txt')
    inputs['group_2'] = os.path.join(result_dir(), 'group_2.txt')
    return inputs


# diff_iso_filt_cp is used because rmats_long_plots will overwrite the
# original file. plot_gene_sets can safely use the copy.
rule rmats_long:
    input:
        unpack(rmats_long_input),
    output:
        coeff=os.path.join(result_dir(), 'rmats_long_{mode}', 'coeff.tsv'),
        count_tsv=os.path.join(result_dir(), 'rmats_long_{mode}', 'count.tsv'),
        diff_asm=os.path.join(result_dir(), 'rmats_long_{mode}',
                              'differential_asms.tsv'),
        diff_iso=os.path.join(result_dir(), 'rmats_long_{mode}',
                              'differential_isoforms.tsv'),
        diff_iso_filt_cp=os.path.join(
            result_dir(), 'rmats_long_{mode}',
            'differential_isoforms_filtered_copy.tsv'),
        error=os.path.join(result_dir(), 'rmats_long_{mode}', 'error.tsv'),
        lrtp=os.path.join(result_dir(), 'rmats_long_{mode}', 'lrtp.tsv'),
        warn=os.path.join(result_dir(), 'rmats_long_{mode}', 'warning.tsv'),
    log:
        out=os.path.join(result_dir(), 'rmats_long_{mode}_log.out'),
        err=os.path.join(result_dir(), 'rmats_long_{mode}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts', 'rmats_long.py'),
        out_dir=rmats_long_out_dir,
        adj_pvalue=adj_pvalue_param(),
        delta_proportion=delta_proportion_param(),
        min_cpm_per_asm=min_cpm_per_asm_param(),
        covar_tsv=covar_tsv_param(),
        min_isoform_reads=min_isoform_reads_param(),
        top_n_isoforms=limit_asm_to_top_n_isoforms_param(),
        avg_reads_per_group=average_reads_per_group_param(),
        avg_cpm_per_group=average_cpm_per_group_param(),
        min_cpm_per_group=min_cpm_per_group_param(),
        asm_prop_of_gene=asm_proportion_of_gene_param(),
        diff_iso_filt=rmats_long_diff_iso_filt_param,
    threads: config['rmats_long_threads']
    resources:
        mem_mb=config['rmats_long_mem_gb'] * 1024,
        time_hours=config['rmats_long_time_hr'],
    shell:
        '{params.conda_wrapper} python'
        ' {params.script}'
        ' --align-dir {input.align_dir}'
        ' --asm-counts {input.asm_counts}'
        ' --event-dir {input.event_dir}'
        ' --group-1 {input.group_1}'
        ' --group-2 {input.group_2}'
        ' --out-dir {params.out_dir}'
        ' --process-top-n 0'  # instead make plots in rmats_long_plots
        ' {params.adj_pvalue}'
        ' {params.delta_proportion}'
        ' {params.min_cpm_per_asm}'
        ' {params.covar_tsv}'
        ' {params.min_isoform_reads}'
        ' {params.top_n_isoforms}'
        ' {params.avg_reads_per_group}'
        ' {params.avg_cpm_per_group}'
        ' {params.min_cpm_per_group}'
        ' {params.asm_prop_of_gene}'
        ' --num-threads {threads}'
        ' 1> {log.out}'
        ' 2> {log.err}'
        ' && cp {params.diff_iso_filt} {output.diff_iso_filt_cp}'
        ' 1>> {log.out}'
        ' 2>> {log.err}'


def rmats_long_plots_out_dir(wildcards):
    return os.path.join(result_dir(), 'rmats_long_{}'.format(wildcards.mode))


def group_1_name_param():
    name = config.get('group_1_name')
    if not name:
        return ''

    return '--group-1-name {}'.format(name)


def group_2_name_param():
    name = config.get('group_2_name')
    if not name:
        return ''

    return '--group-2-name {}'.format(name)


def plot_file_type_param():
    plot_file_type = config.get('plot_file_type')
    if not plot_file_type:
        return ''

    return '--plot-file-type {}'.format(plot_file_type)


def intron_scaling_param():
    intron_scaling = config.get('intron_scaling')
    if not intron_scaling:
        return ''

    return '--intron-scaling {}'.format(intron_scaling)


def max_transcripts_param():
    max_transcripts = config.get('max_transcripts')
    if not max_transcripts:
        return ''

    return '--max-transcripts {}'.format(max_transcripts)


def no_splice_graph_param(wildcards):
    if wildcards.mode == 'asm':
        return ''

    return '--no-splice-graph-plot'


def rmats_long_plots_summary_path_param(wildcards):
    rmats_long_dir_name = 'rmats_long_{}'.format(wildcards.mode)
    return os.path.join(result_dir(), rmats_long_dir_name, 'summary.txt')


def rmats_long_plots_input(wildcards):
    inputs = dict()
    inputs['gtf_dir'] = os.path.join(result_dir(), 'annotation')
    inputs['align_dir'] = os.path.join(result_dir(), 'alignments')
    event_dir = '{}_events'.format(wildcards.mode)
    inputs['event_dir'] = os.path.join(result_dir(), event_dir)
    inputs['gtf'] = get_combined_or_ref_gtf()
    inputs['group_1'] = os.path.join(result_dir(), 'group_1.txt')
    inputs['group_2'] = os.path.join(result_dir(), 'group_2.txt')
    rmats_long_dir = 'rmats_long_{}'.format(wildcards.mode)
    inputs['diff_iso'] = os.path.join(result_dir(), rmats_long_dir,
                                      'differential_isoforms.tsv')
    return inputs


rule rmats_long_plots:
    input:
        unpack(rmats_long_plots_input),
    output:
        by_gene=directory(os.path.join(result_dir(), 'rmats_long_{mode}',
                                       'results_by_gene')),
    log:
        out=os.path.join(result_dir(), 'rmats_long_plots_{mode}_log.out'),
        err=os.path.join(result_dir(), 'rmats_long_plots_{mode}_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts', 'rmats_long.py'),
        out_dir=rmats_long_plots_out_dir,
        group_1_name=group_1_name_param(),
        group_2_name=group_2_name_param(),
        adj_pvalue=adj_pvalue_param(),
        delta_proportion=delta_proportion_param(),
        plot_file_type=plot_file_type_param(),
        intron_scaling=intron_scaling_param(),
        max_transcripts=max_transcripts_param(),
        no_splice_graph=no_splice_graph_param,
        summary_txt=rmats_long_plots_summary_path_param,
        avg_reads_per_group=average_reads_per_group_param(),
        avg_cpm_per_group=average_cpm_per_group_param(),
        min_cpm_per_group=min_cpm_per_group_param(),
        asm_prop_of_gene=asm_proportion_of_gene_param(),
    threads: config['rmats_long_plots_threads']
    resources:
        mem_mb=config['rmats_long_plots_mem_gb'] * 1024,
        time_hours=config['rmats_long_plots_time_hr'],
    shell:
        # The summary from the non-plots run will be replaced
        'rm -f {params.summary_txt}'
        ' && {params.conda_wrapper} python'
        ' {params.script}'
        ' --diff-transcripts {input.diff_iso}'
        ' --gtf-dir {input.gtf_dir}'
        ' --align-dir {input.align_dir}'
        ' --event-dir {input.event_dir}'
        ' --gencode-gtf {input.gtf}'
        ' --group-1 {input.group_1}'
        ' --group-2 {input.group_2}'
        ' --out-dir {params.out_dir}'
        ' {params.group_1_name}'
        ' {params.group_2_name}'
        ' {params.adj_pvalue}'
        ' {params.delta_proportion}'
        ' {params.plot_file_type}'
        ' {params.intron_scaling}'
        ' {params.max_transcripts}'
        ' {params.no_splice_graph}'
        ' {params.avg_reads_per_group}'
        ' {params.avg_cpm_per_group}'
        ' {params.min_cpm_per_group}'
        ' {params.asm_prop_of_gene}'
        ' --num-threads {threads}'
        ' 1> {log.out}'
        ' 2> {log.err}'

# The differential_isoforms_filtered_copy.tsv files won't be overwritten
# by rmats_long_plots.
rule plot_gene_sets:
    input:
        full_genes=os.path.join(result_dir(), 'rmats_long_full_length',
                                'differential_isoforms_filtered_copy.tsv'),
        asm_genes=os.path.join(result_dir(), 'rmats_long_asm',
                               'differential_isoforms_filtered_copy.tsv'),
    output:
        venn=os.path.join(result_dir(), 'significant_genes_venn.png'),
    log:
        out=os.path.join(result_dir(), 'plot_gene_sets_log.out'),
        err=os.path.join(result_dir(), 'plot_gene_sets_log.err'),
    params:
        conda_wrapper=config['conda_wrapper'],
        script=os.path.join('scripts', 'plot_gene_sets.R'),
        set_1_name='Full Length',
        set_2_name='ASM',
    resources:
        mem_mb=config['plot_gene_sets_mem_gb'] * 1024,
        time_hours=config['plot_gene_sets_time_hr'],
    shell:
        '{params.conda_wrapper} Rscript'
        ' {params.script}'
        ' {input.full_genes}'
        ' {input.asm_genes}'
        ' "{params.set_1_name}"'
        ' "{params.set_2_name}"'
        ' {output.venn}'
        ' 1> {log.out}'
        ' 2> {log.err}'
