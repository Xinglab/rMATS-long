import argparse
import datetime
import math
import multiprocessing
import os
import os.path
import subprocess
import sys
import tempfile

import threadpoolctl

import rpy2
import rpy2.rinterface_lib
import rpy2.rinterface_lib.embedded

import rmats_long_utils


# `import rpy2.robjects` initializes an instance of R.
# The import is not done at the top of the file so that each
# process can initialize its own R instance.
def initialize_rpy2_instance():
    import rpy2.robjects
    import rpy2.robjects.packages
    # The mblogit() function uses BLAS library calls which may try
    # to use all available CPUs. Since run_stat_model.py already creates
    # threads, limit each BLAS call to only 1 thread.
    threadpoolctl.threadpool_limits(limits=1)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run a statistical model on isoform counts')
    compat_group_str = 'Read compatibility'
    count_group_str = 'Transcript abundance'
    use_either_str = 'Use either {} or {}'
    compat_group = parser.add_argument_group(
        compat_group_str,
        use_either_str.format(compat_group_str, count_group_str))
    count_group = parser.add_argument_group(
        count_group_str,
        use_either_str.format(count_group_str, compat_group_str))
    required_group = parser.add_argument_group('Required')
    optional_group = parser.add_argument_group('Optional')
    # TODO change name (--compat-tsv)
    compat_group.add_argument('--counts-tsv',
                              help='The input read isoform compatibility')
    count_group.add_argument(
        '--abundance', help='The path to the transcript abundance by sample')
    required_group.add_argument('--out-dir',
                                required=True,
                                help='The directory to write output to')
    required_group.add_argument(
        '--group-1',
        required=True,
        help=('The path to a file listing the sample names for group 1. The'
              ' file should have a single line with the sample names as a'
              ' comma separated list.'))
    required_group.add_argument(
        '--group-2',
        required=True,
        help='The path to a file listing the sample names for group 2')
    optional_group.add_argument(
        '--num-threads',
        type=int,
        default=1,
        help='How many threads to use (default %(default)s)')
    optional_group.add_argument(
        '--group-1-name',
        default='group_1',
        help='A name for group 1 (default %(default)s)')
    optional_group.add_argument(
        '--group-2-name',
        default='group_2',
        help='A name for group 2 (default %(default)s)')
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
              ' (default: %(default)s)'))
    optional_group.add_argument(
        '--min-cpm-per-asm',
        type=float,
        default=0,
        help=('Only consider ASMs where at least 1 sample has at least'
              ' this CPM of reads assigned to the ASM.'
              ' (default %(default)s)'))
    optional_group.add_argument(
        '--sample-read-total-tsv',
        required=True,
        help=('A .tsv file with two columns: sample and total.'
              ' The 1st line is the header'))
    optional_group.add_argument(
        '--em-tolerance',
        type=float,
        default=0.001,
        help=('Stop performing EM iterations when the change is at or below'
              ' this value (default %(default)s)'))
    optional_group.add_argument(
        '--em-max-iter',
        type=int,
        default=100,
        help='Perform at most this many EM iterations (default %(default)s)')
    optional_group.add_argument(
        '--random-seed',
        type=int,
        default=123,
        help='Passed to R base::set.seed() (default %(default)s)')
    optional_group.add_argument(
        '--progress-every-n',
        type=int,
        default=100,
        help=('Print a status message after a certain number of ASMs'
              ' (default %(default)s)'))
    optional_group.add_argument(
        '--sort-buffer-size',
        default='2G',
        help=('Used for the --buffer-size argument of sort.'
              ' Default: %(default)s'))
    optional_group.add_argument(
        '--covar-tsv',
        help=('A .tsv with 1 line per sample. The first line has the column'
              ' names. The first column is sample_id. Each additional column'
              ' is a covariate.'))

    args = parser.parse_args()
    has_counts = args.counts_tsv is not None
    has_abun = args.abundance is not None
    if (has_counts and has_abun) or not (has_counts or has_abun):
        parser.error('Exactly one of --counts-tsv or --abundance is required')

    return args


def import_r_packages():
    packages = dict()
    packages['mclogit'] = rpy2.robjects.packages.importr('mclogit')
    packages['base'] = rpy2.robjects.packages.importr('base')
    packages['stats'] = rpy2.robjects.packages.importr('stats')
    return packages


def create_append_to_list_func(target):
    def append_to(x):
        target.append(x)

    return append_to


def capture_r_output():
    result = dict()
    result['prev_out_callback'] = (
        rpy2.rinterface_lib.callbacks.consolewrite_print)
    result['prev_err_callback'] = (
        rpy2.rinterface_lib.callbacks.consolewrite_warnerror)
    out_list = list()
    err_list = list()
    result['out'] = out_list
    result['err'] = err_list
    rpy2.rinterface_lib.callbacks.consolewrite_print = (
        create_append_to_list_func(out_list))
    rpy2.rinterface_lib.callbacks.consolewrite_warnerror = (
        create_append_to_list_func(err_list))

    return result


def restore_r_output_callbacks(result):
    rpy2.rinterface_lib.callbacks.consolewrite_print = (
        result['prev_out_callback'])
    rpy2.rinterface_lib.callbacks.consolewrite_warnerror = (
        result['prev_err_callback'])


def print_with_timestamp(message, flush=True):
    format_str = '%Y-%m-%dT%H:%M:%S'
    now = datetime.datetime.now()
    time_str = now.strftime(format_str)
    formatted = '[{}] {}'.format(time_str, message)
    print(formatted, flush=flush)


def parse_group_file(path):
    groups = list()
    with open(path, 'rt') as handle:
        for line_i, line in enumerate(handle):
            if line_i != 0 and line != '':
                raise Exception(
                    '{} should only have groups on the 1st line'.format(path))

            groups = line.rstrip('\n').split(',')

    return groups


def determine_covar_type(values):
    try:
        for value in values:
            int(value)

        return 'int'
    except ValueError:
        pass

    try:
        for value in values:
            float(value)

        return 'float'
    except ValueError:
        pass

    return 'string'


def try_to_convert_covar_strings_to_int_or_float(names, by_sample):
    for covar in names:
        covar_values = list()
        for sample_dict in by_sample.values():
            value = sample_dict[covar]
            covar_values.append(value)

        covar_type = determine_covar_type(covar_values)
        if covar_type == 'string':
            continue

        if covar_type == 'int':
            for sample_dict in by_sample.values():
                sample_dict[covar] = int(sample_dict[covar])

        if covar_type == 'float':
            for sample_dict in by_sample.values():
                sample_dict[covar] = float(sample_dict[covar])


def parse_covar_tsv(path):
    by_sample = dict()
    names = list()
    covars = {'by_sample': by_sample, 'names': names}
    if not path:
        return covars

    sample_header = 'sample_id'
    with open(path, 'rt') as handle:
        for line_i, line in enumerate(handle):
            columns = rmats_long_utils.read_tsv_line(line)
            if line_i == 0:
                headers = columns
                if sample_header not in headers:
                    raise Exception('Expected column name {} in {}'.format(
                        sample_header, path))

                # TODO add a prefix or suffix to avoid any conflict with
                # other columns (sample_id, isoform_id, count, group)
                other_headers = [x for x in headers if x != sample_header]
                names.extend(other_headers)
                continue

            row = dict(zip(headers, columns))
            sample = row.pop(sample_header)
            by_sample[sample] = row

    try_to_convert_covar_strings_to_int_or_float(names, by_sample)
    return covars


def convert_covar_to_vector(values):
    if len(values) == 0:
        return rpy2.robjects.StrVector(values)

    value = values[0]
    if isinstance(value, int):
        return rpy2.robjects.IntVector(values)

    if isinstance(value, float):
        return rpy2.robjects.FloatVector(values)

    return rpy2.robjects.StrVector(values)


def write_output_headers(out_handles):
    rmats_long_utils.write_tsv_line(
        out_handles['count'],
        ['asm_id', 'gene_id', 'sample_id', 'isoform_id', 'count', 'prop'])
    rmats_long_utils.write_tsv_line(
        out_handles['coeff'],
        ['asm_id', 'gene_id', 'isoform_id', 'lr', 'df', 'p_value'])
    rmats_long_utils.write_tsv_line(
        out_handles['lrtp'],
        ['asm_id', 'gene_id', 'lrt_p', 'lr', 'df', 'converge', 'iter', 'eps'])
    rmats_long_utils.write_tsv_line(out_handles['error'],
                                    ['asm_id', 'gene_id', 'error'])
    rmats_long_utils.write_tsv_line(out_handles['warn'],
                                    ['asm_id', 'gene_id', 'warning'])


def get_isoforms_for_read(rows, all_isoforms):
    sample = None
    isoforms = list()
    for row in rows:
        sample = row['sample_id']
        isoform = row['isoform_id']
        isoforms.append(isoform)
        all_isoforms.add(isoform)

    return {'sample_id': sample, 'isoforms': isoforms}


def update_single_isoform(isoform, sample, single_isoform):
    for_sample = single_isoform.get(sample)
    if not for_sample:
        for_sample = dict()
        single_isoform[sample] = for_sample

    old_count = for_sample.get(isoform, 0)
    for_sample[isoform] = old_count + 1


def update_multi_isoform(isoforms, sample, multi_isoform):
    for_sample = multi_isoform.get(sample)
    if not for_sample:
        for_sample = dict()
        multi_isoform[sample] = for_sample

    isoform_key = tuple(sorted(isoforms))
    old_count = for_sample.get(isoform_key, 0)
    for_sample[isoform_key] = old_count + 1


def create_read_counts_initial_grouping():
    return {
        'isoforms': set(),
        'single': dict(),
        'multi': dict(),
    }


# TODO maybe use a mapping of isoforms -> int
def update_read_counts_grouping(read_rows, read_grouping):
    if not read_rows:
        return

    all_isoforms = read_grouping['isoforms']
    read_result = get_isoforms_for_read(read_rows, all_isoforms)
    sample = read_result['sample_id']
    isoforms = read_result['isoforms']

    single_isoform = read_grouping['single']
    multi_isoform = read_grouping['multi']
    if len(isoforms) == 1:
        update_single_isoform(isoforms[0], sample, single_isoform)
    else:
        update_multi_isoform(isoforms, sample, multi_isoform)


def get_min_samples_per_group(sample_to_group):
    samples_by_group = dict()
    for group in sample_to_group.values():
        old_count = samples_by_group.get(group, 0)
        samples_by_group[group] = old_count + 1

    return min(samples_by_group.values())


def create_grouping_from_isoform_counts(by_sample_by_isoform, sample_to_group):
    min_samples_per_group = get_min_samples_per_group(sample_to_group)
    single_isoform = dict()
    multi_isoform = dict()
    num_samples_with_abun_by_isoform = dict()
    for sample, by_isoform in by_sample_by_isoform.items():
        for_sample = dict()
        single_isoform[sample] = for_sample
        for isoform, abundance in by_isoform.items():
            for_sample[isoform] = abundance
            if abundance >= 1:
                old_count = num_samples_with_abun_by_isoform.get(isoform, 0)
                num_samples_with_abun_by_isoform[isoform] = old_count + 1

    # Keep only isoforms with at least 1 abundance in at least
    # min_samples_per_group.
    # Based on a recommendation in the documentation for DRIMSeq::dmFilter
    kept_isoforms = set()
    for isoform, num_samples in num_samples_with_abun_by_isoform.items():
        if num_samples >= min_samples_per_group:
            kept_isoforms.add(isoform)

    for for_sample in single_isoform.values():
        to_remove = set(for_sample.keys()).difference(kept_isoforms)
        for isoform in to_remove:
            del for_sample[isoform]

    return {
        'isoforms': kept_isoforms,
        'single': single_isoform,
        'multi': multi_isoform
    }


def intersect_with_column_names(values, data_frame, r_packages):
    result = list()
    col_names = set(r_packages['base'].colnames(data_frame))
    for value in values:
        if value in col_names:
            result.append(value)

    return result


def create_formulas(var, covars, data_frame, response, r_packages):
    var_vector = rpy2.robjects.StrVector(var)
    covar_names = covars['names']
    used_covars = intersect_with_column_names(covar_names, data_frame,
                                              r_packages)
    if used_covars:
        covar_vector = rpy2.robjects.StrVector(used_covars)
        full = r_packages['stats'].reformulate(var_vector + covar_vector,
                                               response=response)
        null = r_packages['stats'].reformulate(covar_vector, response=response)
    else:
        full = r_packages['stats'].reformulate(var_vector, response=response)
        null = rpy2.robjects.r('{} ~ 1'.format(response))

    return {'full': full, 'null': null}


def initialize_single_isoform_totals(single_isoform, total_single_by_isoform,
                                     counts_by_sample):
    for sample, by_isoform in single_isoform.items():
        count_by_isoform = rmats_long_utils.try_get_or_set_default(
            counts_by_sample, sample, dict())
        for isoform, count in by_isoform.items():
            old_total = total_single_by_isoform.get(isoform, 0)
            total_single_by_isoform[isoform] = old_total + count
            old_count = count_by_isoform.get(isoform, 0)
            count_by_isoform[isoform] = old_count + count


def initialize_multi_isoform_totals(multi_isoform, total_multi_by_isoform,
                                    counts_by_sample):
    for sample, by_isoform in multi_isoform.items():
        count_by_isoform = rmats_long_utils.try_get_or_set_default(
            counts_by_sample, sample, dict())
        for isoforms, count in by_isoform.items():
            num_isoforms = len(isoforms)
            count_per = count / num_isoforms
            for isoform in isoforms:
                old_total = total_multi_by_isoform.get(isoform, 0)
                total_multi_by_isoform[isoform] = old_total + count_per
                old_count = count_by_isoform.get(isoform, 0)
                count_by_isoform[isoform] = old_count + count_per


def initialize_counts(single_isoform, multi_isoform):
    total_single_by_isoform = dict()
    total_multi_by_isoform = dict()
    counts_by_sample = dict()
    initialize_single_isoform_totals(single_isoform, total_single_by_isoform,
                                     counts_by_sample)
    initialize_multi_isoform_totals(multi_isoform, total_multi_by_isoform,
                                    counts_by_sample)

    totals = {
        'total_single_by_isoform': total_single_by_isoform,
        'total_multi_by_isoform': total_multi_by_isoform
    }
    return {'counts_by_sample': counts_by_sample, 'totals': totals}


def distribute_counts_with_all_estimated(isoforms, combo_count, prop_total,
                                         prop_by_isoform,
                                         new_count_by_isoform):
    for isoform in isoforms:
        prop = prop_by_isoform[isoform]
        ratio = prop / prop_total
        count = combo_count * ratio
        old_count = new_count_by_isoform.get(isoform, 0)
        new_count_by_isoform[isoform] = old_count + count


def distribute_counts_with_some_zeros(isoforms, combo_count, count_total,
                                      count_by_isoform, new_count_by_isoform):
    # At least 1 isoform doesn't have an estimate.
    # Add 1/num_isoforms_in_combo to each isoform so that
    # all isoforms have some proportion.
    by_n = 1 / len(isoforms)
    count_total += 1
    for isoform in isoforms:
        count_for_ratio = count_by_isoform.get(isoform, 0)
        count_for_ratio += by_n
        ratio = count_for_ratio / count_total
        count = combo_count * ratio
        old_count = new_count_by_isoform.get(isoform, 0)
        new_count_by_isoform[isoform] = old_count + count


def distribute_counts_for_multis(multi_isoform, counts_by_sample,
                                 props_by_sample, new_counts_by_sample):
    for sample, by_combo in multi_isoform.items():
        count_by_isoform = counts_by_sample.get(sample, dict())
        prop_by_isoform = props_by_sample.get(sample, dict())
        new_count_by_isoform = dict()
        new_counts_by_sample[sample] = new_count_by_isoform
        for isoforms, combo_count in by_combo.items():
            count_total = 0
            prop_total = 0
            any_zero = False
            for isoform in isoforms:
                count = count_by_isoform.get(isoform, 0)
                count_total += count
                prop = prop_by_isoform.get(isoform, 0)
                prop_total += prop
                if count == 0 or prop == 0:
                    any_zero = True

            if any_zero:
                distribute_counts_with_some_zeros(isoforms, combo_count,
                                                  count_total,
                                                  count_by_isoform,
                                                  new_count_by_isoform)
            else:
                distribute_counts_with_all_estimated(isoforms, combo_count,
                                                     prop_total,
                                                     prop_by_isoform,
                                                     new_count_by_isoform)


def distribute_counts_for_singles(single_isoform, new_counts_by_sample):
    for sample, by_isoform in single_isoform.items():
        new_count_by_isoform = rmats_long_utils.try_get_or_set_default(
            new_counts_by_sample, sample, dict())
        for isoform, count in by_isoform.items():
            old_count = new_count_by_isoform.get(isoform, 0)
            new_count_by_isoform[isoform] = old_count + count


def distribute_counts(single_isoform, multi_isoform, counts_and_props):
    new_counts_by_sample = dict()
    counts_by_sample = counts_and_props['count_by_sample']
    props_by_sample = counts_and_props['prop_by_sample']

    distribute_counts_for_multis(multi_isoform, counts_by_sample,
                                 props_by_sample, new_counts_by_sample)
    distribute_counts_for_singles(single_isoform, new_counts_by_sample)
    counts_and_props['count_by_sample'] = new_counts_by_sample


def select_top_n_isoforms_by_proportion(props_by_sample,
                                        limit_asm_to_top_n_isoforms):
    total_by_isoform = dict()
    for prop_by_isoform in props_by_sample.values():
        for isoform, prop in prop_by_isoform.items():
            old_isoform_total = total_by_isoform.get(isoform, 0)
            total_by_isoform[isoform] = old_isoform_total + prop

    sorted_isoforms_and_props = sorted(total_by_isoform.items(),
                                       key=lambda pair: pair[1],
                                       reverse=True)
    top_isoforms = set()
    n = 0
    nth_highest_prop = None
    for isoform, prop in sorted_isoforms_and_props:
        if not limit_asm_to_top_n_isoforms:
            top_isoforms.add(isoform)
            continue

        if n < limit_asm_to_top_n_isoforms:
            n += 1
            top_isoforms.add(isoform)
            nth_highest_prop = prop
            continue

        # Keep all isoforms that have the nth highest prop
        if prop == nth_highest_prop:
            top_isoforms.add(isoform)
        else:
            break

    return top_isoforms


def create_model_data_frame(min_count, limit_asm_to_top_n_isoforms,
                            sample_to_group, covars, isoform_order,
                            counts_and_props, r_packages):
    samples = list()
    isoforms = list()
    counts = list()
    groups = list()
    covar_vectors = dict()
    covar_names = covars['names']
    covar_by_sample = covars['by_sample']
    for name in covar_names:
        covar_vectors[name] = list()

    filtered_by_top_n = set()
    props_by_sample = counts_and_props['prop_by_sample']
    top_n_isoforms = select_top_n_isoforms_by_proportion(
        props_by_sample, limit_asm_to_top_n_isoforms)
    counts_by_sample = counts_and_props['count_by_sample']
    for sample, count_by_isoform in counts_by_sample.items():
        group = sample_to_group[sample]
        sample_covars = covar_by_sample.get(sample)
        for isoform, count in count_by_isoform.items():
            if (count == 0) or (count < min_count):
                continue

            if isoform not in top_n_isoforms:
                filtered_by_top_n.add(isoform)
                continue

            samples.append(sample)
            isoforms.append(isoform)
            counts.append(count)
            groups.append(group)
            if not sample_covars:
                continue

            for name, value in sample_covars.items():
                covar_vectors[name].append(value)

    warnings = list()
    if filtered_by_top_n:
        warning = '{} isoforms filtered by --limit-asm-to-top-n-isoforms'.format(
            len(filtered_by_top_n))
        warnings.append(warning)

    isoform_set = set(isoforms)
    ref_isoform = None
    for isoform in isoform_order:
        if isoform in isoform_set:
            ref_isoform = isoform
            break

    if ref_isoform is None:
        return {'error': 'No isoforms after count filter', 'warn': warnings}

    counts_vector = rpy2.robjects.FloatVector(counts)
    samples_vector = rpy2.robjects.StrVector(samples)
    isoforms_vector = rpy2.robjects.StrVector(isoforms).factor()
    # Set which isoform will be used as the reference in the model
    isoforms_vector = r_packages['stats'].relevel(isoforms_vector,
                                                  ref=ref_isoform)
    isoform_levels = r_packages['base'].levels(isoforms_vector)
    groups_vector = rpy2.robjects.StrVector(groups)
    data_dict = {
        'sample_id': samples_vector,
        'isoform_id': isoforms_vector,
        'count': counts_vector,
        'group': groups_vector
    }
    for name, values in covar_vectors.items():
        num_unique = len(set(values))
        if num_unique >= 2:
            vector = convert_covar_to_vector(values)
            data_dict[name] = vector

    data = rpy2.robjects.DataFrame(data_dict)

    return {
        'error': None,
        'warn': warnings,
        'data': data,
        'counts': counts_vector,
        'samples': samples,
        'isoforms': isoforms,
        'isoform_levels': isoform_levels,
        'groups': groups
    }


def run_mblogit(isoform_levels,
                groups,
                formula,
                data,
                weights,
                r_packages,
                dispersion=True):
    num_isoforms = len(isoform_levels)
    if num_isoforms < 2:
        return {'error': 'Less than 2 isoforms: {}'.format(num_isoforms)}

    num_groups = len(set(groups))
    if num_groups < 2:
        return {'error': 'Less than 2 groups: {}'.format(num_groups)}

    # suppress printing from mblogit
    capture_details = capture_r_output()
    try:
        fit = r_packages['mclogit'].mblogit(formula,
                                            data=data,
                                            weights=weights,
                                            dispersion=dispersion)
    except rpy2.rinterface_lib.embedded.RRuntimeError as e:
        return {'error': str(e)}
    finally:
        restore_r_output_callbacks(capture_details)

    return {'error': None, 'fit': fit}


def fit_model(min_count, limit_asm_to_top_n_isoforms, sample_to_group, var,
              covars, isoform_order, counts_and_props, r_packages):
    model_df = create_model_data_frame(min_count, limit_asm_to_top_n_isoforms,
                                       sample_to_group, covars, isoform_order,
                                       counts_and_props, r_packages)
    if model_df['error']:
        return {'error': model_df['error'], 'warn': model_df['warn']}

    data = model_df['data']
    counts = model_df['counts']
    isoform_levels = model_df['isoform_levels']
    groups = model_df['groups']
    formulas = create_formulas(var, covars, data, 'isoform_id', r_packages)
    mblogit_result = run_mblogit(isoform_levels, groups, formulas['full'],
                                 data, counts, r_packages)
    if mblogit_result['error']:
        return {'error': mblogit_result['error'], 'warn': model_df['warn']}

    fit = mblogit_result['fit']
    return {
        'error': None,
        'warn': model_df['warn'],
        'fit': fit,
        'model_df': model_df
    }


def calculate_p_value(var, covars, fit, model_df, r_packages):
    data = model_df['data']
    counts = model_df['counts']
    isoform_levels = model_df['isoform_levels']
    groups = model_df['groups']

    formulas = create_formulas(var, covars, data, 'isoform_id', r_packages)
    mblogit_result = run_mblogit(isoform_levels, groups, formulas['null'],
                                 data, counts, r_packages)
    if mblogit_result['error']:
        return {'error': mblogit_result['error']}

    null_fit = mblogit_result['fit']
    stats = get_stats_from_full_and_null_fits(fit, null_fit, r_packages)
    return {
        'error': None,
        'p_value': stats['p_value'],
        'lr': stats['lr'],
        'df': stats['df']
    }


def calculate_props(counts_by_sample, props_by_sample):
    max_delta = 0
    for sample, count_by_isoform in counts_by_sample.items():
        sample_total = sum(count_by_isoform.values())
        prop_by_isoform = dict()
        old_prop_by_isoform = props_by_sample.get(sample, dict())
        props_by_sample[sample] = prop_by_isoform
        if sample_total == 0:
            continue

        for isoform, count in count_by_isoform.items():
            prop = count / sample_total
            prop_by_isoform[isoform] = prop
            old_prop = old_prop_by_isoform.get(isoform)
            if old_prop is not None:
                delta = abs(prop - old_prop)
                max_delta = max(max_delta, delta)

    return max_delta


def init_props(counts_by_sample):
    props_by_sample = dict()
    calculate_props(counts_by_sample, props_by_sample)
    return props_by_sample


def update_props(counts_and_props):
    counts_by_sample = counts_and_props['count_by_sample']
    props_by_sample = counts_and_props['prop_by_sample']
    delta = calculate_props(counts_by_sample, props_by_sample)
    return delta


def em_iteration(single_isoform, multi_isoform, counts_and_props):
    distribute_counts(single_isoform, multi_isoform, counts_and_props)
    delta = update_props(counts_and_props)
    return delta


def get_total_by_sample_from_model_counts(model_counts, model_samples):
    total_by_sample = dict()
    for count_i, count in enumerate(model_counts):
        sample = model_samples[count_i]
        old_total = total_by_sample.get(sample, 0)
        total_by_sample[sample] = old_total + count

    return total_by_sample


def get_or_initialize_isoform_details(covar_names, details_by_isoform,
                                      isoform):
    details = details_by_isoform.get(isoform)
    if details:
        return details

    details = {
        'samples': list(),
        'groups': list(),
        'counts': list(),
        'is_isoform': list()
    }
    covar_vectors = dict()
    for name in covar_names:
        covar_vectors[name] = list()

    details['covars'] = covar_vectors
    details_by_isoform[isoform] = details
    return details


def get_or_set_sample_default(default_by_sample, sample, group, total):
    sample_default = default_by_sample.get(sample)
    if sample_default:
        return sample_default

    sample_default = {
        'samples': [sample],
        'groups': [group],
        'counts': [total],
        'is_isoform': ['n'],
        'covars': dict()
    }
    default_by_sample[sample] = sample_default
    return sample_default


def get_isoform_details_and_set_counts(covar_names, covar_by_sample,
                                       model_counts, model_samples,
                                       model_isoforms, model_groups,
                                       total_by_sample, details_by_isoform,
                                       counts, default_by_sample):
    for count_i, count in enumerate(model_counts):
        sample = model_samples[count_i]
        isoform = model_isoforms[count_i]
        group = model_groups[count_i]
        total = total_by_sample[sample]
        if total == 0:
            continue

        prop = count / total
        remaining_count = total - count
        counts.append({
            'sample_id': sample,
            'isoform_id': isoform,
            'count': count,
            'prop': prop
        })

        details = get_or_initialize_isoform_details(covar_names,
                                                    details_by_isoform,
                                                    isoform)
        sample_default = get_or_set_sample_default(default_by_sample, sample,
                                                   group, total)

        details['samples'].extend([sample, sample])
        details['groups'].extend([group, group])
        details['counts'].extend([count, remaining_count])
        details['is_isoform'].extend(['y', 'n'])
        covar_vectors = details['covars']
        sample_covars = covar_by_sample.get(sample)
        if not sample_covars:
            continue

        for name, value in sample_covars.items():
            covar_vectors[name].extend([value, value])
            sample_default['covars'][name] = [value]


def create_isoform_data_frames(details_by_isoform, default_by_sample):
    df_by_isoform = dict()
    for isoform, details in details_by_isoform.items():
        for sample, default in default_by_sample.items():
            if sample not in details['samples']:
                details['samples'].extend(default['samples'])
                details['groups'].extend(default['groups'])
                details['counts'].extend(default['counts'])
                details['is_isoform'].extend(default['is_isoform'])
                for name, values in details['covars'].items():
                    values.extend(default['covars'][name])

        samples_vector = rpy2.robjects.StrVector(details['samples'])
        groups_vector = rpy2.robjects.StrVector(details['groups'])
        counts_vector = rpy2.robjects.FloatVector(details['counts'])
        is_isoform_vector = rpy2.robjects.StrVector(
            details['is_isoform']).factor()
        data_dict = {
            'sample_id': samples_vector,
            'group': groups_vector,
            'count': counts_vector,
            'is_isoform': is_isoform_vector
        }
        for name, values in details['covars'].items():
            num_unique = len(set(values))
            if num_unique >= 2:
                vector = convert_covar_to_vector(values)
                data_dict[name] = vector

        data = rpy2.robjects.DataFrame(data_dict)
        df_by_isoform[isoform] = data

    return df_by_isoform


def create_isoform_data_frames_and_set_counts(model_df, counts, covars):
    model_counts = model_df['counts']
    model_samples = model_df['samples']
    model_groups = model_df['groups']
    model_isoforms = model_df['isoforms']
    covar_names = covars['names']
    covar_by_sample = covars['by_sample']
    total_by_sample = get_total_by_sample_from_model_counts(
        model_counts, model_samples)

    details_by_isoform = dict()
    default_by_sample = dict()
    get_isoform_details_and_set_counts(covar_names, covar_by_sample,
                                       model_counts, model_samples,
                                       model_isoforms, model_groups,
                                       total_by_sample, details_by_isoform,
                                       counts, default_by_sample)

    df_by_isoform = create_isoform_data_frames(details_by_isoform,
                                               default_by_sample)
    return df_by_isoform


def get_stats_from_full_and_null_fits(fit, null_fit, r_packages):
    test = r_packages['stats'].anova(null_fit, fit, test='Chisq')
    p_value = test.rx2('Pr(>Chi)').rx(2)
    p_value = p_value[0]  # convert to a float
    lr = test.rx2('Deviance').rx(2)
    lr = lr[0]
    df = test.rx2('Df').rx(2)
    df = df[0]
    return {'p_value': p_value, 'lr': lr, 'df': df}


def test_isoform_proportions(var, covars, model_df, counts, coefficients,
                             r_packages):
    df_by_isoform = create_isoform_data_frames_and_set_counts(
        model_df, counts, covars)
    for isoform, data_frame in df_by_isoform.items():
        error_result = {'isoform_id': isoform, 'lr': 0, 'df': 0, 'p_value': 1}
        formulas = create_formulas(var, covars, data_frame, 'is_isoform',
                                   r_packages)
        isoform_levels = data_frame.rx2('is_isoform')
        groups = data_frame.rx2('group')
        weights = data_frame.rx2('count')
        full_mblogit_result = run_mblogit(isoform_levels, groups,
                                          formulas['full'], data_frame,
                                          weights, r_packages)
        if full_mblogit_result['error']:
            coefficients.append(error_result)
            continue

        null_mblogit_result = run_mblogit(isoform_levels, groups,
                                          formulas['null'], data_frame,
                                          weights, r_packages)
        if null_mblogit_result['error']:
            coefficients.append(error_result)
            continue

        fit = full_mblogit_result['fit']
        null_fit = null_mblogit_result['fit']
        stats = get_stats_from_full_and_null_fits(fit, null_fit, r_packages)
        coefficients.append({
            'isoform_id': isoform,
            'lr': stats['lr'],
            'df': stats['df'],
            'p_value': stats['p_value']
        })


def run_em_iterations(tolerance, max_iter, single_isoform, multi_isoform,
                      counts_and_props, converge):
    em_i = 0
    while True:
        delta = em_iteration(single_isoform, multi_isoform, counts_and_props)
        em_i += 1
        if delta <= tolerance:
            converge['converged'] = True
            converge['delta'] = delta
            break

        if em_i == max_iter:
            converge['converged'] = False
            converge['delta'] = delta
            break

    converge['iter'] = em_i


def determine_isoform_order(totals):
    single_by_isoform = totals['total_single_by_isoform']
    multi_by_isoform = totals['total_multi_by_isoform']
    # Put isoforms with higher uniquely assigned read counts first
    isoform_order = sorted(single_by_isoform.keys(),
                           key=lambda k: single_by_isoform[k],
                           reverse=True)
    multi_order = sorted(multi_by_isoform.keys(),
                         key=lambda k: multi_by_isoform[k],
                         reverse=True)

    # Then order any isoforms without uniquely assigned reads
    from_single = set(isoform_order)
    for isoform in multi_order:
        if isoform not in from_single:
            isoform_order.append(isoform)

    return isoform_order


def check_min_asm_reads(single_isoform, multi_isoform,
                        min_asm_reads_by_sample):
    fails = list()
    samples = set(single_isoform.keys()).union(set(multi_isoform.keys()))
    samples = sorted(samples)
    for sample in samples:
        min_count = min_asm_reads_by_sample.get(sample, 0)
        sample_total = 0
        single_by_isoform = single_isoform.get(sample, dict())
        multi_by_isoform = multi_isoform.get(sample, dict())
        for count in single_by_isoform.values():
            sample_total += count

        for count in multi_by_isoform.values():
            sample_total += count

        if sample_total >= min_count:
            return None

        fails.append('{}: {}<{:.2g}'.format(sample, sample_total, min_count))

    return 'check_min_asm_reads: {}'.format(','.join(fails))


def run_stat_model_on_asm(asm_id, gene_id, read_grouping, min_count,
                          limit_asm_to_top_n_isoforms, min_asm_reads_by_sample,
                          tolerance, max_iter, seed, sample_to_group, covars,
                          r_packages):
    var = ['group']
    rpy2.robjects.r('base::set.seed({})'.format(seed))

    counts = list()
    coefficients = list()
    converge = {'iter': None, 'converged': None, 'delta': None}
    warnings = list()
    results = {
        'asm_id': asm_id,
        'gene_id': gene_id,
        'error': None,
        'warn': warnings,
        'count': counts,
        'coeff': coefficients,
        'LRT_p': {
            'p_value': None,
            'lr': None,
            'df': None
        },
        'converge': converge
    }

    isoforms = sorted(read_grouping['isoforms'])
    single_isoform = read_grouping['single']
    multi_isoform = read_grouping['multi']
    num_isoforms = len(isoforms)
    if num_isoforms < 2:
        results['error'] = 'Less than 2 isoforms: {}'.format(num_isoforms)
        return results

    min_reads_error = check_min_asm_reads(single_isoform, multi_isoform,
                                          min_asm_reads_by_sample)
    if min_reads_error:
        results['error'] = min_reads_error
        return results

    init_result = initialize_counts(single_isoform, multi_isoform)
    counts_by_sample = init_result['counts_by_sample']
    init_totals = init_result['totals']
    isoform_order = determine_isoform_order(init_totals)
    props_by_sample = init_props(counts_by_sample)
    counts_and_props = {
        'count_by_sample': counts_by_sample,
        'prop_by_sample': props_by_sample
    }
    if multi_isoform:
        run_em_iterations(tolerance, max_iter, single_isoform, multi_isoform,
                          counts_and_props, converge)

    fit_results = fit_model(min_count, limit_asm_to_top_n_isoforms,
                            sample_to_group, var, covars, isoform_order,
                            counts_and_props, r_packages)
    warnings.extend(fit_results.get('warn', list()))
    if fit_results['error']:
        results['error'] = fit_results['error']
        return results

    fit = fit_results['fit']
    model_df = fit_results['model_df']
    p_value_result = calculate_p_value(var, covars, fit, model_df, r_packages)
    if p_value_result['error']:
        results['error'] = p_value_result['error']
        return results

    p_value = p_value_result['p_value']
    lr = p_value_result['lr']
    df = p_value_result['df']
    results['LRT_p']['p_value'] = p_value
    results['LRT_p']['lr'] = lr
    results['LRT_p']['df'] = df

    test_isoform_proportions(var, covars, model_df, counts, coefficients,
                             r_packages)
    return results


def parse_float(value):
    if value == 'NA':
        return None

    return float(value)


def format_bool(value):
    if value is None:
        return 'NA'

    if value:
        return 'true'

    return 'false'


def escape_newlines(string):
    return string.replace('\n', '\\n')


def write_results(results, out_handles):
    asm_id = results['asm_id']
    gene_id = results['gene_id']
    if results.get('warn'):
        formatted = ';'.join(results['warn'])
        rmats_long_utils.write_tsv_line(out_handles['warn'],
                                        [asm_id, gene_id, formatted])

    if results['error']:
        escaped = escape_newlines(results['error'])
        rmats_long_utils.write_tsv_line(out_handles['error'],
                                        [asm_id, gene_id, escaped])
        return

    count_handle = out_handles['count']
    for count_row in results['count']:
        rmats_long_utils.write_tsv_line(count_handle, [
            asm_id, gene_id, count_row['sample_id'], count_row['isoform_id'],
            rmats_long_utils.format_float(count_row['count']),
            rmats_long_utils.format_float(count_row['prop'])
        ])

    coeff_handle = out_handles['coeff']
    for coeff_row in results['coeff']:
        rmats_long_utils.write_tsv_line(coeff_handle, [
            asm_id, gene_id, coeff_row['isoform_id'],
            rmats_long_utils.format_float(coeff_row['lr']),
            rmats_long_utils.format_float(coeff_row['df']),
            rmats_long_utils.format_float(coeff_row['p_value'])
        ])

    lrtp_handle = out_handles['lrtp']
    conv_results = results['converge']
    lrtp_results = results['LRT_p']
    p_value = lrtp_results['p_value']
    lr = lrtp_results['lr']
    df = lrtp_results['df']
    rmats_long_utils.write_tsv_line(lrtp_handle, [
        asm_id, gene_id,
        rmats_long_utils.format_float(p_value),
        rmats_long_utils.format_float(lr),
        rmats_long_utils.format_float(df),
        format_bool(conv_results['converged']),
        rmats_long_utils.format_float(conv_results['iter']),
        rmats_long_utils.format_float(conv_results['delta'])
    ])


def make_sample_to_group(group_1, group_2, group_1_name, group_2_name):
    sample_to_group = dict()
    for sample in group_1:
        sample_to_group[sample] = group_1_name

    for sample in group_2:
        sample_to_group[sample] = group_2_name

    return sample_to_group


def create_threads_and_queues(num_threads, min_count,
                              limit_asm_to_top_n_isoforms,
                              min_asm_reads_by_sample, tolerance, max_iter,
                              seed, sample_to_group, covars):
    # Every thread should be able to have something on the queue.
    # queue_size_mult allows some extra room.
    queue_size_mult = 10
    threads = list()
    result = {
        'threads': threads,
        'in_queue': None,
        'out_queue': None,
        'signal_queue': None
    }
    if num_threads == 1:
        return result

    num_workers = num_threads - 1
    signal_queue = multiprocessing.Queue(num_workers)
    queue_size = num_workers * queue_size_mult
    in_queue = multiprocessing.Queue(queue_size)
    out_queue = multiprocessing.Queue(queue_size)
    result['in_queue'] = in_queue
    result['out_queue'] = out_queue
    result['signal_queue'] = signal_queue
    for _ in range(num_workers):
        thread = multiprocessing.Process(
            target=run_stat_model_thread,
            args=(min_count, limit_asm_to_top_n_isoforms,
                  min_asm_reads_by_sample, tolerance, max_iter, seed,
                  sample_to_group, covars, in_queue, out_queue, signal_queue))
        threads.append(thread)
        thread.start()

    return result


def run_stat_model_thread(min_count, limit_asm_to_top_n_isoforms,
                          min_asm_reads_by_sample, tolerance, max_iter, seed,
                          sample_to_group, covars, in_queue, out_queue,
                          signal_queue):
    initialize_rpy2_instance()
    r_packages = import_r_packages()

    while True:
        signal = rmats_long_utils.try_get_from_queue_without_wait(signal_queue)
        if signal is not None:
            return

        args = rmats_long_utils.try_get_from_queue_with_short_wait(in_queue)
        if args is None:
            continue

        asm_id = args['asm_id']
        gene_id = args['gene_id']
        read_grouping = args['read_grouping']
        results = run_stat_model_on_asm(asm_id, gene_id, read_grouping,
                                        min_count, limit_asm_to_top_n_isoforms,
                                        min_asm_reads_by_sample, tolerance,
                                        max_iter, seed, sample_to_group,
                                        covars, r_packages)
        out_queue.put(results)


def create_queue_asm_handler(in_queue, out_queue, threads, out_handles):
    def handler(asm_id, gene_id, read_grouping, status_counts, status_handler):
        # Attempt to queue this for the worker threads
        args = {
            'asm_id': asm_id,
            'gene_id': gene_id,
            'read_grouping': read_grouping
        }
        while True:
            was_put = rmats_long_utils.try_put_to_queue_without_wait(
                in_queue, args)
            if was_put:
                status_counts['pending'] += 1
                return

            # The in_queue could be full due to a full output queue.
            # Clear the output queue.
            while True:
                results = rmats_long_utils.try_get_from_queue_with_short_wait(
                    out_queue)
                if results is None:
                    break

                write_results(results, out_handles)
                status_counts['pending'] -= 1
                status_counts['completed'] += 1
                status_handler(status_counts)

            # The in_queue could be full because of errors in other threads.
            rmats_long_utils.raise_exception_if_thread_exited_early(threads)

    return handler


def create_process_asm_handler(min_count, limit_asm_to_top_n_isoforms,
                               min_asm_reads_by_sample, tolerance, max_iter,
                               seed, sample_to_group, covars, r_packages,
                               out_handles):
    def handler(asm_id, gene_id, read_grouping, status_counts, status_handler):
        results = run_stat_model_on_asm(asm_id, gene_id, read_grouping,
                                        min_count, limit_asm_to_top_n_isoforms,
                                        min_asm_reads_by_sample, tolerance,
                                        max_iter, seed, sample_to_group,
                                        covars, r_packages)
        write_results(results, out_handles)
        status_counts['completed'] += 1
        status_handler(status_counts)

    return handler


def create_status_handler(progress_every_n):
    def handler(status_counts, force=False):
        completed = status_counts['completed']
        old_next = status_counts.get('next', progress_every_n)
        should_print = force
        if completed == old_next:
            should_print = True
            status_counts['next'] = old_next + progress_every_n

        if should_print:
            print_with_timestamp('completed {} ASMs'.format(completed))

    return handler


def run_stat_model_with_out_handles(counts_tsv, abundance, covars,
                                    min_asm_reads_by_sample, group_1, group_2,
                                    group_1_name, group_2_name, num_threads,
                                    min_count, limit_asm_to_top_n_isoforms,
                                    tolerance, max_iter, seed,
                                    progress_every_n, out_handles):
    all_samples = group_1 + group_2
    sample_to_group = make_sample_to_group(group_1, group_2, group_1_name,
                                           group_2_name)
    thread_details = create_threads_and_queues(num_threads, min_count,
                                               limit_asm_to_top_n_isoforms,
                                               min_asm_reads_by_sample,
                                               tolerance, max_iter, seed,
                                               sample_to_group, covars)
    threads = thread_details['threads']
    in_queue = thread_details['in_queue']
    out_queue = thread_details['out_queue']
    signal_queue = thread_details['signal_queue']

    initialize_rpy2_instance()
    r_packages = import_r_packages()

    write_output_headers(out_handles)
    if threads:
        on_asm_handler = create_queue_asm_handler(in_queue, out_queue, threads,
                                                  out_handles)
    else:
        on_asm_handler = create_process_asm_handler(
            min_count, limit_asm_to_top_n_isoforms, min_asm_reads_by_sample,
            tolerance, max_iter, seed, sample_to_group, covars, r_packages,
            out_handles)

    try:
        run_stat_model_main_thread(counts_tsv, abundance, all_samples,
                                   sample_to_group, progress_every_n,
                                   out_queue, threads, out_handles,
                                   on_asm_handler)
    finally:
        cleanup_threads(threads, in_queue, out_queue, signal_queue)


def cleanup_threads(threads, in_queue, out_queue, signal_queue):
    if not threads:
        return

    # signal threads to stop
    for _ in threads:
        signal_queue.put(True)

    rmats_long_utils.drain_queue(in_queue)
    rmats_long_utils.drain_queue(out_queue)

    unjoined_threads = list()
    non_zero_exit_codes = list()
    for thread in threads:
        thread.join(1)  # wait at most 1 second for the thread to join
        exit_code = thread.exitcode
        if thread.exitcode is None:
            unjoined_threads.append(thread)
        elif exit_code != 0:
            non_zero_exit_codes.append(exit_code)

    messages = list()
    if unjoined_threads:
        messages.append('{} thread(s) were not joined'.format(
            len(unjoined_threads)))
        for thread in unjoined_threads:
            thread.terminate()

    if non_zero_exit_codes:
        num_non_zero = len(non_zero_exit_codes)
        messages.append('{} thread(s) had non-zero exit codes: {}'.format(
            num_non_zero, non_zero_exit_codes))

    if messages:
        raise Exception('. '.join(messages))


def process_counts_tsv(counts_tsv, all_samples, status_counts, on_asm_handler,
                       status_handler):
    # counts_tsv is sorted by asm_id, then read_id
    expected_headers = [
        'asm_id', 'gene_id', 'read_id', 'sample_id', 'isoform_id'
    ]
    gene_id = None
    current_asm = None
    read_rows = list()
    read_grouping = create_read_counts_initial_grouping()
    with open(counts_tsv, 'rt') as handle:
        for line_i, line in enumerate(handle):
            columns = rmats_long_utils.read_tsv_line(line)
            if line_i == 0:
                if columns != expected_headers:
                    raise Exception('Expected {} to have columns: {}'.format(
                        counts_tsv, expected_headers))

                continue

            row = dict(zip(expected_headers, columns))
            asm_id = row['asm_id']
            if (current_asm is None) or (asm_id != current_asm):
                if current_asm is not None:
                    update_read_counts_grouping(read_rows, read_grouping)
                    read_rows = list()
                    on_asm_handler(current_asm, gene_id, read_grouping,
                                   status_counts, status_handler)

                current_asm = asm_id
                read_grouping = create_read_counts_initial_grouping()

            row_gene_id = row['gene_id']
            if row_gene_id:
                gene_id = row_gene_id

            sample = row['sample_id']
            if sample not in all_samples:
                raise Exception('Sample {} not found in {}'.format(
                    sample, all_samples))

            row_read_id = row['read_id']
            if (not read_rows) or (row_read_id == read_rows[0]['read_id']):
                read_rows.append(row)
            else:
                update_read_counts_grouping(read_rows, read_grouping)
                read_rows = [row]

        if current_asm is not None:
            update_read_counts_grouping(read_rows, read_grouping)
            on_asm_handler(current_asm, gene_id, read_grouping, status_counts,
                           status_handler)


def read_abundance_file(abundance, all_samples):
    by_gene_by_sample_by_isoform = dict()
    expected_headers = ['transcript_ID', 'transcript_name', 'gene_ID']
    num_expected = len(expected_headers)
    headers = list()
    with open(abundance, 'rt') as handle:
        for line_i, line in enumerate(handle):
            columns = rmats_long_utils.read_tsv_line(line)
            if line_i == 0:
                found_headers = columns[:num_expected]
                if found_headers != expected_headers:
                    raise Exception(
                        'Expected {} to start with columns: {}'.format(
                            abundance, expected_headers))

                headers = columns
                found_samples = set(headers[num_expected:])
                all_samples_set = set(all_samples)
                extra_samples = found_samples - all_samples_set
                if extra_samples:
                    print(
                        'warning: skipping samples not assigned to a group: {}'
                        .format(sorted(extra_samples)),
                        file=sys.stderr)

                missing_samples = all_samples_set - found_samples
                if missing_samples:
                    raise Exception('Some samples not found in {}: {}'.format(
                        abundance, sorted(missing_samples)))

                continue

            row = dict(zip(headers, columns))
            gene = row['gene_ID']
            if gene == 'NA':
                continue

            isoform = row['transcript_ID']
            genes = gene.split(',')
            for gene in genes:
                by_sample_by_isoform = rmats_long_utils.try_get_or_set_default(
                    by_gene_by_sample_by_isoform, gene, dict())
                for sample in all_samples:
                    by_isoform = rmats_long_utils.try_get_or_set_default(
                        by_sample_by_isoform, sample, dict())
                    sample_abun = row[sample]
                    by_isoform[isoform] = float(sample_abun)

    return by_gene_by_sample_by_isoform


def process_abundance_tsv(abundance, all_samples, sample_to_group,
                          status_counts, on_asm_handler, status_handler):
    by_gene_by_sample_by_isoform = read_abundance_file(abundance, all_samples)
    for gene, by_sample_by_isoform in by_gene_by_sample_by_isoform.items():
        read_grouping = create_grouping_from_isoform_counts(
            by_sample_by_isoform, sample_to_group)
        on_asm_handler(gene, gene, read_grouping, status_counts,
                       status_handler)


def run_stat_model_main_thread(counts_tsv, abundance, all_samples,
                               sample_to_group, progress_every_n, out_queue,
                               threads, out_handles, on_asm_handler):
    status_counts = {'pending': 0, 'completed': 0}
    status_handler = create_status_handler(progress_every_n)
    if abundance is not None:
        process_abundance_tsv(abundance, all_samples, sample_to_group,
                              status_counts, on_asm_handler, status_handler)
    else:
        process_counts_tsv(counts_tsv, all_samples, status_counts,
                           on_asm_handler, status_handler)

    while status_counts['pending'] > 0:
        results = rmats_long_utils.try_get_from_queue_with_short_wait(
            out_queue)
        if results is None:
            rmats_long_utils.raise_exception_if_thread_exited_early(threads)
            continue

        write_results(results, out_handles)
        status_counts['pending'] -= 1
        status_counts['completed'] += 1
        status_handler(status_counts)

    status_handler(status_counts, force=True)


def adjust_p_values(p_values, r_packages):
    adjusted = r_packages['stats'].p_adjust(p_values, method='BH')
    return adjusted


def adjust_asm_p_values(lrtp_path, differential_asms_path, r_packages):
    with open(lrtp_path, 'rt') as lrtp_handle:
        gene_id_values = list()
        p_values = list()
        lr_values = list()
        df_values = list()
        asm_ids = list()
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(
                lrtp_handle):
            gene_id = row['gene_id']
            asm_id = row['asm_id']
            p_value = parse_float(row['lrt_p'])
            lr = parse_float(row['lr'])
            df = parse_float(row['df'])
            if p_value is None:
                p_value = math.nan

            gene_id_values.append(gene_id)
            p_values.append(p_value)
            lr_values.append(lr)
            df_values.append(df)
            asm_ids.append(asm_id)

    adjusted_pvalues = adjust_p_values(p_values, r_packages)
    with open(differential_asms_path, 'wt') as asm_handle:
        asm_headers = ['asm_id', 'gene_id', 'lr', 'df', 'pvalue', 'adj_pvalue']
        rmats_long_utils.write_tsv_line(asm_handle, asm_headers)
        for asm_id_i, asm_id in enumerate(asm_ids):
            gene_id = gene_id_values[asm_id_i]
            p_value = p_values[asm_id_i]
            adj_p_value = adjusted_pvalues[asm_id_i]
            lr = lr_values[asm_id_i]
            df = df_values[asm_id_i]
            rmats_long_utils.write_tsv_line(asm_handle, [
                asm_id, gene_id,
                rmats_long_utils.format_float(lr),
                rmats_long_utils.format_float(df),
                rmats_long_utils.format_float(p_value),
                rmats_long_utils.format_float(adj_p_value)
            ])


def adjust_isoform_p_values(coeff_path, differential_isoforms_path,
                            r_packages):
    with open(coeff_path, 'rt') as coeff_handle:
        gene_id_values = list()
        p_values = list()
        lr_values = list()
        df_values = list()
        isoform_ids = list()
        asm_ids = list()
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(
                coeff_handle):
            gene_id = row['gene_id']
            isoform_id = row['isoform_id']
            asm_id = row['asm_id']
            p_value = parse_float(row['p_value'])
            lr = parse_float(row['lr'])
            df = parse_float(row['df'])
            if p_value is None:
                p_value = math.nan

            gene_id_values.append(gene_id)
            p_values.append(p_value)
            lr_values.append(lr)
            df_values.append(df)
            isoform_ids.append(isoform_id)
            asm_ids.append(asm_id)

    adjusted_pvalues = adjust_p_values(p_values, r_packages)
    with open(differential_isoforms_path, 'wt') as isoform_handle:
        isoform_headers = [
            'asm_id', 'gene_id', 'isoform_id', 'lr', 'df', 'pvalue',
            'adj_pvalue'
        ]
        rmats_long_utils.write_tsv_line(isoform_handle, isoform_headers)
        for isoform_id_i, isoform_id in enumerate(isoform_ids):
            gene_id = gene_id_values[isoform_id_i]
            p_value = p_values[isoform_id_i]
            adj_p_value = adjusted_pvalues[isoform_id_i]
            lr = lr_values[isoform_id_i]
            df = df_values[isoform_id_i]
            asm_id = asm_ids[isoform_id_i]
            rmats_long_utils.write_tsv_line(isoform_handle, [
                asm_id, gene_id, isoform_id,
                rmats_long_utils.format_float(lr),
                rmats_long_utils.format_float(df),
                rmats_long_utils.format_float(p_value),
                rmats_long_utils.format_float(adj_p_value)
            ])


def calculate_adjusted_p_values(coeff_path, lrtp_path, differential_asms_path,
                                differential_isoforms_path):
    initialize_rpy2_instance()
    r_packages = import_r_packages()
    adjust_asm_p_values(lrtp_path, differential_asms_path, r_packages)
    adjust_isoform_p_values(coeff_path, differential_isoforms_path, r_packages)


def add_split_asm_id_colums(in_path, out_path):
    old_headers = None
    with open(in_path, 'rt') as in_handle:
        with open(out_path, 'wt') as out_handle:
            asm_id_i = None
            for line_i, line in enumerate(in_handle):
                columns = rmats_long_utils.read_tsv_line(line)
                if line_i == 0:
                    old_headers = columns
                    asm_id_i = old_headers.index('asm_id')
                    continue

                asm_id = columns[asm_id_i]
                chr_i, asm_i = asm_id.split('_')
                new_columns = [chr_i, asm_i] + columns
                rmats_long_utils.write_tsv_line(out_handle, new_columns)

    return old_headers


def sort_by_split_asm_id_columns(sort_buffer_size, in_path, out_path,
                                 temp_dir):
    chr_i_key_arg = '-k1,1n'
    asm_i_key_arg = '-k2,2n'
    env = {'LC_ALL': 'C'}  # to ensure sort order
    command = [
        'sort', '--buffer-size', sort_buffer_size, '--temporary-directory',
        temp_dir, chr_i_key_arg, asm_i_key_arg, '--output', out_path, in_path
    ]
    subprocess.run(command, env=env, check=True)


def remove_split_asm_id_columns(old_headers, in_path, out_path):
    with open(in_path, 'rt') as in_handle:
        with open(out_path, 'wt') as out_handle:
            rmats_long_utils.write_tsv_line(out_handle, old_headers)
            for line in in_handle:
                columns = rmats_long_utils.read_tsv_line(line)
                new_columns = columns[2:]
                rmats_long_utils.write_tsv_line(out_handle, new_columns)


# asm_id values are formatted like {chr_i}_{asm_i}.
# First add 2 new columns {chr_i} and {asm_i}.
# Then sort by those columns.
# Finally remove those columns.
def sort_by_asm_id(sort_buffer_size, in_file, temp_dir):
    split_path = os.path.join(temp_dir, 'split')
    sort_path = os.path.join(temp_dir, 'sort')
    old_headers = add_split_asm_id_colums(in_file, split_path)
    sort_by_split_asm_id_columns(sort_buffer_size, split_path, sort_path,
                                 temp_dir)
    remove_split_asm_id_columns(old_headers, sort_path, in_file)


def sort_output_files(sort_buffer_size, has_asm_ids, count_path, coeff_path,
                      lrtp_path, error_path, warn_path, differential_asms_path,
                      differential_isoforms_path, out_dir):
    if not has_asm_ids:
        return  # let detect_differential_isoforms.py sort

    with tempfile.TemporaryDirectory(suffix='_tmp',
                                     prefix='sort_',
                                     dir=out_dir) as temp_dir:
        sort_by_asm_id(sort_buffer_size, count_path, temp_dir)
        sort_by_asm_id(sort_buffer_size, coeff_path, temp_dir)
        sort_by_asm_id(sort_buffer_size, lrtp_path, temp_dir)
        sort_by_asm_id(sort_buffer_size, error_path, temp_dir)
        sort_by_asm_id(sort_buffer_size, warn_path, temp_dir)
        sort_by_asm_id(sort_buffer_size, differential_asms_path, temp_dir)
        sort_by_asm_id(sort_buffer_size, differential_isoforms_path, temp_dir)


def get_min_asm_reads_by_sample(sample_read_total_tsv, min_cpm_per_asm):
    min_asm_reads_by_sample = dict()
    if not (sample_read_total_tsv and min_cpm_per_asm):
        return min_asm_reads_by_sample

    total_by_sample = rmats_long_utils.parse_sample_totals_file(
        sample_read_total_tsv)
    for sample, total in total_by_sample.items():
        reads_per_cpm = total / 1e6
        min_reads_per_asm = reads_per_cpm * min_cpm_per_asm
        min_asm_reads_by_sample[sample] = min_reads_per_asm

    return min_asm_reads_by_sample


def run_stat_model(counts_tsv, abundance, covar_tsv, sample_read_total_tsv,
                   min_cpm_per_asm, group_1_path, group_2_path, group_1_name,
                   group_2_name, num_threads, min_count,
                   limit_asm_to_top_n_isoforms, tolerance, max_iter, seed,
                   progress_every_n, sort_buffer_size, out_dir):
    group_1 = parse_group_file(group_1_path)
    group_2 = parse_group_file(group_2_path)
    covars = parse_covar_tsv(covar_tsv)
    min_asm_reads_by_sample = get_min_asm_reads_by_sample(
        sample_read_total_tsv, min_cpm_per_asm)
    out_handles = dict()
    count_path = os.path.join(out_dir, 'count.tsv')
    coeff_path = os.path.join(out_dir, 'coeff.tsv')
    lrtp_path = os.path.join(out_dir, 'lrtp.tsv')
    error_path = os.path.join(out_dir, 'error.tsv')
    warn_path = os.path.join(out_dir, 'warning.tsv')
    try:
        out_handles['count'] = open(count_path, 'wt')
        out_handles['coeff'] = open(coeff_path, 'wt')
        out_handles['lrtp'] = open(lrtp_path, 'wt')
        out_handles['error'] = open(error_path, 'wt')
        out_handles['warn'] = open(warn_path, 'wt')
        run_stat_model_with_out_handles(
            counts_tsv, abundance, covars, min_asm_reads_by_sample, group_1,
            group_2, group_1_name, group_2_name, num_threads, min_count,
            limit_asm_to_top_n_isoforms, tolerance, max_iter, seed,
            progress_every_n, out_handles)
    finally:
        for handle in out_handles.values():
            if handle:
                handle.close()

    differential_asms_path = os.path.join(out_dir, 'differential_asms.tsv')
    differential_isoforms_path = os.path.join(out_dir,
                                              'differential_isoforms.tsv')
    calculate_adjusted_p_values(coeff_path, lrtp_path, differential_asms_path,
                                differential_isoforms_path)
    has_asm_ids = counts_tsv is not None
    sort_output_files(sort_buffer_size, has_asm_ids, count_path, coeff_path,
                      lrtp_path, error_path, warn_path, differential_asms_path,
                      differential_isoforms_path, out_dir)


def main():
    args = parse_args()
    rmats_long_utils.create_output_dir(args.out_dir)
    run_stat_model(args.counts_tsv, args.abundance, args.covar_tsv,
                   args.sample_read_total_tsv, args.min_cpm_per_asm,
                   args.group_1, args.group_2, args.group_1_name,
                   args.group_2_name, args.num_threads, args.min_isoform_reads,
                   args.limit_asm_to_top_n_isoforms, args.em_tolerance,
                   args.em_max_iter, args.random_seed, args.progress_every_n,
                   args.sort_buffer_size, args.out_dir)
    print_with_timestamp('run_stat_model.py finished')
    capture_r_output()  # ignore R exit messages


if __name__ == '__main__':
    main()
