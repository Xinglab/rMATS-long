import argparse
import os
import os.path

from rmats_long import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Check pairs of ASMs for splicing coupling')
    parser.add_argument('--pair-count-dir',
                        required=True,
                        help=('The directory of ASM pair isoform counts from'
                              ' count_reads_for_asm_pairs.py'))
    parser.add_argument('--diff-isoforms',
                        required=True,
                        help=('The differential_isoforms.tsv from'
                              ' detect_differential_isoforms.py'))
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
    parser.add_argument('--out-tsv',
                        required=True,
                        help='Where to write coupling results')

    args = parser.parse_args()
    return args


def get_most_abundant_isoform_by_asm(diff_isoforms_path, group_1_samples,
                                     group_2_samples):
    by_asm = dict()
    samples = group_1_samples + group_2_samples
    with open(diff_isoforms_path, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            asm = row['asm_id']
            isoform = row['isoform_id']
            total_prop = 0
            for sample in samples:
                key = '{}_proportion'.format(sample)
                prop_str = row[key]
                if prop_str == 'NA':
                    prop = 0
                else:
                    prop = float(prop_str)

                total_prop += prop

            prev = by_asm.get(asm)
            if (not prev) or (total_prop > prev['prop']):
                by_asm[asm] = {'isoform': isoform, 'prop': total_prop}

    return by_asm


def is_chr_id_file_name(name):
    prefix = 'chr_id_'
    suffix = '.tsv'
    return (name.startswith(prefix) and name.endswith(suffix))


def process_asm_pair_values(asm_1, asm_2, asm_pair_values, most_abun_by_asm,
                            out_handle):
    main_iso_1_details = most_abun_by_asm.get(asm_1)
    main_iso_2_details = most_abun_by_asm.get(asm_2)
    if None in [main_iso_1_details, main_iso_2_details]:
        return

    main_iso_1 = main_iso_1_details['isoform']
    main_iso_2 = main_iso_2_details['isoform']
    # separate counts for groups 1 and 2
    main_main = [0, 0]
    main_other = [0, 0]
    other_main = [0, 0]
    other_other = [0, 0]
    unclear = [0, 0]
    only_1 = [0, 0]
    only_2 = [0, 0]
    for iso_key, group_totals in asm_pair_values.items():
        isos_1, isos_2 = iso_key
        group_1_total, group_2_total = group_totals
        is_main_1 = main_iso_1 in isos_1
        is_main_2 = main_iso_2 in isos_2
        is_both_1 = is_main_1 and (len(isos_1) > 1)
        is_both_2 = is_main_2 and (len(isos_2) > 1)
        to_increment = None
        if isos_1 == ('', ):
            to_increment = only_2
        elif isos_2 == ('', ):
            to_increment = only_1
        elif is_both_1 or is_both_2:
            to_increment = unclear
        elif is_main_1 and is_main_2:
            to_increment = main_main
        elif is_main_1 and (not is_main_2):
            to_increment = main_other
        elif (not is_main_1) and is_main_2:
            to_increment = other_main
        else:
            to_increment = other_other

        to_increment[0] += group_1_total
        to_increment[1] += group_2_total

    for group_i in [0, 1]:
        m_m = main_main[group_i]
        m_o = main_other[group_i]
        o_m = other_main[group_i]
        o_o = other_other[group_i]
        unc = unclear[group_i]
        o_1 = only_1[group_i]
        o_2 = only_2[group_i]
        columns = [
            main_iso_1, main_iso_2, group_i, m_m, m_o, o_m, o_o, unc, o_1, o_2
        ]
        rmats_long_utils.write_tsv_line(out_handle, [str(x) for x in columns])


def get_group_totals_from_row(row, group_1_samples, group_2_samples):
    group_1_total = 0
    for sample in group_1_samples:
        group_1_total += int(row[sample])

    group_2_total = 0
    for sample in group_2_samples:
        group_2_total += int(row[sample])

    return group_1_total, group_2_total


def check_asm_coupling(pair_count_dir, group_1_samples, group_2_samples,
                       most_abun_by_asm, out_handle):
    names = os.listdir(pair_count_dir)
    for name in names:
        path = os.path.join(pair_count_dir, name)
        if not is_chr_id_file_name(name):
            continue

        prev_asm_1 = None
        prev_asm_2 = None
        asm_pair_values = dict()
        with open(path, 'rt') as in_handle:
            for row in rmats_long_utils.row_iterator_for_tsv_with_header(
                    in_handle):
                asm_pair = row['asm_pair']
                asm_1, asm_2 = asm_pair.split(',')
                if asm_pair_values and ((asm_1 != prev_asm_1) or
                                        (asm_2 != prev_asm_2)):
                    process_asm_pair_values(prev_asm_1, prev_asm_2,
                                            asm_pair_values, most_abun_by_asm,
                                            out_handle)
                    asm_pair_values = dict()

                prev_asm_1 = asm_1
                prev_asm_2 = asm_2
                asm_1_isos_str = row['asm_1_isoforms']
                asm_2_isos_str = row['asm_2_isoforms']
                asm_1_isos = asm_1_isos_str.split(',')
                asm_2_isos = asm_2_isos_str.split(',')
                group_1_total, group_2_total = get_group_totals_from_row(
                    row, group_1_samples, group_2_samples)
                key = tuple((tuple(asm_1_isos), tuple(asm_2_isos)))
                asm_pair_values[key] = (group_1_total, group_2_total)

        if asm_pair_values:
            process_asm_pair_values(prev_asm_1, prev_asm_2, asm_pair_values,
                                    most_abun_by_asm, out_handle)


def main():
    args = parse_args()
    group_1_samples = rmats_long_utils.parse_group_file(args.group_1)
    group_2_samples = rmats_long_utils.parse_group_file(args.group_2)
    most_abun_by_asm = get_most_abundant_isoform_by_asm(
        args.diff_isoforms, group_1_samples, group_2_samples)
    with open(args.out_tsv, 'wt') as out_handle:
        out_headers = [
            'asm_1_main', 'asm_2_main', 'group_i', 'main_main', 'main_other',
            'other_main', 'other_other', 'unclear', 'only_1', 'only_2'
        ]
        rmats_long_utils.write_tsv_line(out_handle, out_headers)
        check_asm_coupling(args.pair_count_dir, group_1_samples,
                           group_2_samples, most_abun_by_asm, out_handle)

    print('check_asm_coupling.py finished')


if __name__ == '__main__':
    main()
