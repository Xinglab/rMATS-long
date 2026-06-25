import argparse
import os
import os.path

from rmats_long import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Check pairs of ASMs for splicing coupling')
    parser.add_argument(
        '--event-dir',
        required=True,
        help='The output directory from detect_splicing_events.py')
    parser.add_argument(
        '--asm-counts-dir',
        required=True,
        help='The output directory from count_reads_for_asms.py')
    parser.add_argument('--count-tsv',
                        required=True,
                        help='The count.tsv from rmats_long.py')
    parser.add_argument(
        '--group-1',
        help=('The path to a file listing the sample names for group 1. The'
              ' file should have a single line with the sample names as a'
              ' comma separated list.'))
    parser.add_argument(
        '--group-2',
        help='The path to a file listing the sample names for group 2')
    parser.add_argument('--group-1-name',
                        default='group 1',
                        help='A name for --group-1 (default %(default)s)')
    parser.add_argument('--group-2-name',
                        default='group 2',
                        help='A name for --group-2 (default %(default)s)')
    parser.add_argument(
        '--group-tsv',
        help=('A .tsv with the headers "sample" and "group".'
              ' Can be used instead of --group-1 and --group-2 to'
              ' define more than two groups'))
    parser.add_argument('--out-tsv',
                        required=True,
                        help='Where to write coupling results')
    parser.add_argument(
        '--isoform-out-tsv',
        required=True,
        help='Where to write the main isoform determined for each ASM')

    args = parser.parse_args()
    group_error = 'Use either --group-tsv or (--group-1 and --group-2)'
    if args.group_tsv:
        if args.group_1 or args.group_2:
            parser.error(group_error)
    elif not (args.group_1 and args.group_2):
        parser.error(group_error)

    return args


def get_sample_to_group_from_tsv(group_tsv, sample_to_group):
    with open(group_tsv, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            sample = row['sample']
            group = row['group']
            sample_to_group[sample] = group


def get_sample_to_group_from_group_files(group_1, group_2, group_1_name,
                                         group_2_name, sample_to_group):
    group_1_samples = rmats_long_utils.parse_group_file(group_1)
    group_2_samples = rmats_long_utils.parse_group_file(group_2)
    for sample in group_1_samples:
        sample_to_group[sample] = group_1_name

    for sample in group_2_samples:
        sample_to_group[sample] = group_2_name


def get_sample_grouping(group_1, group_2, group_1_name, group_2_name,
                        group_tsv):
    sample_to_group = dict()
    if group_tsv:
        get_sample_to_group_from_tsv(group_tsv, sample_to_group)
    else:
        get_sample_to_group_from_group_files(group_1, group_2, group_1_name,
                                             group_2_name, sample_to_group)

    num_samples_by_group = dict()
    for group in sample_to_group.values():
        old_count = num_samples_by_group.get(group, 0)
        num_samples_by_group[group] = old_count + 1

    return {
        'sample_to_group': sample_to_group,
        'num_samples_by_group': num_samples_by_group
    }


def write_out_headers(out_handle):
    headers = [
        'sample_id', 'group', 'asm_id_1', 'asm_id_2', 'high1_high2',
        'high1_low2', 'low1_high2', 'low1_low2'
    ]
    rmats_long_utils.write_tsv_line(out_handle, headers)


def write_isoform_out_headers(out_handle):
    headers = ['asm_id', 'isoform_id']
    rmats_long_utils.write_tsv_line(out_handle, headers)


def write_main_isoform_line(asm, main_isoform, out_handle):
    columns = [asm, main_isoform]
    rmats_long_utils.write_tsv_line(out_handle, columns)


def asm_is_after_target(target_asm, asm):
    asm_chr_id = int(asm['chr'])
    target_chr_id = int(target_asm['chr'])
    chr_is_after = asm_chr_id > target_chr_id
    chr_is_same = asm_chr_id == target_chr_id
    event_is_after = asm['event_i'] > target_asm['event_i']
    return chr_is_after or (chr_is_same and event_is_after)


def update_main_isoform_from_props(current_asm, num_samples_by_group,
                                   props_by_isoform_by_group,
                                   main_isoform_by_asm):
    top_average = 0
    for isoform, props_by_group in props_by_isoform_by_group.items():
        isoform_average = 0
        for group, props in props_by_group.items():
            num_samples = num_samples_by_group[group]
            group_average = sum(props) / num_samples
            isoform_average += group_average

        if isoform_average > top_average:
            top_average = isoform_average
            main_isoform_by_asm[current_asm] = isoform


def get_row_from_line_and_headers(line, headers):
    columns = rmats_long_utils.read_tsv_line(line)
    return dict(zip(headers, columns))


def get_main_isoform_for_asms(asms, sample_to_group, num_samples_by_group,
                              count_headers, count_handle_and_next_line):
    final_asm = asms[-1]
    parsed_final_asm = rmats_long_utils.parse_asm_id(final_asm)
    main_isoform_by_asm = dict()
    current_asm = None
    props_by_isoform_by_group = dict()
    count_line = count_handle_and_next_line.read_line()
    while count_line:
        count_row = get_row_from_line_and_headers(count_line, count_headers)
        asm_id = count_row['asm_id']
        parsed_asm = rmats_long_utils.parse_asm_id(asm_id)
        if asm_is_after_target(parsed_final_asm, parsed_asm):
            break

        count_line = count_handle_and_next_line.clear_and_read_line()
        if asm_id not in asms:
            continue

        sample_id = count_row['sample_id']
        group = sample_to_group[sample_id]
        isoform_id = count_row['isoform_id']
        prop_str = count_row['prop']
        prop = float(prop_str)
        if asm_id != current_asm:
            update_main_isoform_from_props(current_asm, num_samples_by_group,
                                           props_by_isoform_by_group,
                                           main_isoform_by_asm)
            current_asm = None
            props_by_isoform_by_group = dict()

        current_asm = asm_id
        props_by_group = rmats_long_utils.try_get_or_set_default(
            props_by_isoform_by_group, isoform_id, dict())
        props = rmats_long_utils.try_get_or_set_default(
            props_by_group, group, list())
        props.append(prop)

    if current_asm is not None:
        update_main_isoform_from_props(current_asm, num_samples_by_group,
                                       props_by_isoform_by_group,
                                       main_isoform_by_asm)
        current_asm = None
        props_by_isoform_by_group = dict()

    return main_isoform_by_asm


def seek_to_asm(gene_asms, asm, asm_to_offset, headers, handle):
    offset = asm_to_offset.get(asm)
    if offset is not None:
        handle.seek(offset)
        return True

    offset = handle.tell()
    line = handle.readline()
    current_asm = None
    while line:
        row = get_row_from_line_and_headers(line, headers)
        asm_id = row['asm_id']
        if asm_id != current_asm:
            current_asm = asm_id
            if asm_id in gene_asms:
                asm_to_offset[asm_id] = offset

            if asm_id == asm:
                handle.seek(offset)
                return True

        offset = handle.tell()
        line = handle.readline()

    return False


def get_row_for_key_value(key, value, read_count_headers,
                          handle_and_next_line):
    line = handle_and_next_line.read_line()
    if not line:
        return None

    row = get_row_from_line_and_headers(line, read_count_headers)
    found_value = row[key]
    if found_value != value:
        return None

    return row


def get_row_for_asm(asm, read_count_headers, handle_and_next_line):
    return get_row_for_key_value('asm_id', asm, read_count_headers,
                                 handle_and_next_line)


def get_row_for_read(read_id, read_count_headers, handle_and_next_line):
    read_id_str = str(read_id)
    return get_row_for_key_value('read_id', read_id_str, read_count_headers,
                                 handle_and_next_line)


def get_next_shared_read(asm_1, asm_2, read_count_headers,
                         handle_and_next_line_1, handle_and_next_line_2):
    row_1 = get_row_for_asm(asm_1, read_count_headers, handle_and_next_line_1)
    row_2 = get_row_for_asm(asm_2, read_count_headers, handle_and_next_line_2)
    while row_1 and row_2:
        read_1 = int(row_1['read_id'])
        read_2 = int(row_2['read_id'])
        if read_1 < read_2:
            handle_and_next_line_1.clear_line()
            row_1 = get_row_for_asm(asm_1, read_count_headers,
                                    handle_and_next_line_1)
            continue
        if read_2 < read_1:
            handle_and_next_line_2.clear_line()
            row_2 = get_row_for_asm(asm_2, read_count_headers,
                                    handle_and_next_line_2)
            continue

        sample_id = row_1['sample_id']
        return {'read_id': read_1, 'sample_id': sample_id}

    return None


def get_isoforms_for_read(read_id, read_count_headers, handle_and_next_line):
    isoforms = list()
    row = get_row_for_read(read_id, read_count_headers, handle_and_next_line)
    while row:
        isoform = row['isoform_id']
        isoforms.append(isoform)
        handle_and_next_line.clear_line()
        row = get_row_for_read(read_id, read_count_headers,
                               handle_and_next_line)

    return isoforms


def update_details_by_sample_from_isoforms(sample_id, main_isoform_1,
                                           main_isoform_2, isoforms_1,
                                           isoforms_2, details_by_sample):
    num_isos_1 = len(isoforms_1)
    num_isos_2 = len(isoforms_2)
    if main_isoform_1 in isoforms_1:
        p_high_1 = 1 / num_isos_1
    else:
        p_high_1 = 0

    if main_isoform_2 in isoforms_2:
        p_high_2 = 1 / num_isos_2
    else:
        p_high_2 = 0

    p_low_1 = 1 - p_high_1
    p_low_2 = 1 - p_high_2
    high_high = p_high_1 * p_high_2
    high_low = p_high_1 * p_low_2
    low_high = p_low_1 * p_high_2
    low_low = p_low_1 * p_low_2
    sample_details = details_by_sample.get(sample_id)
    if not sample_details:
        sample_details = {
            'high_high': 0,
            'high_low': 0,
            'low_high': 0,
            'low_low': 0
        }
        details_by_sample[sample_id] = sample_details

    sample_details['high_high'] += high_high
    sample_details['high_low'] += high_low
    sample_details['low_high'] += low_high
    sample_details['low_low'] += low_low


def check_coupling_for_asm_pair(asm_1, asm_2, main_isoform_1, main_isoform_2,
                                sample_to_group, read_count_headers,
                                read_count_handle_1, read_count_handle_2,
                                out_handle):
    details_by_sample = dict()
    handle_and_next_line_1 = rmats_long_utils.HandleAndNextLine(
        read_count_handle_1)
    handle_and_next_line_2 = rmats_long_utils.HandleAndNextLine(
        read_count_handle_2)
    while True:
        read_details = get_next_shared_read(asm_1, asm_2, read_count_headers,
                                            handle_and_next_line_1,
                                            handle_and_next_line_2)
        if not read_details:
            break

        read_id = read_details['read_id']
        sample_id = read_details['sample_id']
        isoforms_1 = get_isoforms_for_read(read_id, read_count_headers,
                                           handle_and_next_line_1)
        isoforms_2 = get_isoforms_for_read(read_id, read_count_headers,
                                           handle_and_next_line_2)
        update_details_by_sample_from_isoforms(sample_id, main_isoform_1,
                                               main_isoform_2, isoforms_1,
                                               isoforms_2, details_by_sample)

    for sample, details in details_by_sample.items():
        group = sample_to_group[sample]
        high_high = rmats_long_utils.format_float(details['high_high'])
        high_low = rmats_long_utils.format_float(details['high_low'])
        low_high = rmats_long_utils.format_float(details['low_high'])
        low_low = rmats_long_utils.format_float(details['low_low'])
        columns = [
            sample, group, asm_1, asm_2, high_high, high_low, low_high, low_low
        ]
        rmats_long_utils.write_tsv_line(out_handle, columns)

    offset_1 = handle_and_next_line_1.offset
    offset_2 = handle_and_next_line_2.offset
    next_offset = max(offset_1, offset_2)
    return next_offset


def check_asm_coupling_for_gene_asms(gene_asms, sample_to_group,
                                     num_samples_by_group, count_headers,
                                     count_handle_and_next_line, asm_to_offset,
                                     read_count_headers, read_count_handle_1,
                                     read_count_handle_2, out_handle,
                                     isoform_out_handle):
    if not gene_asms:
        return

    main_isoform_by_asm = get_main_isoform_for_asms(
        gene_asms, sample_to_group, num_samples_by_group, count_headers,
        count_handle_and_next_line)
    asms_with_main = list()
    for asm in gene_asms:
        main_isoform = main_isoform_by_asm.get(asm)
        if main_isoform:
            asms_with_main.append(asm)
            write_main_isoform_line(asm, main_isoform, isoform_out_handle)
        else:
            print('No main isoform for {}'.format(asm))

    next_offset = None
    for asm_1_i, asm_1 in enumerate(asms_with_main):
        main_isoform_1 = main_isoform_by_asm[asm_1]
        for asm_2 in asms_with_main[asm_1_i + 1:]:
            main_isoform_2 = main_isoform_by_asm[asm_2]
            found_asm_1 = seek_to_asm(gene_asms, asm_1, asm_to_offset,
                                      read_count_headers, read_count_handle_1)
            if not found_asm_1:
                print('No compatible read counts found for {}'.format(asm_1))
                continue

            found_asm_2 = seek_to_asm(gene_asms, asm_2, asm_to_offset,
                                      read_count_headers, read_count_handle_2)
            if not found_asm_2:
                print('No compatible read counts found for {}'.format(asm_2))
                continue

            next_offset = check_coupling_for_asm_pair(
                asm_1, asm_2, main_isoform_1, main_isoform_2, sample_to_group,
                read_count_headers, read_count_handle_1, read_count_handle_2,
                out_handle)

    asm_to_offset.clear()
    if next_offset:
        read_count_handle_1.seek(next_offset)
        read_count_handle_2.seek(next_offset)


def check_asm_coupling_with_handles(sample_to_group, num_samples_by_group,
                                    counts_chr_path, event_handle,
                                    count_headers, count_handle_and_next_line,
                                    out_handle, isoform_out_handle):
    asm_to_offset = dict()
    # Open two handles, one for each asm.
    # Separate handles may have better performance than constantly seeking.
    with open(counts_chr_path, 'rt') as read_count_handle_1:
        read_count_header_line = read_count_handle_1.readline()
        read_count_headers = rmats_long_utils.read_tsv_line(
            read_count_header_line)
        with open(counts_chr_path, 'rt') as read_count_handle_2:
            current_gene_i = None
            gene_asms = list()
            for event_row in rmats_long_utils.row_iterator_for_tsv_with_header(
                    event_handle):
                asm_id = event_row['asm_id']
                gene_i = event_row['gene_i']
                is_strict = event_row['is_strict'] == 'True'
                if gene_i != current_gene_i:
                    check_asm_coupling_for_gene_asms(
                        gene_asms, sample_to_group, num_samples_by_group,
                        count_headers, count_handle_and_next_line,
                        asm_to_offset, read_count_headers, read_count_handle_1,
                        read_count_handle_2, out_handle, isoform_out_handle)
                    current_gene_i = gene_i
                    gene_asms = list()

                if not is_strict:
                    continue

                gene_asms.append(asm_id)

            if current_gene_i is not None:
                check_asm_coupling_for_gene_asms(
                    gene_asms, sample_to_group, num_samples_by_group,
                    count_headers, count_handle_and_next_line, asm_to_offset,
                    read_count_headers, read_count_handle_1,
                    read_count_handle_2, out_handle, isoform_out_handle)


def check_asm_coupling(event_dir, asm_counts_dir, count_tsv, group_1, group_2,
                       group_1_name, group_2_name, group_tsv, out_tsv,
                       isoform_out_tsv):
    grouping_details = get_sample_grouping(group_1, group_2, group_1_name,
                                           group_2_name, group_tsv)
    sample_to_group = grouping_details['sample_to_group']
    num_samples_by_group = grouping_details['num_samples_by_group']
    chr_id_mapping = rmats_long_utils.parse_chr_mapping(event_dir)
    id_to_chr = chr_id_mapping['id_to_chr']
    chr_ids = sorted(id_to_chr)
    with open(out_tsv, 'wt') as out_handle:
        write_out_headers(out_handle)
        with open(isoform_out_tsv, 'wt') as isoform_out_handle:
            write_isoform_out_headers(isoform_out_handle)
            with open(count_tsv, 'rt') as count_handle:
                count_handle_and_next_line = rmats_long_utils.HandleAndNextLine(
                    count_handle)
                count_header_line = count_handle_and_next_line.read_line()
                count_headers = rmats_long_utils.read_tsv_line(
                    count_header_line)
                count_handle_and_next_line.clear_line()
                for chr_id in chr_ids:
                    event_chr_path = rmats_long_utils.get_chr_file_path(
                        event_dir, chr_id)
                    counts_chr_path = rmats_long_utils.get_chr_file_path(
                        asm_counts_dir, chr_id)
                    has_events = os.path.exists(event_chr_path)
                    has_counts = os.path.exists(counts_chr_path)
                    if not (has_events and has_counts):
                        continue

                    with open(event_chr_path, 'rt') as event_handle:
                        check_asm_coupling_with_handles(
                            sample_to_group, num_samples_by_group,
                            counts_chr_path, event_handle, count_headers,
                            count_handle_and_next_line, out_handle,
                            isoform_out_handle)


def main():
    args = parse_args()
    check_asm_coupling(args.event_dir, args.asm_counts_dir, args.count_tsv,
                       args.group_1, args.group_2, args.group_1_name,
                       args.group_2_name, args.group_tsv, args.out_tsv,
                       args.isoform_out_tsv)

    print('check_asm_coupling.py finished')


if __name__ == '__main__':
    main()
