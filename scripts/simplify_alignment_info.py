import argparse
import os
import os.path
import shutil
import subprocess
import tempfile

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Read a SAM or BAM file and output the'
                     ' splice junctions for each read'))
    parser.add_argument('--in-file',
                        required=True,
                        help='The path to a .sam or .bam file')
    parser.add_argument('--out-tsv',
                        required=True,
                        help='The path to write the output reads')
    parser.add_argument('--sort-buffer-size',
                        default='2G',
                        help=('Used for the --buffer-size argument of sort.'
                              ' Default: %(default)s'))

    return parser.parse_args()


def create_sjs_string(sjs):
    strings = list()
    for start, end in sjs:
        sj = '{}-{}'.format(start, end)
        strings.append(sj)

    return ','.join(strings)


def parse_sam_flag(flag_str):
    flag = int(flag_str)
    # bin returns a string that starts with 0b
    bit_str = bin(flag)[2:]
    # fill in 0's on the left
    padded = '{:0>12}'.format(bit_str)
    is_supplementary = padded[0] == '1'
    # is_duplicate = padded[1] == '1'
    # failed_filter = padded[2] == '1'
    is_secondary = padded[3] == '1'
    # is_last_segment = padded[4] == '1'
    # is_first_segment = padded[5] == '1'
    # next_segment_is_reversed = padded[6] == '1'
    is_reversed = padded[7] == '1'
    # next_segment_unmapped = padded[8] == '1'
    segment_unmapped = padded[9] == '1'
    # all_segments_aligned = padded[10] == '1'
    # has_multiple_segments = padded[11] == '1'

    is_primary = not (is_supplementary or is_secondary)
    strand = '-' if is_reversed else '+'
    return {
        'is_mapped': not segment_unmapped,
        'is_primary': is_primary,
        'strand': strand
    }


def get_end_and_sjs_from_cigar(start, cigar):
    pos = start
    count_str = list()
    sjs = list()
    for c in cigar:
        if c.isdigit():
            count_str.append(c)
            continue

        count = int(''.join(count_str))
        count_str = list()
        operation = c
        if operation in '=XMD':
            pos += count
        elif operation in 'N':
            sj_start = pos - 1
            pos += count
            sjs.append((sj_start, pos))
        elif operation in 'ISHP':
            pass
        else:
            raise Exception('Unrecognized cigar operation {} {}'.format(
                operation, cigar))

    # (pos - 1) because pos is where the next match would start
    return {'end': pos - 1, 'sjs': sjs}


def simplify_alignment(line):
    if line.startswith('@'):
        return None

    columns = rmats_long_utils.read_tsv_line(line)
    # read_id = columns[0]
    flag_str = columns[1]
    chr_name = columns[2]
    start_str = columns[3]
    # mapq = columns[4]
    cigar = columns[5]
    flag_details = parse_sam_flag(flag_str)
    if not (flag_details['is_mapped'] and flag_details['is_primary']):
        return None

    strand = flag_details['strand']
    start = int(start_str)
    coord_details = get_end_and_sjs_from_cigar(start, cigar)
    end = coord_details['end']
    sjs = coord_details['sjs']
    sjs_string = create_sjs_string(sjs)
    return [chr_name, str(start), str(end), sjs_string, strand]


def simplify_alignment_info(in_path, out_path):
    with open(out_path, 'wt') as out_handle:
        command = ['samtools', 'view', in_path]
        with subprocess.Popen(command, stdout=subprocess.PIPE,
                              text=True) as process:
            while True:
                line = process.stdout.readline()
                if not line:
                    break

                simplified = simplify_alignment(line)
                if simplified is None:
                    continue

                rmats_long_utils.write_tsv_line(out_handle, simplified)

            process.communicate()
            if process.returncode != 0:
                raise Exception('{} exited with return {}'.format(
                    command, process.returncode))


def sort_by_chr(in_path, buffer_size):
    dir_path = os.path.dirname(os.path.abspath(in_path))
    with tempfile.TemporaryDirectory(suffix='_tmp',
                                     prefix='sort_by_chr',
                                     dir=dir_path) as temp_dir:
        tmp_path = os.path.join(temp_dir, 'sorted.tmp')
        chr_key_arg = '-k1,1'
        start_key_arg = '-k2,2g'
        end_key_arg = '-k3,3g'
        env = {'LC_ALL': 'C'}  # to ensure sort order
        command = [
            'sort', '--buffer-size', buffer_size, chr_key_arg, start_key_arg,
            end_key_arg, '--output', tmp_path, in_path
        ]
        subprocess.run(command, env=env, check=True)
        shutil.move(tmp_path, in_path)


def main():
    args = parse_args()
    simplify_alignment_info(args.in_file, args.out_tsv)
    sort_by_chr(args.out_tsv, args.sort_buffer_size)
    print('simplify_alignment_info.py finished')


if __name__ == '__main__':
    main()
