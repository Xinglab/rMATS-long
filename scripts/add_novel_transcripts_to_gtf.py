import argparse
import shutil

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add ESPRESSO novel transcripts to a gtf")
    parser.add_argument('--espresso-gtf',
                        required=True,
                        help='the updated.gtf from ESPRESSO')
    parser.add_argument('--other-gtf',
                        required=True,
                        help='a gtf to add novel transcripts to')
    parser.add_argument('--out-gtf',
                        required=True,
                        help='where to write the output gtf')

    return parser.parse_args()


def add_novel_transcripts_to_gtf(in_handle, out_handle):
    for line in in_handle:
        if line.startswith('#'):
            continue

        columns = rmats_long_utils.read_tsv_line(line)
        source = columns[1]
        if source == 'novel_isoform':
            out_handle.write(line)


def main():
    args = parse_args()
    shutil.copy(args.other_gtf, args.out_gtf)
    with open(args.espresso_gtf, 'rt') as in_handle:
        with open(args.out_gtf, 'at') as out_handle:
            add_novel_transcripts_to_gtf(in_handle, out_handle)


if __name__ == '__main__':
    main()
