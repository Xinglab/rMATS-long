import argparse
import os
import os.path

from rmats_long import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Count the number of transcripts compatible with each ASM')
    parser.add_argument(
        '--event-dir',
        required=True,
        help='The output directory from detect_splicing_events.py for ASMs')
    parser.add_argument(
        '--transcript-compat-dir',
        required=True,
        help=('The output directory from count_reads_for_asms.py using'
              ' ASMs as events and transcripts as alignments'))
    parser.add_argument('--out-tsv',
                        required=True,
                        help='Where to write output counts')

    return parser.parse_args()


def load_asm_path_counts(event_dir):
    asm_to_path_count = dict()
    chr_id_mapping = rmats_long_utils.parse_chr_mapping(event_dir)
    id_to_chr = chr_id_mapping['id_to_chr']
    chr_ids = sorted(id_to_chr)
    for chr_id in chr_ids:
        tsv_path = rmats_long_utils.get_chr_file_path(event_dir, chr_id)
        if not os.path.exists(tsv_path):
            continue

        with open(tsv_path, 'rt') as handle:
            for row in rmats_long_utils.row_iterator_for_tsv_with_header(
                    handle):
                asm_id = row['asm_id']
                isoform_ids_str = row['isoform_ids']
                isoform_ids = rmats_long_utils.parse_isoform_ids_str(
                    isoform_ids_str)
                num_paths = len(isoform_ids)
                asm_to_path_count[asm_id] = num_paths

    return asm_to_path_count


def load_asm_compatibility(asm_to_path_count, asm_compat_dir):
    asm_to_compat = dict()
    chr_id_mapping = rmats_long_utils.parse_chr_mapping(asm_compat_dir)
    id_to_chr = chr_id_mapping['id_to_chr']
    chr_ids = sorted(id_to_chr)
    for chr_id in chr_ids:
        tsv_path = rmats_long_utils.get_chr_file_path(asm_compat_dir, chr_id)
        if not os.path.exists(tsv_path):
            continue

        with open(tsv_path, 'rt') as handle:
            for row in rmats_long_utils.row_iterator_for_tsv_with_header(
                    handle):
                transcript = row['read_id']
                asm_iso = row['isoform_id']
                asm_iso_details = rmats_long_utils.parse_asm_isoform_id(
                    asm_iso)
                asm_id = asm_iso_details['asm']
                if asm_id not in asm_to_path_count:
                    print('unknown ASM: {}'.format(asm_id))
                    continue

                compat_transcripts = asm_to_compat.get(asm_id)
                if not compat_transcripts:
                    compat_transcripts = set()
                    asm_to_compat[asm_id] = compat_transcripts

                compat_transcripts.add(transcript)

    return asm_to_compat


def asm_to_key(asm):
    parsed = rmats_long_utils.parse_asm_id(asm)
    chr_id_str = parsed['chr']
    chr_id = int(chr_id_str)
    event_i = parsed['event_i']
    return [chr_id, event_i]


def output_results(asm_to_path_count, asm_to_compat, out_path):
    with open(out_path, 'wt') as handle:
        headers = [
            'asm', 'num_asm_paths', 'num_transcripts_compatible_with_asm'
        ]
        rmats_long_utils.write_tsv_line(handle, headers)
        asms = sorted(asm_to_path_count, key=asm_to_key)
        for asm in asms:
            path_count = asm_to_path_count[asm]
            compat_transcripts = asm_to_compat.get(asm, set())
            columns = [asm, path_count, len(compat_transcripts)]
            columns = [str(x) for x in columns]
            rmats_long_utils.write_tsv_line(handle, columns)


def main():
    args = parse_args()
    asm_to_path_count = load_asm_path_counts(args.event_dir)
    asm_to_compat = load_asm_compatibility(asm_to_path_count,
                                           args.transcript_compat_dir)
    output_results(asm_to_path_count, asm_to_compat, args.out_tsv)
    print('count_transcripts_for_asms.py finished')


if __name__ == '__main__':
    main()
