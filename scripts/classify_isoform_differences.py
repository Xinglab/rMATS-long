import argparse
import os
import os.path
import tempfile

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Compare the structures of isoforms within a gene'))
    parser.add_argument(
        '--main-transcript-id',
        required=True,
        help='The transcript_id of the main isoform in the .gtf file')
    parser.add_argument('--updated-gtf',
                        required=True,
                        help='The path to the updated.gtf file from ESPRESSO')
    parser.add_argument(
        '--gencode-gtf',
        help=('The path to a gencode annotation.gtf file. Can be used to'
              ' compare against isoforms not detected by ESPRESSO'))
    parser.add_argument('--out-tsv',
                        required=True,
                        help='The path of the output file')
    parser.add_argument(
        '--second-transcript-id',
        help='If given, only compare the main transcript to this transcript')

    return parser.parse_args()


def get_temp_files(out_tsv):
    temp_files = dict()
    out_tsv_abs_path = os.path.abspath(out_tsv)
    out_dir = os.path.dirname(out_tsv_abs_path)
    out_basename = os.path.basename(out_tsv_abs_path)
    out_prefix, out_extension = os.path.splitext(out_basename)
    name_suffix_pairs = [('isoform_gtf', '_tmp_isoform.gtf'),
                         ('isoform_tsv', '_tmp_isoform.tsv')]
    for name, suffix in name_suffix_pairs:
        with tempfile.NamedTemporaryFile(suffix=suffix,
                                         prefix=out_prefix,
                                         dir=out_dir,
                                         delete=False) as handle:
            temp_files[name] = handle.name

    return temp_files


def classify_isoform_differences(args):
    temp_files = dict()
    try:
        temp_files = get_temp_files(args.out_tsv)
        classify_isoform_differences_with_temp_files(
            temp_files, args.main_transcript_id, args.second_transcript_id,
            args.updated_gtf, args.gencode_gtf, args.out_tsv)
    finally:
        for temp_file in temp_files.values():
            if os.path.exists(temp_file):
                os.remove(temp_file)


def get_alt_script_path():
    script_dir = rmats_long_utils.get_script_dir()
    alt_script_path = os.path.join(script_dir, 'FindAltTSEvents.py')
    return alt_script_path


def filter_gtf_to_gene(updated_gtf, gencode_gtf, main_transcript_id):
    gene_id = find_gene_id_from_gtf(main_transcript_id, updated_gtf)
    if gene_id is None:
        gene_id = find_gene_id_from_gtf(main_transcript_id, gencode_gtf)

    if gene_id is None:
        raise Exception(
            'Unable to find gene_id for {}'.format(main_transcript_id))

    gene_ids = gene_id.split(',')
    lines_by_transcript = get_gtf_lines_by_transcript(updated_gtf, gene_ids)
    if gencode_gtf:
        gencode_lines_by_transcript = get_gtf_lines_by_transcript(
            gencode_gtf, gene_ids)
        for transcript, lines in gencode_lines_by_transcript.items():
            if transcript not in lines_by_transcript:
                lines_by_transcript[transcript] = lines

    return lines_by_transcript


def find_gene_id_from_gtf(transcript_id, gtf_path):
    with open(gtf_path, 'rt') as gtf_handle:
        for line in gtf_handle:
            parsed = rmats_long_utils.parse_gtf_line(line)
            if not parsed:
                continue

            found_transcript_id = parsed['attributes'].get('transcript_id')
            if found_transcript_id == transcript_id:
                gene_id = parsed['attributes'].get('gene_id')
                return gene_id

    return None


def get_gtf_lines_by_transcript(gtf_path, gene_ids):
    lines_by_transcript = dict()
    with open(gtf_path, 'rt') as gtf_handle:
        for line in gtf_handle:
            parsed = rmats_long_utils.parse_gtf_line(line)
            if not parsed:
                continue

            line_gene_id_str = parsed['attributes'].get('gene_id')
            line_gene_ids = line_gene_id_str.split(',')
            transcript_id = parsed['attributes'].get('transcript_id')
            matched_genes = set(line_gene_ids).intersection(set(gene_ids))
            if transcript_id and matched_genes:
                transcript_lines = lines_by_transcript.get(transcript_id)
                if not transcript_lines:
                    transcript_lines = list()
                    lines_by_transcript[transcript_id] = transcript_lines

                transcript_lines.append(parsed)

    return lines_by_transcript


def write_gtf_lines(file_name, main_lines, other_lines):
    with open(file_name, 'wt') as handle:
        for parsed in main_lines + other_lines:
            handle.write('{}\n'.format(parsed['line']))


def classify_isoform_differences_with_temp_files(temp_files,
                                                 main_transcript_id,
                                                 second_transcript_id,
                                                 updated_gtf, gencode_gtf,
                                                 out_tsv):
    python = rmats_long_utils.get_python_executable()
    script = get_alt_script_path()
    isoform_gtf = temp_files['isoform_gtf']
    isoform_tsv = temp_files['isoform_tsv']
    lines_by_transcript = filter_gtf_to_gene(updated_gtf, gencode_gtf,
                                             main_transcript_id)
    main_transcript_lines = lines_by_transcript[main_transcript_id]
    if second_transcript_id:
        second_transcript_lines = lines_by_transcript[second_transcript_id]
        write_gtf_lines(isoform_gtf, main_transcript_lines,
                        second_transcript_lines)
        command = [python, script, '-i', isoform_gtf, '-o', out_tsv]
        rmats_long_utils.run_command(command)
        return

    command = [python, script, '-i', isoform_gtf, '-o', isoform_tsv]
    headers = ['transcript1', 'transcript2', 'event', 'coordinates']
    with open(out_tsv, 'wt') as combined_handle:
        rmats_long_utils.write_tsv_line(combined_handle, headers)
        for transcript, lines in lines_by_transcript.items():
            if transcript == main_transcript_id:
                continue

            write_gtf_lines(isoform_gtf, main_transcript_lines, lines)
            rmats_long_utils.run_command(command)
            with open(isoform_tsv, 'rt') as single_handle:
                for line_i, line in enumerate(single_handle):
                    if line_i == 0:
                        continue  # skip header

                    combined_handle.write(line)


def main():
    args = parse_args()
    classify_isoform_differences(args)


if __name__ == '__main__':
    main()
