import argparse
import os
import os.path
import subprocess
import sys
import tempfile


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Compare the structures of isoforms within a gene'))
    parser.add_argument(
        '--main-transcript-id',
        required=True,
        help='The transcript_id of the main isoform in the .gtf file')
    parser.add_argument('--gtf',
                        required=True,
                        help='The path to a .gtf describing the isoforms')
    parser.add_argument('--out-tsv',
                        required=True,
                        help='The path of the output file')

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
        classify_isoform_differences_with_temp_files(temp_files,
                                                     args.main_transcript_id,
                                                     args.gtf, args.out_tsv)
    finally:
        for temp_file in temp_files.values():
            if os.path.exists(temp_file):
                os.remove(temp_file)


def get_python_executable():
    python_executable = sys.executable
    if not python_executable:
        python_executable = 'python'

    return python_executable


def get_alt_script_path():
    py_script_rel_path = sys.argv[0]
    py_script_abs_path = os.path.abspath(py_script_rel_path)
    script_dir = os.path.dirname(py_script_abs_path)
    alt_script_path = os.path.join(script_dir, 'FindAltTSEvents.py')
    return alt_script_path


def parse_gtf_line(line):
    line = line.strip()
    if line.startswith('#'):
        # skip comment lines
        return None

    columns = line.split('\t')
    chr_name = columns[0]
    source = columns[1]
    feature = columns[2]
    start_str = columns[3]
    start = int(start_str)
    end_str = columns[4]
    end = int(end_str)
    score = columns[5]
    strand = columns[6]
    frame = columns[7]
    attributes_str = columns[8]
    attributes = parse_gtf_attributes(attributes_str)
    return {
        'chr': chr_name,
        'source': source,
        'feature': feature,
        'start': start,
        'end': end,
        'score': score,
        'strand': strand,
        'frame': frame,
        'attributes': attributes,
        'line': line
    }


def parse_gtf_attributes(attributes_str):
    attributes = dict()
    attribute_pairs = attributes_str.split(';')
    for attribute_pair in attribute_pairs:
        attribute_pair = attribute_pair.strip()
        first_space = attribute_pair.find(' ')
        if first_space <= 0:
            continue

        key = attribute_pair[:first_space]
        value = attribute_pair[first_space + 1:]
        if value[0] == '"' and value[-1] == '"':
            # remove quotes
            value = value[1:-1]

        attributes[key] = value

    return attributes


def filter_gtf_to_gene(gtf, main_transcript_id):
    gene_id = None
    with open(gtf, 'rt') as gtf_handle:
        for line in gtf_handle:
            parsed = parse_gtf_line(line)
            if not parsed:
                continue

            transcript_id = parsed['attributes'].get('transcript_id')
            if transcript_id == main_transcript_id:
                gene_id = parsed['attributes'].get('gene_id')
                break

    if gene_id is None:
        raise Exception(
            'Unable to find gene_id for {}'.format(main_transcript_id))

    lines_by_transcript = dict()
    with open(gtf, 'rt') as gtf_handle:
        for line in gtf_handle:
            parsed = parse_gtf_line(line)
            if not parsed:
                continue

            line_gene_id = parsed['attributes'].get('gene_id')
            transcript_id = parsed['attributes'].get('transcript_id')
            if transcript_id and line_gene_id == gene_id:
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


def write_tsv_line(handle, columns):
    handle.write('{}\n'.format('\t'.join(columns)))


def classify_isoform_differences_with_temp_files(temp_files,
                                                 main_transcript_id, gtf,
                                                 out_tsv):
    python = get_python_executable()
    script = get_alt_script_path()
    isoform_gtf = temp_files['isoform_gtf']
    isoform_tsv = temp_files['isoform_tsv']
    lines_by_transcript = filter_gtf_to_gene(gtf, main_transcript_id)

    main_transcript_lines = lines_by_transcript[main_transcript_id]
    headers = ['transcript1', 'transcript2', 'event', 'coordinates']
    with open(out_tsv, 'wt') as combined_handle:
        write_tsv_line(combined_handle, headers)
        for transcript, lines in lines_by_transcript.items():
            if transcript == main_transcript_id:
                continue

            write_gtf_lines(isoform_gtf, main_transcript_lines, lines)
            command = [python, script, '-i', isoform_gtf, '-o', isoform_tsv]
            print('running: {}'.format(command))
            subprocess.run(command, check=True)
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
