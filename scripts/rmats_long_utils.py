import math
import os
import os.path
import subprocess
import sys


def get_python_executable():
    python_executable = sys.executable
    if not python_executable:
        python_executable = 'python'

    return python_executable


def get_script_dir():
    py_script_rel_path = sys.argv[0]
    py_script_abs_path = os.path.abspath(py_script_rel_path)
    script_dir = os.path.dirname(py_script_abs_path)
    return script_dir


def create_output_dir(dir_path):
    if os.path.exists(dir_path):
        if not os.path.isdir(dir_path):
            raise Exception('{} already exists and is a file'.format(dir_path))

        return

    os.makedirs(dir_path)


def parse_group_file(file_path):
    groups = None
    with open(file_path, 'rt') as handle:
        for line_i, line in enumerate(handle):
            if line and line_i != 0:
                print('expected only 1 line in {}'.format(file_path))

            groups = line.rstrip('\n').split(',')

    if not groups:
        raise Exception('could not parse groups from: {}'.format(file_path))

    return groups


def parse_gtf_line(line):
    if line.strip().startswith('#'):
        # skip comment lines
        return None

    columns = line.rstrip('\n').split('\t')
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

        existing_attributes = attributes.get(key)
        if existing_attributes:
            if isinstance(existing_attributes, list):
                existing_attributes.append(value)
            else:
                attributes[key] = [existing_attributes, value]
        else:
            attributes[key] = value

    return attributes


def get_gene_name_and_canonical_transcript_from_gtf(gene_id, gtf_path):
    found_gene_name = None
    canonical_transcript = None
    with open(gtf_path, 'rt') as gtf_handle:
        for line in gtf_handle:
            parsed = parse_gtf_line(line)
            if parsed is None:
                continue

            found_gene = parsed['attributes'].get('gene_id')
            if found_gene != gene_id:
                continue

            current_gene_name = parsed['attributes'].get('gene_name')
            if current_gene_name:
                found_gene_name = current_gene_name

            tags = parsed['attributes'].get('tag')
            if tags and 'Ensembl_canonical' in tags:
                transcript_id = parsed['attributes'].get('transcript_id')
                if transcript_id:
                    canonical_transcript = transcript_id

    return found_gene_name, canonical_transcript


def select_significant_transcripts(gene_id, diff_transcripts_path):
    most_significant_positive = None
    most_significant_negative = None
    with open(diff_transcripts_path, 'rt') as transcripts_handle:
        for row in row_iterator_for_tsv_with_header(transcripts_handle):
            if row['gene_id'] != gene_id:
                continue

            transcript = row['feature_id']
            pvalue = parse_float(row['adj_pvalue'])
            delta = parse_float(row['delta_isoform_proportion'])
            if math.isnan(pvalue) or math.isnan(delta):
                continue

            parsed = {'id': transcript, 'pvalue': pvalue, 'delta': delta}
            if delta > 0:
                if ((most_significant_positive is None)
                        or (pvalue < most_significant_positive['pvalue'])):
                    most_significant_positive = parsed
            if delta < 0:
                if ((most_significant_negative is None)
                        or (pvalue < most_significant_negative['pvalue'])):
                    most_significant_negative = parsed

    if most_significant_positive and most_significant_negative:
        pos_pvalue = most_significant_positive['pvalue']
        neg_pvalue = most_significant_negative['pvalue']
        if pos_pvalue < neg_pvalue:
            return [
                most_significant_positive['id'],
                most_significant_negative['id']
            ]

        return [
            most_significant_negative['id'], most_significant_positive['id']
        ]
    if most_significant_positive:
        return [most_significant_positive['id']]
    if most_significant_negative:
        return [most_significant_negative['id']]

    return list()


def row_iterator_for_tsv_with_header(in_handle):
    for line_i, line in enumerate(in_handle):
        columns = line.rstrip('\n').split('\t')
        if line_i == 0:
            headers = columns
            continue

        yield dict(zip(headers, columns))


def parse_abundance_file(abundance_path):
    counts_by_gene_by_transcript_by_sample = dict()
    total_by_gene_by_sample = dict()
    total_by_sample = dict()
    expected_headers = ['transcript_ID', 'transcript_name', 'gene_ID']
    with open(abundance_path, 'rt') as in_handle:
        for line_i, line in enumerate(in_handle):
            in_columns = line.rstrip('\n').split('\t')
            if line_i == 0:
                in_headers = in_columns
                first_headers = in_headers[:len(expected_headers)]
                if first_headers != expected_headers:
                    raise Exception('expected headers in {} to start with {}'
                                    ' but found {}'.format(
                                        abundance_path, expected_headers,
                                        first_headers))

                sample_names = in_headers[len(expected_headers):]
                continue

            row = dict(zip(in_headers, in_columns))
            gene = row['gene_ID']
            transcript = row['transcript_ID']
            if gene == 'NA':
                continue

            gene_total_by_sample = total_by_gene_by_sample.get(gene)
            if not gene_total_by_sample:
                gene_total_by_sample = dict()
                total_by_gene_by_sample[gene] = gene_total_by_sample

            counts_by_transcript_by_sample = counts_by_gene_by_transcript_by_sample.get(
                gene)
            if not counts_by_transcript_by_sample:
                counts_by_transcript_by_sample = dict()
                counts_by_gene_by_transcript_by_sample[gene] = (
                    counts_by_transcript_by_sample)

            counts_by_sample = counts_by_transcript_by_sample.get(transcript)
            if not counts_by_sample:
                counts_by_sample = dict()
                counts_by_transcript_by_sample[transcript] = counts_by_sample

            for sample in sample_names:
                count = parse_float(row[sample])
                old_total = total_by_sample.get(sample, 0)
                total_by_sample[sample] = old_total + count
                old_gene_total = gene_total_by_sample.get(sample, 0)
                gene_total_by_sample[sample] = old_gene_total + count
                counts_by_sample[sample] = count

    return {
        'counts_by_gene_by_transcript_by_sample':
        counts_by_gene_by_transcript_by_sample,
        'total_by_gene_by_sample': total_by_gene_by_sample,
        'total_by_sample': total_by_sample,
        'sample_names': sample_names
    }


def write_tsv_line(handle, columns):
    handle.write('{}\n'.format('\t'.join(columns)))


def parse_float(string):
    if string == 'NA':
        return math.nan

    return float(string)


def format_float(value, places=4):
    if math.isnan(value):
        return 'NA'

    return str(round(value, places))


def round_float_string(string):
    value = parse_float(string)
    return format_float(value)


def run_command(command):
    print('running: {}'.format(command))
    subprocess.run(command, check=True)
