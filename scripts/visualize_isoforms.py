# Based on code from Yang Xu
import argparse
import os
import os.path
import subprocess
import sys


def transcript_colors():
    return [
        '#FF5F42',  # light red
        '#003E7F',  # dark blue
        '#0068AF',  # medium blue
        '#5495E1',  # light blue
        '#A1C2E8',  # lighter blue
        '#D9D9D9',  # grey
    ]


# The last color is used for "other"
MAX_TRANSCRIPTS = len(transcript_colors()) - 1


def group_colors():
    return [
        "#FF0018",  # bright red
        "#E68E19",  # orange
    ]


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Visualize the structure and abundance of isoforms'))
    parser.add_argument('--gene-id',
                        required=True,
                        help='The gene_id to visualize')
    parser.add_argument(
        '--abundance',
        required=True,
        help='The path to the abundance.esp file from ESPRESSO')
    parser.add_argument('--updated-gtf',
                        required=True,
                        help='The path to the updated.gtf file from ESPRESSO')
    parser.add_argument('--out-dir',
                        required=True,
                        help='The path to use as the output directory')
    parser.add_argument(
        '--plot-file-type',
        choices=['.pdf', '.png'],
        default='.pdf',
        help='The file type for output plots (default %(default)s))')
    parser.add_argument(
        '--main-transcript-id',
        required=True,
        help='The transcript_id of the main transcript to plot')
    parser.add_argument(
        '--max-transcripts',
        type=int,
        default=5,
        help='How many transcripts to plot individually.'
        ' The remaining transcripts in the gene will be grouped together'
        ' (max {}, default %(default)s)'.format(MAX_TRANSCRIPTS))
    # TODO intron scaling: find regions that only have introns -> scale them as 1/intron_scaling
    # is_exon = np.array(size)
    # for transcript -> for exon -> is_exon[coord] = True
    # coord_translation = np.array(size)
    # new_i = 0; for i, is_exon in is_exon: -> (only advance new_i at most 1 for every intron_scaling if not is_exon
    parser.add_argument(
        '--group-1',
        help=('The path to a file listing the sample names for group 1. The'
              ' file should have a single line with the sample names as a'
              ' comma separated list. The sample names should match with the'
              ' ESPRESSO abundance column names.'))
    parser.add_argument(
        '--group-2',
        help='The path to a file listing the sample names for group 2.')
    parser.add_argument('--group-1-name',
                        default='group 1',
                        help='A name for group 1 (default %(default)s)')
    parser.add_argument('--group-2-name',
                        default='group 2',
                        help='A name for group 2 (default %(default)s)')

    args = parser.parse_args()
    if args.max_transcripts > MAX_TRANSCRIPTS:
        parser.error('--max_transcripts was {} but can be at most {}'
                     ' (due to pre-defined colors)'.format(
                         args.max_transcripts, MAX_TRANSCRIPTS))

    if args.group_1 or args.group_2:
        if not args.group_1:
            parser.error('--group-2 given without --group-1')
        if not args.group_2:
            parser.error('--group-1 given without --group-2')

        args.group_1 = parse_group_file(args.group_1)
        args.group_2 = parse_group_file(args.group_2)

    return args


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


def read_abundance_file(abundance_path):
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
                count = float(row[sample])
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


def calculate_transcript_rows(counts_by_transcript_by_sample,
                              gene_total_by_sample, total_by_sample, group_1,
                              group_2, group_1_name, group_2_name):
    rows_by_transcript = dict()
    for transcript, counts_by_sample in counts_by_transcript_by_sample.items():
        counts_by_sample = counts_by_transcript_by_sample[transcript]
        for sample_i, sample in enumerate(group_1 + group_2):
            is_group_1 = sample_i < len(group_1)
            gene_total = gene_total_by_sample[sample]
            sample_total = total_by_sample[sample]
            count = counts_by_sample[sample]
            if sample_total == 0:
                cpm = 'NA'
                proportion = 'NA'
            else:
                cpm = (count * 1e6) / sample_total
                if gene_total == 0:
                    proportion = 'NA'
                else:
                    proportion = count / gene_total

            transcript_rows = rows_by_transcript.get(transcript)
            if not transcript_rows:
                transcript_rows = list()
                rows_by_transcript[transcript] = transcript_rows

            if is_group_1:
                group_name = group_1_name
            else:
                group_name = group_2_name

            transcript_rows.append({
                'sample': sample,
                'group': group_name,
                'cpm': cpm,
                'proportion': proportion,
            })

    return rows_by_transcript


# Sort transcripts by average proportion across samples.
# The sort order determines the color.
# The main_transcript_id is always first in the sort order.
def sort_transcripts(rows_by_transcript, main_transcript_id):
    average_proportion_by_transcript = dict()
    for transcript, rows in rows_by_transcript.items():
        if transcript == main_transcript_id:
            continue

        total_prop = 0
        num_samples_with_data = 0
        for row in rows:
            prop = row['proportion']
            if prop == 'NA':
                continue

            total_prop += prop
            num_samples_with_data += 1

        if num_samples_with_data == 0:
            average_proportion_by_transcript[transcript] = -1
        else:
            average = total_prop / num_samples_with_data
            average_proportion_by_transcript[transcript] = average

    def sort_key_func(transcript):
        avg_prop = average_proportion_by_transcript[transcript]
        # break ties with the transcript_id
        return (avg_prop, transcript)

    sorted_transcripts = sorted(average_proportion_by_transcript,
                                key=sort_key_func,
                                reverse=True)
    return [main_transcript_id] + sorted_transcripts


def calc_proportion_and_cpm(abundance_path, gene_id, main_transcript_id,
                            group_1, group_2, group_1_name, group_2_name,
                            cpm_and_proportion_path):
    abundance_details = read_abundance_file(abundance_path)
    counts_by_gene_by_transcript_by_sample = (
        abundance_details['counts_by_gene_by_transcript_by_sample'])
    total_by_gene_by_sample = abundance_details['total_by_gene_by_sample']
    total_by_sample = abundance_details['total_by_sample']
    sample_names = abundance_details['sample_names']
    if group_1 is None and group_2 is None:
        group_1 = sample_names
        group_2 = list()

    counts_by_transcript_by_sample = (
        counts_by_gene_by_transcript_by_sample[gene_id])
    gene_total_by_sample = total_by_gene_by_sample[gene_id]
    rows_by_transcript = calculate_transcript_rows(
        counts_by_transcript_by_sample, gene_total_by_sample, total_by_sample,
        group_1, group_2, group_1_name, group_2_name)
    sorted_transcripts = sort_transcripts(rows_by_transcript,
                                          main_transcript_id)
    with open(cpm_and_proportion_path, 'wt') as out_handle:
        out_headers = [
            'gene', 'transcript', 'sample', 'group', 'cpm', 'proportion'
        ]
        write_tsv_line(out_handle, out_headers)
        for transcript in sorted_transcripts:
            transcript_rows = rows_by_transcript[transcript]
            for row in transcript_rows:
                out_columns = [
                    gene_id, transcript, row['sample'], row['group'],
                    na_float_to_str(row['cpm']),
                    na_float_to_str(row['proportion'])
                ]
                write_tsv_line(out_handle, out_columns)


def na_float_to_str(na_float):
    if na_float == 'NA':
        return na_float

    return str(round(na_float, 4))


def write_tsv_line(handle, columns):
    handle.write('{}\n'.format('\t'.join(columns)))


def run_command(command):
    print('running: {}'.format(command))
    subprocess.run(command, check=True)


def plot_abundance(out_path, cpm_and_proportion_path, max_transcripts):
    py_script_rel_path = sys.argv[0]
    py_script_abs_path = os.path.abspath(py_script_rel_path)
    script_dir = os.path.dirname(py_script_abs_path)
    r_script_path = os.path.join(script_dir, 'visualize_isoforms.R')
    transcript_colors_string = ','.join(transcript_colors())
    group_colors_string = ','.join(group_colors())
    command = [
        'Rscript', r_script_path, cpm_and_proportion_path,
        str(max_transcripts), transcript_colors_string, group_colors_string,
        out_path
    ]
    run_command(command)


def plot_structure(out_path, proportion_path, gtf_path, gene_id,
                   main_transcript_id, max_transcripts):
    pass  # TODO


def visualize_isoforms(args):
    create_output_dir(args.out_dir)
    cpm_and_proportion_file = os.path.join(
        args.out_dir, '{}_cpm_and_proportion.tsv'.format(args.gene_id))
    abundance_plot_file = os.path.join(
        args.out_dir, '{}_abundance{}'.format(args.gene_id,
                                              args.plot_file_type))
    structure_plot_file = os.path.join(
        args.out_dir, '{}_structure{}'.format(args.gene_id,
                                              args.plot_file_type))

    calc_proportion_and_cpm(args.abundance, args.gene_id,
                            args.main_transcript_id, args.group_1,
                            args.group_2, args.group_1_name, args.group_2_name,
                            cpm_and_proportion_file)
    plot_abundance(abundance_plot_file, cpm_and_proportion_file,
                   args.max_transcripts)
    plot_structure(structure_plot_file, cpm_and_proportion_file,
                   args.updated_gtf, args.gene_id, args.main_transcript_id,
                   args.max_transcripts)


def create_output_dir(dir_path):
    if os.path.exists(dir_path):
        if not os.path.isdir(dir_path):
            raise Exception('{} already exists and is a file'.format(dir_path))

        return

    os.makedirs(dir_path)


def main():
    args = parse_args()
    visualize_isoforms(args)


if __name__ == '__main__':
    main()
