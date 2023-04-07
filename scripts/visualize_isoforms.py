# Based on code from Yang Xu
import argparse
import os
import os.path
import sys

import matplotlib
# Select a non-interctive backend immediately after import
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

import rmats_long_utils


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
        '--gene-name',
        help=('The name for the gene (used as plot title). If not given then'
              ' the gene_name from --gencode-gtf will be used. If no other'
              ' name is found then --gene-id is used as a default'))
    parser.add_argument(
        '--abundance',
        required=True,
        help='The path to the abundance.esp file from ESPRESSO')
    parser.add_argument('--updated-gtf',
                        required=True,
                        help='The path to the updated.gtf file from ESPRESSO')
    parser.add_argument(
        '--gencode-gtf',
        help=('The path to a gencode annotation.gtf file. Can be used to'
              ' identify the gene_name and Ensemble canonical isoform'))
    parser.add_argument(
        '--diff-transcripts',
        help=('The path to the differential transcript results. Can be used to'
              ' determine --main-transcript-ids'))
    parser.add_argument('--out-dir',
                        required=True,
                        help='The path to use as the output directory')
    parser.add_argument(
        '--plot-file-type',
        choices=['.pdf', '.png'],
        default='.pdf',
        help='The file type for output plots (default %(default)s))')
    parser.add_argument(
        '--main-transcript-ids',
        help=(
            'A comma separated list of transcript IDs to plot as the main'
            ' transcripts. If not given then the most significant isoform'
            ' from --diff-transcripts, a second significant isoform with a'
            ' delta proportion in the opposite direction, and the Ensembl'
            ' canonical isoform from --gencode-gtf will be used if possible'))
    parser.add_argument(
        '--max-transcripts',
        type=int,
        default=5,
        help='How many transcripts to plot individually.'
        ' The remaining transcripts in the gene will be grouped together'
        ' (max {}, default %(default)s)'.format(MAX_TRANSCRIPTS))
    parser.add_argument(
        '--intron-scaling',
        type=int,
        default=1,
        help=('The factor to use to reduce intron length in the plot.'
              ' A value of 2 would reduce introns to 1/2 of the'
              ' original plot length (default %(default)s)'))
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

        args.group_1 = rmats_long_utils.parse_group_file(args.group_1)
        args.group_2 = rmats_long_utils.parse_group_file(args.group_2)

    if args.main_transcript_ids:
        args.main_transcript_ids = args.main_transcript_ids.split(',')

    return args


def calculate_transcript_rows(counts_by_transcript_by_sample,
                              gene_total_by_sample, total_by_sample, group_1,
                              group_2, group_1_name, group_2_name):
    rows_by_transcript = dict()
    for transcript, counts_by_sample in counts_by_transcript_by_sample.items():
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
# The main_transcript_ids are first in the sort order.
def sort_transcripts(rows_by_transcript, main_transcript_ids):
    average_proportion_by_transcript = dict()
    for transcript, rows in rows_by_transcript.items():
        if transcript in main_transcript_ids:
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
    return main_transcript_ids + sorted_transcripts


def calc_proportion_and_cpm(abundance_path, gene_id, main_transcript_ids,
                            group_1, group_2, group_1_name, group_2_name,
                            cpm_and_proportion_path):
    abundance_details = rmats_long_utils.parse_abundance_file(abundance_path)
    counts_by_gene_by_transcript_by_sample = (
        abundance_details['counts_by_gene_by_transcript_by_sample'])
    total_by_gene_by_sample = abundance_details['total_by_gene_by_sample']
    total_by_sample = abundance_details['total_by_sample']
    sample_names = abundance_details['sample_names']
    if group_1 is None and group_2 is None:
        group_1 = sample_names
        group_2 = list()

    if sample_names[0] in group_1:
        group_name_of_first_sample = group_1_name
    else:
        group_name_of_first_sample = group_2_name

    counts_by_transcript_by_sample = (
        counts_by_gene_by_transcript_by_sample[gene_id])
    gene_total_by_sample = total_by_gene_by_sample[gene_id]
    rows_by_transcript = calculate_transcript_rows(
        counts_by_transcript_by_sample, gene_total_by_sample, total_by_sample,
        group_1, group_2, group_1_name, group_2_name)
    sorted_transcripts = sort_transcripts(rows_by_transcript,
                                          main_transcript_ids)
    with open(cpm_and_proportion_path, 'wt') as out_handle:
        out_headers = [
            'gene', 'transcript', 'sample', 'group', 'cpm', 'proportion'
        ]
        rmats_long_utils.write_tsv_line(out_handle, out_headers)
        for transcript in sorted_transcripts:
            transcript_rows = rows_by_transcript.get(transcript)
            if not transcript_rows:
                # The main_transcript_ids might not have any read counts.
                # Add a line to the output with 0 values.
                out_columns = [
                    gene_id, transcript, sample_names[0],
                    group_name_of_first_sample, '0', '0'
                ]
                rmats_long_utils.write_tsv_line(out_handle, out_columns)
                continue

            for row in transcript_rows:
                out_columns = [
                    gene_id, transcript, row['sample'], row['group'],
                    rmats_long_utils.round_float_string(row['cpm']),
                    rmats_long_utils.round_float_string(row['proportion'])
                ]
                rmats_long_utils.write_tsv_line(out_handle, out_columns)


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
    rmats_long_utils.run_command(command)


# The proportion file has the transcripts in sorted order.
def read_transcript_ids(proportion_path):
    transcript_ids = list()
    with open(proportion_path, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            transcript = row['transcript']
            if transcript in transcript_ids:
                continue

            transcript_ids.append(transcript)

    return transcript_ids


def read_transcript_details(gtf_path, transcript_ids):
    details_by_transcript = dict()
    with open(gtf_path, 'rt') as handle:
        for line in handle:
            parsed = rmats_long_utils.parse_gtf_line(line)
            if parsed is None:
                continue

            transcript_id = parsed['attributes'].get('transcript_id')
            if transcript_id not in transcript_ids:
                continue

            details = details_by_transcript.get(transcript_id)
            if not details:
                details = dict()
                details_by_transcript[transcript_id] = details

            strand = parsed['strand']
            if strand != '.':
                details['strand'] = strand

            if parsed['feature'] != 'exon':
                continue

            exons = details.get('exons')
            if not exons:
                exons = list()
                details['exons'] = exons

            exons.append((parsed['start'], parsed['end']))

    return details_by_transcript


def add_transcript_details(transcript_details, gtf_path, transcript_ids):
    new_details = read_transcript_details(gtf_path, transcript_ids)
    for transcript in transcript_ids:
        if transcript in transcript_details:
            continue

        found_details = new_details.get(transcript)
        if found_details:
            transcript_details[transcript] = found_details


def get_exon_length(exon):
    return (exon[1] - exon[0]) + 1


def get_intron_length(exon_1, exon_2):
    return (exon_2[0] - exon_1[1]) - 1


def get_min_and_max_coords(transcript_details):
    min_coord = None
    max_coord = None
    for details in transcript_details.values():
        exons = details['exons']
        exons.sort()
        transcript_min_coord = exons[0][0]
        transcript_max_coord = exons[-1][1]
        if min_coord is None:
            min_coord = transcript_min_coord
        else:
            min_coord = min(min_coord, transcript_min_coord)

        if max_coord is None:
            max_coord = transcript_max_coord
        else:
            max_coord = max(max_coord, transcript_max_coord)

    return min_coord, max_coord


def get_plot_coords_by_transcript(transcript_details, min_coord,
                                  region_length):
    plot_coords_by_transcript = dict()
    for transcript, details in transcript_details.items():
        plot_coords = list()
        is_minus_strand = details['strand'] == '-'
        exons = details['exons']
        transcript_min_coord = exons[0][0]
        start_pos = transcript_min_coord - min_coord
        current_pos = start_pos
        for exon_i in range(len(exons) - 1):
            exon_length = get_exon_length(exons[exon_i])
            intron_length = get_intron_length(exons[exon_i], exons[exon_i + 1])
            exon_start = current_pos
            intron_start = current_pos + exon_length
            exon_end = intron_start - 1
            intron_end = exon_end + intron_length
            plot_coords.append((exon_start, exon_end))
            plot_coords.append((intron_start, intron_end))
            current_pos = intron_end + 1

        # last_exon
        exon_start = current_pos
        exon_length = get_exon_length(exons[-1])
        exon_end = exon_start + (exon_length - 1)
        plot_coords.append((exon_start, exon_end))
        # Flip the coordinates by subtracting from the region_length.
        # This puts the 5' end on the left.
        if is_minus_strand:
            adjusted_plot_coords = list()
            for start, end in reversed(plot_coords):
                adjusted_start = (region_length - 1) - start
                adjusted_end = (region_length - 1) - end
                adjusted_plot_coords.append((adjusted_end, adjusted_start))

            plot_coords = adjusted_plot_coords

        plot_coords_by_transcript[transcript] = plot_coords

    return plot_coords_by_transcript


def plot_transcripts(ax, transcripts_to_plot, plot_coords_by_transcript,
                     region_length, line_space, exon_height, colors,
                     exon_edge_color, exon_line_width, intron_line_color,
                     intron_line_width, font_size):
    exon_z_order = 100
    five_prime_x_val = -0.03 * region_length
    colors_to_plot = colors[:len(transcripts_to_plot)]
    reversed_transcripts = list(reversed(transcripts_to_plot))
    reversed_colors = list(reversed(colors_to_plot))
    for transcript_i, transcript_id in enumerate(reversed_transcripts):
        color = reversed_colors[transcript_i]
        plot_coords = plot_coords_by_transcript[transcript_id]
        mid_y_val = line_space * (transcript_i + 1)
        top_y_val = mid_y_val + (exon_height / 2)
        bottom_y_val = mid_y_val - (exon_height / 2)
        text_y_val = mid_y_val - (0.1 * line_space) + exon_height
        five_prime_y_val = mid_y_val - 2
        is_exon = True
        for region in plot_coords:
            start, end = region
            if is_exon:
                is_exon = False

                exon_x = np.array([start, end])
                exon_y_top = np.array([top_y_val] * 2)
                exon_y_bottom = np.array([bottom_y_val] * 2)
                ax.fill_between(exon_x,
                                exon_y_top,
                                exon_y_bottom,
                                facecolor=color,
                                edgecolor=exon_edge_color,
                                linewidth=exon_line_width,
                                zorder=exon_z_order)
            else:
                is_exon = True
                intron_x = np.array([start, end])
                intron_y = np.array([mid_y_val] * 2)
                plt.plot(intron_x,
                         intron_y,
                         color=intron_line_color,
                         linewidth=intron_line_width)

        plt.text(0, text_y_val, transcript_id, fontsize=font_size)
        plt.text(five_prime_x_val,
                 five_prime_y_val,
                 "5'",
                 fontsize=font_size + 1)


def get_exon_regions(plot_coords_by_transcript, region_length):
    coord_is_exon = np.full(region_length, False)
    for plot_coords in plot_coords_by_transcript.values():
        is_exon = True
        for region in plot_coords:
            if not is_exon:
                is_exon = True
                continue

            start, end = region
            coord_is_exon[start:end + 1] = True
            is_exon = False

    return coord_is_exon


def get_intron_coord_translation(coord_is_exon, region_length, intron_scaling):
    num_removed_units = 0
    coord_translation = np.full(region_length, 0)
    new_coord = 0
    in_intron = False
    current_scaling_skip_count = 0
    max_to_skip = intron_scaling - 1
    for i, is_exon in enumerate(coord_is_exon):
        coord_translation[i] = new_coord
        # Leave exon regions alone
        if is_exon:
            in_intron = False
            new_coord += 1
            continue

        # Keep the first coord of an intron region
        if not in_intron:
            in_intron = True
            current_scaling_skip_count = 0
            new_coord += 1
            continue

        # Skip max_to_skip before using more space for an intron region
        if current_scaling_skip_count == max_to_skip:
            current_scaling_skip_count = 0
            new_coord += 1
            continue

        # skip
        current_scaling_skip_count += 1
        num_removed_units += 1

    return coord_translation, num_removed_units


def translate_coords(plot_coords_by_transcript, coord_translation):
    for plot_coords in plot_coords_by_transcript.values():
        for i, region in enumerate(plot_coords):
            start, end = region
            new_start = coord_translation[start]
            new_end = coord_translation[end]
            plot_coords[i] = (new_start, new_end)


# Regions that only have introns will be scaled to 1/intron_scaling
def apply_intron_scaling(plot_coords_by_transcript, intron_scaling,
                         region_length):
    if intron_scaling == 1:
        return region_length

    coord_is_exon = get_exon_regions(plot_coords_by_transcript, region_length)
    coord_translation, num_removed_units = get_intron_coord_translation(
        coord_is_exon, region_length, intron_scaling)
    translate_coords(plot_coords_by_transcript, coord_translation)

    return region_length - num_removed_units


def plot_structure(out_path, proportion_path, gtf_path, gencode_gtf_path,
                   gene_name, main_transcript_ids, max_transcripts,
                   intron_scaling):
    colors = transcript_colors()
    proportion_transcript_ids = read_transcript_ids(proportion_path)
    all_transcript_ids = set(proportion_transcript_ids).union(
        set(main_transcript_ids))
    transcript_details = read_transcript_details(gtf_path, all_transcript_ids)
    if gencode_gtf_path:
        add_transcript_details(transcript_details, gencode_gtf_path,
                               all_transcript_ids)

    found_transcripts = set(transcript_details)
    if found_transcripts != all_transcript_ids:
        raise Exception('Did not find all expected transcripts: {}.'
                        ' Missing: {}'.format(
                            all_transcript_ids,
                            all_transcript_ids.difference(found_transcripts)))

    transcripts_to_plot = main_transcript_ids
    for transcript_id in proportion_transcript_ids:
        if transcript_id not in transcripts_to_plot:
            transcripts_to_plot.append(transcript_id)

    transcripts_to_plot = transcripts_to_plot[:max_transcripts]
    num_transcripts = len(transcripts_to_plot)

    fig = plt.figure(figsize=(12, 4), dpi=300)
    ax = fig.add_subplot(111)
    exon_height = 6
    line_space = 13
    intron_line_width = 0.5
    intron_line_color = 'black'
    exon_line_width = 0.5
    exon_edge_color = 'black'
    font_size = 13
    min_coord, max_coord = get_min_and_max_coords(transcript_details)
    region_length = (max_coord - min_coord) + 1
    plot_coords_by_transcript = get_plot_coords_by_transcript(
        transcript_details, min_coord, region_length)
    region_length = apply_intron_scaling(plot_coords_by_transcript,
                                         intron_scaling, region_length)
    plot_transcripts(ax, transcripts_to_plot, plot_coords_by_transcript,
                     region_length, line_space, exon_height, colors,
                     exon_edge_color, exon_line_width, intron_line_color,
                     intron_line_width, font_size)

    low_x_lim = -0.1 * region_length
    high_x_lim = 1.1 * region_length
    low_y_lim = 0
    high_y_lim = line_space * (num_transcripts + 1.5)
    mid_x_lim = (high_x_lim + low_x_lim) / 2
    title_x_val = mid_x_lim
    title_y_val = line_space * (num_transcripts + 1)
    plt.text(title_x_val,
             title_y_val,
             gene_name,
             fontsize=font_size + 1,
             horizontalalignment='center')
    ax.set_xlim(low_x_lim, high_x_lim)
    ax.set_ylim(low_y_lim, high_y_lim)
    plt.axis('off')
    # Let the plot take up the whole figure.
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0, hspace=0)
    plt.savefig(out_path, dpi=300, pad_inches=0)


def determine_main_transcripts_and_gene_name(gene_id, orig_gene_name,
                                             orig_main_transcript_ids,
                                             gencode_gtf_path,
                                             diff_transcripts_path):
    gene_name = None
    canonical_transcript = None
    if gencode_gtf_path:
        gene_name, canonical_transcript = (
            rmats_long_utils.get_gene_name_and_canonical_transcript_from_gtf(
                gene_id, gencode_gtf_path))

    if orig_gene_name is not None:
        gene_name = orig_gene_name
    elif gene_name is None:
        gene_name = gene_id

    if orig_main_transcript_ids:
        return orig_main_transcript_ids, gene_name

    main_transcript_ids = list()
    if diff_transcripts_path:
        main_transcript_ids.extend(
            rmats_long_utils.select_significant_transcripts(
                gene_id, diff_transcripts_path))

    if ((canonical_transcript
         and (canonical_transcript not in main_transcript_ids))):
        main_transcript_ids.append(canonical_transcript)

    return main_transcript_ids, gene_name


def visualize_isoforms(args):
    rmats_long_utils.create_output_dir(args.out_dir)
    main_transcript_ids, gene_name = determine_main_transcripts_and_gene_name(
        args.gene_id, args.gene_name, args.main_transcript_ids,
        args.gencode_gtf, args.diff_transcripts)

    cpm_and_proportion_file = os.path.join(
        args.out_dir, '{}_cpm_and_proportion.tsv'.format(args.gene_id))
    abundance_plot_file = os.path.join(
        args.out_dir, '{}_abundance{}'.format(args.gene_id,
                                              args.plot_file_type))
    structure_plot_file = os.path.join(
        args.out_dir, '{}_structure{}'.format(args.gene_id,
                                              args.plot_file_type))

    calc_proportion_and_cpm(args.abundance, args.gene_id, main_transcript_ids,
                            args.group_1, args.group_2, args.group_1_name,
                            args.group_2_name, cpm_and_proportion_file)
    plot_abundance(abundance_plot_file, cpm_and_proportion_file,
                   args.max_transcripts)
    plot_structure(structure_plot_file, cpm_and_proportion_file,
                   args.updated_gtf, args.gencode_gtf, gene_name,
                   main_transcript_ids, args.max_transcripts,
                   args.intron_scaling)


def main():
    args = parse_args()
    visualize_isoforms(args)


if __name__ == '__main__':
    main()
