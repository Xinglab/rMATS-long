# Based on code from Yang Xu
import argparse
import os
import os.path
import stat
import sys

import matplotlib
# Select a non-interactive backend immediately after import
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

import rmats_long_utils

FORCE_775 = False  # workaround for a filesystem issue
COLORS = {
    'beige': '#E6DAA6',
    'blue': '#0072B2',
    'bluishgreen': '#009E73',
    'bright red': '#FF0018',
    'orange': '#E69F00',
    'white': '#FFFFFF',
}


def get_transcript_colors(sorted_transcripts, main_transcript_ids,
                          max_transcripts):
    num_main = len(main_transcript_ids)
    num_total = len(sorted_transcripts)
    num_non_main = num_total - num_main
    colors = main_transcript_colors(num_main)
    num_selected_non_main = max(0, max_transcripts - num_main)
    colors.extend(other_transcript_colors(num_selected_non_main))
    if num_non_main > num_selected_non_main:
        colors.append(collapsed_transcripts_color())

    return colors


def main_transcript_colors(num_colors):
    colors = [
        COLORS['bright red'],
        COLORS['blue'],
    ]
    if num_colors > len(colors):
        raise Exception('Only {} colors defined, but {} requested'.format(
            len(colors), num_colors))
    return colors[:num_colors]


def other_transcript_colors(num_colors):
    # equally spaced grey values
    min_value = 0
    max_value = 255
    increment = max_value / (num_colors + 1)
    colors = list()
    value = min_value
    for _ in range(num_colors):
        value += increment
        rounded_value = round(value)
        # Always 2 characters. Zero padded if needed.
        hex_code = '{:02X}'.format(rounded_value)
        colors.append('#{x}{x}{x}'.format(x=hex_code))

    return colors


def group_colors():
    return [
        COLORS['bluishgreen'],
        COLORS['orange'],
    ]


def collapsed_transcripts_color():
    return COLORS['beige']


def parse_args():
    parser = argparse.ArgumentParser(
        description='Visualize the structure and abundance of isoforms')
    parser.add_argument('--gene-id',
                        required=True,
                        help='The gene_id to visualize')
    parser.add_argument(
        '--gene-name',
        help=('The name for the gene (used as plot title). If not given then'
              ' the gene_name from --gencode-gtf will be used. If no other'
              ' name is found then --gene-id is used as a default'))
    parser.add_argument('--asm-id', help='The asm_id to use with --asm-gtf')
    parser.add_argument(
        '--abundance', help='The path to the abundance.esp file from ESPRESSO')
    parser.add_argument('--updated-gtf',
                        required=True,
                        help='The path to the updated.gtf file from ESPRESSO')
    parser.add_argument(
        '--gencode-gtf',
        help=('The path to a gencode annotation.gtf file. Can be used to'
              ' identify the gene_name and Ensembl canonical isoform'))
    parser.add_argument(
        '--asm-gtf',
        help=('The path to a .gtf file with ASM transcripts.'
              ' Can be used to select which gene transcripts to plot'))
    parser.add_argument(
        '--diff-transcripts',
        help=('The path to the differential transcript results. Can be used to'
              ' determine --main-transcript-ids'))
    parser.add_argument('--out-dir',
                        required=True,
                        help='The path to use as the output directory')
    parser.add_argument(
        '--plot-file-type',
        choices=['.pdf', '.png', 'all'],
        default='.png',
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
        default=7,
        help=('How many transcripts to plot individually.'
              ' The remaining transcripts in the gene will be grouped together'
              ' (default %(default)s)'))
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
    parser.add_argument('--is-asm', action='store_true', help='Use ASM data')
    parser.add_argument(
        '--start-coord',
        type=int,
        help='Indicate the start of a region at this coordinate')
    parser.add_argument('--end-coord',
                        type=int,
                        help='Indicate the end of a region at this coordinate')
    parser.add_argument(
        '--graph-file',
        help=('The path to graph_{chr}.txt which has the splice'
              ' graph details for this gene'))

    args = parser.parse_args()
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
# The main_transcript_ids are first in the sort order,
# then the canonical_transcript if it's not a main transcript.
def sort_transcripts_by_proportion(rows_by_transcript, main_transcript_ids,
                                   canonical_transcript):
    average_proportion_by_transcript = dict()
    for transcript, rows in rows_by_transcript.items():
        if (((transcript in main_transcript_ids)
             or (transcript == canonical_transcript))):
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
    final_sort = list(main_transcript_ids)
    if (((canonical_transcript
          and (canonical_transcript not in main_transcript_ids)))):
        final_sort.append(canonical_transcript)

    final_sort.extend(sorted_transcripts)
    return final_sort


def get_score_features_from_exons(exons):
    splice_sites = set()
    splice_junctions = set()
    prev_end = None
    transcript_start = None
    transcript_end = None
    for exon in exons:
        start, end = exon
        if transcript_start is None:
            transcript_start = start

        transcript_end = end
        if prev_end is not None:
            splice_sites.add(prev_end)
            splice_sites.add(start)
            splice_junctions.add((prev_end, start))

        prev_end = end

    exons = set(exons)
    return {
        'start': transcript_start,
        'end': transcript_end,
        'splice_sites': splice_sites,
        'splice_junctions': splice_junctions,
        'exons': exons
    }


def get_score_features_from_parsed_gtf(parsed_details, sig_transcripts):
    splice_sites = set()
    splice_junctions = set()
    exons = set()

    sig_splice_sites = set()
    sig_splice_junctions = set()
    sig_exons = set()
    if parsed_details:
        by_transcript = parsed_details.get('by_transcript', dict())
        for transcript, details in by_transcript.items():
            transcript_features = get_score_features_from_exons(
                details['exons'])
            splice_sites.update(transcript_features['splice_sites'])
            splice_junctions.update(transcript_features['splice_junctions'])
            exons.update(transcript_features['exons'])
            if transcript in sig_transcripts:
                sig_splice_sites.update(transcript_features['splice_sites'])
                sig_splice_junctions.update(
                    transcript_features['splice_junctions'])
                sig_exons.update(transcript_features['exons'])

    return {
        'splice_sites': splice_sites,
        'splice_junctions': splice_junctions,
        'exons': exons,
        'sig_splice_sites': sig_splice_sites,
        'sig_splice_junctions': sig_splice_junctions,
        'sig_exons': sig_exons
    }


def initialize_used_values(earliest_start, latest_end):
    return {
        'earliest_start': earliest_start,
        'latest_end': latest_end,
        'used_earliest_start': False,
        'used_latest_end': False,
        'splice_sites': set(),
        'splice_junctions': set(),
        'exons': set()
    }


def update_used_values(transcript_features, used_values):
    if transcript_features['start'] == used_values['earliest_start']:
        used_values['used_earliest_start'] = True

    if transcript_features['end'] == used_values['latest_end']:
        used_values['used_latest_end'] = True

    used_values['splice_sites'].update(transcript_features['splice_sites'])
    used_values['splice_junctions'].update(
        transcript_features['splice_junctions'])
    used_values['exons'].update(transcript_features['exons'])


def get_earliest_start_and_latest_end(pending_transcripts):
    earliest_start = None
    latest_end = None
    for details in pending_transcripts.values():
        transcript_features = get_score_features_from_exons(details['exons'])
        start = transcript_features['start']
        end = transcript_features['end']
        if (earliest_start is None) or (start < earliest_start):
            earliest_start = start

        if (latest_end is None) or (end > latest_end):
            latest_end = end

    return earliest_start, latest_end


def get_score_from_features(transcript_features, used_values, asm_features):
    start = transcript_features['start']
    end = transcript_features['end']
    splice_sites = transcript_features['splice_sites']
    splice_junctions = transcript_features['splice_junctions']
    exons = transcript_features['exons']

    earliest_start = used_values['earliest_start']
    latest_end = used_values['latest_end']
    used_earliest_start = used_values['used_earliest_start']
    used_latest_end = used_values['used_latest_end']
    used_splice_sites = used_values['splice_sites']
    used_splice_junctions = used_values['splice_junctions']
    used_exons = used_values['exons']

    asm_splice_sites = asm_features['splice_sites']
    asm_splice_junctions = asm_features['splice_junctions']
    asm_exons = asm_features['exons']
    sig_asm_splice_sites = asm_features['sig_splice_sites']
    sig_asm_splice_junctions = asm_features['sig_splice_junctions']
    sig_asm_exons = asm_features['sig_exons']

    sj_score = 1
    exon_score = 2
    ss_score = 4
    endpoint_score = 100
    asm_mult = 10
    sig_asm_mult = 100
    score = 0
    if (not used_earliest_start) and (start == earliest_start):
        score += endpoint_score

    if (not used_latest_end) and (end == latest_end):
        score += endpoint_score

    for splice_site in splice_sites:
        if splice_site in used_splice_sites:
            continue

        if splice_site in sig_asm_splice_sites:
            score += ss_score * sig_asm_mult
        elif splice_site in asm_splice_sites:
            score += ss_score * asm_mult
        else:
            score += ss_score

    for splice_junction in splice_junctions:
        if splice_junction in used_splice_junctions:
            continue

        if splice_junction in sig_asm_splice_junctions:
            score += sj_score * sig_asm_mult
        elif splice_junction in asm_splice_junctions:
            score += sj_score * asm_mult
        else:
            score += sj_score

    for exon in exons:
        if exon in used_exons:
            continue

        if exon in sig_asm_exons:
            score += exon_score * sig_asm_mult
        elif exon in asm_exons:
            score += exon_score * asm_mult
        else:
            score += exon_score

    return score


# Put the main transcripts and the canonical transcript first.
# Then continue adding the transcript with the best "score".
# Transcripts will score based on features not yet found in the sort order.
def sort_transcripts_using_structure(main_transcript_ids, canonical_transcript,
                                     sig_asm_transcripts,
                                     parsed_gencode_details,
                                     parsed_updated_details,
                                     parsed_asm_details):
    asm_features = get_score_features_from_parsed_gtf(parsed_asm_details,
                                                      sig_asm_transcripts)
    pending_transcripts = combine_transcript_details(parsed_gencode_details,
                                                     parsed_updated_details)
    earliest_start, latest_end = get_earliest_start_and_latest_end(
        pending_transcripts)
    used_values = initialize_used_values(earliest_start, latest_end)

    sort_order = list()
    for initial_transcript in main_transcript_ids + [canonical_transcript]:
        details = pending_transcripts.get(initial_transcript)
        if not details:
            continue

        initial_features = get_score_features_from_exons(details['exons'])
        sort_order.append(initial_transcript)
        del pending_transcripts[initial_transcript]
        update_used_values(initial_features, used_values)

    while pending_transcripts:
        top_transcript = None
        top_transcript_features = None
        top_score = 0
        for transcript, details in pending_transcripts.items():
            transcript_features = get_score_features_from_exons(
                details['exons'])
            score = get_score_from_features(transcript_features, used_values,
                                            asm_features)
            if ((top_transcript is None or score > top_score
                 or (score == top_score and transcript < top_transcript))):
                top_transcript = transcript
                top_transcript_features = transcript_features
                top_score = score

        sort_order.append(top_transcript)
        del pending_transcripts[top_transcript]
        update_used_values(top_transcript_features, used_values)

    return sort_order


def calc_proportion_and_cpm(abundance_path, gene_id, main_transcript_ids,
                            canonical_transcript, group_1, group_2,
                            group_1_name, group_2_name, is_asm,
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
    sorted_transcripts = sort_transcripts_by_proportion(
        rows_by_transcript, main_transcript_ids, canonical_transcript)
    with open(cpm_and_proportion_path, 'wt') as out_handle:
        first_header = 'gene'
        if is_asm:
            first_header = 'asm_id'

        out_headers = [
            first_header, 'transcript', 'sample', 'group', 'cpm', 'proportion'
        ]
        rmats_long_utils.write_tsv_line(out_handle, out_headers)
        for transcript in sorted_transcripts:
            transcript_rows = rows_by_transcript.get(transcript)
            if not transcript_rows:
                # The main transcripts and canonical transcript might not
                # have any read counts.
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

    if FORCE_775:
        mode_775 = stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH
        os.chmod(cpm_and_proportion_path, mode_775)

    return sorted_transcripts


def plot_abundance(out_paths, cpm_and_proportion_path, colors,
                   max_transcripts):
    py_script_rel_path = sys.argv[0]
    py_script_abs_path = os.path.abspath(py_script_rel_path)
    script_dir = os.path.dirname(py_script_abs_path)
    r_script_path = os.path.join(script_dir, 'visualize_isoforms.R')
    transcript_colors_string = ','.join(colors)
    group_colors_string = ','.join(group_colors())
    command = [
        'Rscript', r_script_path, cpm_and_proportion_path,
        str(max_transcripts), transcript_colors_string, group_colors_string
    ]
    command.extend(out_paths)
    rmats_long_utils.run_command(command)


def combine_transcript_details(parsed_gencode_details, parsed_updated_details):
    combined = dict()
    if parsed_gencode_details:
        by_transcript = parsed_gencode_details.get('by_transcript', dict())
        for transcript, details in by_transcript.items():
            combined[transcript] = details

    if parsed_updated_details:
        by_transcript = parsed_updated_details.get('by_transcript', dict())
        for transcript, details in by_transcript.items():
            if transcript in combined:
                continue

            combined[transcript] = details

    return combined


def get_exon_length(exon):
    return (exon[1] - exon[0]) + 1


def get_intron_length(exon_1, exon_2):
    return (exon_2[0] - exon_1[1]) - 1


def get_min_and_max_coords(transcript_details, transcripts_to_plot_set):
    min_coord = None
    max_coord = None
    chr_name = None
    for transcript, details in transcript_details.items():
        if transcript not in transcripts_to_plot_set:
            continue

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

        if chr_name is None:
            chr_name = details['chr']

    return {'min': min_coord, 'max': max_coord, 'chr': chr_name}


def get_plot_coords_for_start_and_end(start, end, min_coord, max_coord,
                                      region_length, is_minus_strand):
    if start is None:
        start = min_coord

    if end is None:
        end = max_coord

    start = start - min_coord
    end = end - min_coord
    if is_minus_strand:
        start = (region_length - 1) - start
        end = (region_length - 1) - end

    return start, end


def get_plot_coords_by_transcript(transcript_details, edge_weight_by_node,
                                  min_coord, region_length):
    plot_coords_by_transcript = dict()
    weights_by_transcript = dict()
    any_minus_strand = False
    for transcript, details in transcript_details.items():
        plot_coords = list()
        weights = list()
        is_minus_strand = details['strand'] == '-'
        if is_minus_strand:
            any_minus_strand = True

        exons = details['exons']
        transcript_min_coord = exons[0][0]
        start_pos = transcript_min_coord - min_coord
        current_pos = start_pos
        for exon_i in range(len(exons) - 1):
            exon_coords = exons[exon_i]
            next_exon_coords = exons[exon_i + 1]
            if is_minus_strand:
                exon_start_node = 's{}'.format(exon_coords[1])
                exon_end_node = 'e{}'.format(exon_coords[0])
                prev_exon_end_node = 'e{}'.format(next_exon_coords[0])
                exon_weight = (edge_weight_by_node.get(
                    exon_start_node, dict()).get(exon_end_node))
                junction_weight = (edge_weight_by_node.get(
                    prev_exon_end_node, dict()).get(exon_start_node))
            else:
                exon_start_node = 's{}'.format(exon_coords[0])
                exon_end_node = 'e{}'.format(exon_coords[1])
                next_exon_start_node = 's{}'.format(next_exon_coords[0])
                exon_weight = (edge_weight_by_node.get(
                    exon_start_node, dict()).get(exon_end_node))
                junction_weight = (edge_weight_by_node.get(
                    exon_end_node, dict()).get(next_exon_start_node))

            exon_length = get_exon_length(exon_coords)
            intron_length = get_intron_length(exon_coords, next_exon_coords)
            exon_start = current_pos
            intron_start = current_pos + exon_length
            exon_end = intron_start - 1
            intron_end = exon_end + intron_length
            plot_coords.append((exon_start, exon_end))
            plot_coords.append((intron_start, intron_end))
            weights.append(exon_weight)
            weights.append(junction_weight)
            current_pos = intron_end + 1

        # last_exon
        exon_start = current_pos
        exon_coords = exons[-1]
        if is_minus_strand:
            exon_start_node = 's{}'.format(exon_coords[1])
            exon_end_node = 'e{}'.format(exon_coords[0])
            exon_weight = (edge_weight_by_node.get(exon_start_node,
                                                   dict()).get(exon_end_node))
        else:
            exon_start_node = 's{}'.format(exon_coords[0])
            exon_end_node = 'e{}'.format(exon_coords[1])
            exon_weight = (edge_weight_by_node.get(exon_start_node,
                                                   dict()).get(exon_end_node))

        exon_length = get_exon_length(exon_coords)
        exon_end = exon_start + (exon_length - 1)
        plot_coords.append((exon_start, exon_end))
        weights.append(exon_weight)
        # Flip the coordinates by subtracting from the region_length.
        # This puts the 5' end on the left.
        if is_minus_strand:
            adjusted_plot_coords = list()
            for start, end in reversed(plot_coords):
                adjusted_start = (region_length - 1) - start
                adjusted_end = (region_length - 1) - end
                adjusted_plot_coords.append((adjusted_end, adjusted_start))

            plot_coords = adjusted_plot_coords
            weights = list(reversed(weights))

        plot_coords_by_transcript[transcript] = plot_coords
        weights_by_transcript[transcript] = weights

    return {
        'coords': plot_coords_by_transcript,
        'weights': weights_by_transcript,
        'is_minus_strand': any_minus_strand
    }


def plot_transcripts(ax, transcripts_to_plot, plot_coords_by_transcript,
                     weights_by_transcript, region_length, line_space,
                     exon_height, colors, exon_edge_color, exon_line_width,
                     intron_line_color, intron_line_width, font_size):
    sj_count_details = list()
    exon_z_order = 100
    text_z_order = exon_z_order + 1
    five_prime_x_val = -0.03 * region_length
    colors_to_plot = colors[:len(transcripts_to_plot)]
    reversed_transcripts = list(reversed(transcripts_to_plot))
    reversed_colors = list(reversed(colors_to_plot))
    for transcript_i, transcript_id in enumerate(reversed_transcripts):
        color = reversed_colors[transcript_i]
        plot_coords = plot_coords_by_transcript[transcript_id]
        weights = weights_by_transcript[transcript_id]
        mid_y_val = line_space * (transcript_i + 1)
        top_y_val = mid_y_val + (exon_height / 2)
        bottom_y_val = mid_y_val - (exon_height / 2)
        text_y_val = mid_y_val - (0.1 * line_space) + exon_height
        five_prime_y_val = mid_y_val - 2
        is_exon = True
        for region_i, region in enumerate(plot_coords):
            start, end = region
            weight = weights[region_i]
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
                intron_mid = (start + end) / 2
                plt.plot(intron_x,
                         intron_y,
                         color=intron_line_color,
                         linewidth=intron_line_width)
                if weight is not None:
                    sj_count_details.append({
                        'x': intron_mid,
                        'y': mid_y_val,
                        'weight': weight,
                        'zorder': text_z_order
                    })

        plt.text(0, text_y_val, transcript_id, fontsize=font_size)
        plt.text(five_prime_x_val,
                 five_prime_y_val,
                 "5'",
                 fontsize=font_size + 1)

    return sj_count_details


def plot_start_and_end(start_coord, end_coord, start_was_none, end_was_none,
                       y_min, y_max, color, line_width, is_minus_strand):
    if start_was_none and end_was_none:
        return

    if is_minus_strand:
        start_coord, end_coord = end_coord, start_coord
        start_was_none, end_was_none = end_was_none, start_was_none

    plt.plot([start_coord, start_coord], [y_min, y_max],
             color=color,
             linestyle='--',
             linewidth=line_width)
    if not start_was_none:
        plt.plot(start_coord,
                 y_min,
                 color=color,
                 marker=matplotlib.markers.CARETRIGHTBASE)

    plt.plot([end_coord, end_coord], [y_min, y_max],
             color=color,
             linestyle='--',
             linewidth=line_width)
    if not end_was_none:
        plt.plot(end_coord,
                 y_min,
                 color=color,
                 marker=matplotlib.markers.CARETLEFTBASE)


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
                         region_length, start_coord, end_coord):
    result = {
        'region_length': region_length,
        'start_coord': start_coord,
        'end_coord': end_coord,
    }
    if intron_scaling == 1:
        return result

    coord_is_exon = get_exon_regions(plot_coords_by_transcript, region_length)
    coord_translation, num_removed_units = get_intron_coord_translation(
        coord_is_exon, region_length, intron_scaling)
    translate_coords(plot_coords_by_transcript, coord_translation)

    region_length = region_length - num_removed_units
    result['region_length'] = region_length
    if start_coord is not None:
        result['start_coord'] = coord_translation[start_coord]

    if end_coord is not None:
        result['end_coord'] = coord_translation[end_coord]

    return result


def plot_structure(out_paths, sj_count_out_paths, parsed_gencode_details,
                   parsed_updated_details, gene_name, sorted_transcripts,
                   edge_weight_by_node, colors, max_transcripts,
                   intron_scaling, start_coord, end_coord):
    transcript_details = combine_transcript_details(parsed_gencode_details,
                                                    parsed_updated_details)
    required_transcripts = set(sorted_transcripts)
    found_transcripts = set(transcript_details)
    missing_transcripts = required_transcripts.difference(found_transcripts)
    if missing_transcripts:
        raise Exception('Did not find all expected transcripts: {}.'
                        ' Missing: {}'.format(required_transcripts,
                                              missing_transcripts))

    plot_structure_with_transcripts_and_colors(
        out_paths, sj_count_out_paths, sorted_transcripts, transcript_details,
        edge_weight_by_node, colors, gene_name, max_transcripts,
        intron_scaling, start_coord, end_coord)


def plot_structure_with_transcripts_and_colors(
        out_paths, sj_count_out_paths, transcripts_to_plot, transcript_details,
        edge_weight_by_node, colors, gene_name, max_transcripts,
        intron_scaling, start_coord, end_coord):
    transcripts_to_plot = transcripts_to_plot[:max_transcripts]
    transcripts_to_plot_set = set(transcripts_to_plot)
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
    weight_font_size = 4
    min_max_result = get_min_and_max_coords(transcript_details,
                                            transcripts_to_plot_set)
    min_coord = min_max_result['min']
    max_coord = min_max_result['max']
    chr_name = min_max_result['chr']
    region_length = (max_coord - min_coord) + 1
    plot_coords_result = get_plot_coords_by_transcript(transcript_details,
                                                       edge_weight_by_node,
                                                       min_coord,
                                                       region_length)
    plot_coords_by_transcript = plot_coords_result['coords']
    weights_by_transcript = plot_coords_result['weights']
    is_minus_strand = plot_coords_result['is_minus_strand']
    start_was_none = start_coord is None
    end_was_none = end_coord is None
    start_coord, end_coord = get_plot_coords_for_start_and_end(
        start_coord, end_coord, min_coord, max_coord, region_length,
        is_minus_strand)
    scaling_result = apply_intron_scaling(plot_coords_by_transcript,
                                          intron_scaling, region_length,
                                          start_coord, end_coord)
    region_length = scaling_result['region_length']
    start_coord = scaling_result['start_coord']
    end_coord = scaling_result['end_coord']
    sj_count_details = plot_transcripts(
        ax, transcripts_to_plot, plot_coords_by_transcript,
        weights_by_transcript, region_length, line_space, exon_height, colors,
        exon_edge_color, exon_line_width, intron_line_color, intron_line_width,
        font_size)

    low_x_lim = -0.1 * region_length
    high_x_lim = 1.1 * region_length
    low_y_lim = 0
    high_y_lim = line_space * (num_transcripts + 1.5)
    bot_region_y = line_space - exon_height
    top_region_y = (line_space * num_transcripts) + exon_height
    mid_x_lim = (high_x_lim + low_x_lim) / 2
    plot_start_and_end(start_coord, end_coord, start_was_none, end_was_none,
                       bot_region_y, top_region_y, intron_line_color,
                       intron_line_width, is_minus_strand)

    title_x_val = mid_x_lim
    title_y_val = line_space * (num_transcripts + 1)
    title = '{} {}:{}-{}'.format(gene_name, chr_name, min_coord, max_coord)
    plt.text(title_x_val,
             title_y_val,
             title,
             fontsize=font_size + 1,
             horizontalalignment='center')
    ax.set_xlim(low_x_lim, high_x_lim)
    ax.set_ylim(low_y_lim, high_y_lim)
    plt.axis('off')
    # Let the plot take up the whole figure.
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0, hspace=0)
    for out_path in out_paths:
        if FORCE_775:
            with open(out_path, 'wt') as h:
                pass  # create file

            mode_775 = stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH
            os.chmod(out_path, mode_775)

        plt.savefig(out_path, dpi=300, pad_inches=0)

    if sj_count_details:
        for sj_count_out_path in sj_count_out_paths:
            for details in sj_count_details:
                plt.text(details['x'],
                         details['y'],
                         details['weight'],
                         fontsize=weight_font_size,
                         rotation=90,
                         horizontalalignment='center',
                         verticalalignment='center',
                         zorder=details['zorder'],
                         bbox={
                             'facecolor': COLORS['white'],
                             'linestyle': '',
                             'alpha': 0.5
                         })

            plt.savefig(sj_count_out_path, dpi=300, pad_inches=0)


def determine_main_transcripts_and_gene_name(gene_id, orig_gene_name,
                                             orig_main_transcript_ids,
                                             gencode_gtf_details,
                                             sig_transcripts):
    gene_name = None
    canonical_transcript = None
    if gencode_gtf_details:
        gene_name = gencode_gtf_details['gene_name']
        canonical_transcript = gencode_gtf_details['canonical_transcript']

    if orig_gene_name is not None:
        gene_name = orig_gene_name
    elif gene_name is None:
        gene_name = gene_id

    if orig_main_transcript_ids:
        return {
            'main_transcripts': orig_main_transcript_ids,
            'gene_name': gene_name,
            'canonical_transcript': canonical_transcript
        }

    main_transcript_ids = list()
    if sig_transcripts:
        main_transcript_ids.extend(sig_transcripts)
    elif canonical_transcript:
        main_transcript_ids = [canonical_transcript]

    return {
        'main_transcripts': main_transcript_ids,
        'gene_name': gene_name,
        'canonical_transcript': canonical_transcript
    }


def parse_gtf_details_for_gene_id(gene_id, gtf_path):
    details_by_transcript = dict()
    gene_name = None
    canonical_transcript = None
    if not gtf_path:
        return {
            'by_transcript': details_by_transcript,
            'gene_name': gene_name,
            'canonical_transcript': canonical_transcript
        }

    with open(gtf_path, 'rt') as handle:
        for line in handle:
            parsed = rmats_long_utils.parse_gtf_line(line)
            if not parsed:
                continue

            line_gene_id_str = parsed['attributes'].get('gene_id', '')
            line_gene_ids = line_gene_id_str.split(',')
            if gene_id not in line_gene_ids:
                continue

            line_gene_name = parsed['attributes'].get('gene_name')
            if (not gene_name) and line_gene_name:
                gene_name = line_gene_name

            transcript_id = parsed['attributes'].get('transcript_id')
            if not transcript_id:
                continue

            tags = parsed['attributes'].get('tag', list())
            if (not canonical_transcript) and ('Ensembl_canonical' in tags):
                canonical_transcript = transcript_id

            for_transcript = rmats_long_utils.try_get_or_set_default(
                details_by_transcript, transcript_id, dict())

            for_transcript['chr'] = parsed['chr']
            strand = parsed['strand']
            if strand != '.':
                for_transcript['strand'] = strand

            if parsed['feature'] != 'exon':
                continue

            exons = rmats_long_utils.try_get_or_set_default(
                for_transcript, 'exons', list())
            exons.append((parsed['start'], parsed['end']))

    for details in details_by_transcript.values():
        exons = details.get('exons')
        if not exons:
            details['exons'] = list()
        else:
            exons.sort()

    return {
        'by_transcript': details_by_transcript,
        'gene_name': gene_name,
        'canonical_transcript': canonical_transcript
    }


def get_weight_from_attributes_string(attributes_string):
    pending_chars = list()
    pending_strings = list()
    weight_string = None
    for char in attributes_string:
        if char.isspace() or char in [',', '=', '[']:
            string = ''.join(pending_chars)
            pending_chars = list()
            if string and (pending_strings == ['label', '=']):
                weight_string = string
                break

            if string == 'label':
                pending_strings = ['label']
            elif string:
                pending_strings = list()

            if char == '=':
                if pending_strings == ['label']:
                    pending_strings.append('=')
                    continue

            if not char.isspace():
                pending_strings = list()

            continue

        pending_chars.append(char)

    if not weight_string:
        return None

    if weight_string.startswith('"') and weight_string.endswith('"'):
        weight_string = weight_string[1:-1]

    try:
        return int(weight_string)
    except ValueError:
        return None


def parse_edge_weights_from_graph(gene_id, graph_file):
    edge_weight_by_node = dict()
    graph_lines = rmats_long_utils.find_gene_graph_in_graph_file(
        gene_id, graph_file)
    if not graph_lines:
        return edge_weight_by_node

    for line in graph_lines:
        tokens = line.rstrip('\n').split()
        if not tokens:
            continue

        if tokens[0] == '}':
            break

        if tokens[:2] == ['strict', 'digraph']:
            if tokens[2] != gene_id:
                raise Exception(
                    'Unexpected gene_id {} instead of {} in {}'.format(
                        tokens[2], gene_id, graph_file))

            continue

        if tokens[:2] == ['{', 'rank']:
            continue

        if tokens[0] in ['graph', 'edge']:
            continue

        nodes = list()
        attribute_i = None
        for token_i, token in enumerate(tokens):
            if token.startswith('['):
                attribute_i = token_i
                break

            if token != '->':
                nodes.append(token)

        if len(nodes) != 2:
            continue

        from_node = nodes[0]
        to_node = nodes[1]
        if attribute_i is None:
            continue

        attributes_string = ' '.join(tokens[attribute_i:])
        weight = get_weight_from_attributes_string(attributes_string)
        if weight is None:
            continue

        weight_for_node = rmats_long_utils.try_get_or_set_default(
            edge_weight_by_node, from_node, dict())
        old_weight = weight_for_node.get(to_node, 0)
        weight_for_node[to_node] = old_weight + weight

    return edge_weight_by_node


def visualize_isoforms(args):
    rmats_long_utils.create_output_dir(args.out_dir)
    primary_id = args.gene_id
    if args.is_asm:
        primary_id = args.asm_id

    parsed_gencode_details = parse_gtf_details_for_gene_id(
        primary_id, args.gencode_gtf)
    parsed_updated_details = parse_gtf_details_for_gene_id(
        primary_id, args.updated_gtf)
    edge_weight_by_node = parse_edge_weights_from_graph(
        args.gene_id, args.graph_file)
    parsed_asm_details = parse_gtf_details_for_gene_id(args.asm_id,
                                                       args.asm_gtf)
    sig_asm_transcripts = rmats_long_utils.select_significant_transcripts(
        args.asm_id, args.diff_transcripts)
    if args.is_asm:
        sig_transcripts = sig_asm_transcripts
    else:
        sig_transcripts = rmats_long_utils.select_significant_transcripts(
            args.gene_id, args.diff_transcripts)

    main_transcript_result = determine_main_transcripts_and_gene_name(
        args.gene_id, args.gene_name, args.main_transcript_ids,
        parsed_gencode_details, sig_transcripts)
    main_transcript_ids = main_transcript_result['main_transcripts']
    gene_name = main_transcript_result['gene_name']
    canonical_transcript = main_transcript_result['canonical_transcript']

    cpm_and_proportion_file = os.path.join(
        args.out_dir, '{}_cpm_and_proportion.tsv'.format(primary_id))
    abundance_plot_file_base = os.path.join(args.out_dir,
                                            '{}_abundance'.format(primary_id))
    structure_plot_file_base = os.path.join(args.out_dir,
                                            '{}_structure'.format(primary_id))
    plot_file_types = [args.plot_file_type]
    if args.plot_file_type == 'all':
        plot_file_types = ['.pdf', '.png']

    abundance_plot_files = list()
    structure_plot_files = list()
    sj_count_plot_files = list()
    for plot_file_type in plot_file_types:
        abundance_plot_files.append('{}{}'.format(abundance_plot_file_base,
                                                  plot_file_type))
        structure_plot_files.append('{}{}'.format(structure_plot_file_base,
                                                  plot_file_type))
        sj_count_plot_files.append('{}_sj_counts{}'.format(
            structure_plot_file_base, plot_file_type))

    if not args.abundance:
        sorted_transcripts = sort_transcripts_using_structure(
            main_transcript_ids, canonical_transcript, sig_asm_transcripts,
            parsed_gencode_details, parsed_updated_details, parsed_asm_details)
        colors = get_transcript_colors(sorted_transcripts, main_transcript_ids,
                                       args.max_transcripts)
        plot_structure(structure_plot_files, sj_count_plot_files,
                       parsed_gencode_details, parsed_updated_details,
                       gene_name, sorted_transcripts, edge_weight_by_node,
                       colors, args.max_transcripts, args.intron_scaling,
                       args.start_coord, args.end_coord)
        return

    sorted_transcripts = calc_proportion_and_cpm(
        args.abundance, primary_id, main_transcript_ids, canonical_transcript,
        args.group_1, args.group_2, args.group_1_name, args.group_2_name,
        args.is_asm, cpm_and_proportion_file)
    colors = get_transcript_colors(sorted_transcripts, main_transcript_ids,
                                   args.max_transcripts)
    plot_abundance(abundance_plot_files, cpm_and_proportion_file, colors,
                   args.max_transcripts)
    plot_structure(structure_plot_files, sj_count_plot_files,
                   parsed_gencode_details, parsed_updated_details, gene_name,
                   sorted_transcripts, edge_weight_by_node, colors,
                   args.max_transcripts, args.intron_scaling, args.start_coord,
                   args.end_coord)


def main():
    args = parse_args()
    visualize_isoforms(args)


if __name__ == '__main__':
    main()
