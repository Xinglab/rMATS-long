import argparse
import math

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches

from rmats_long import rmats_long_utils

COLORS = {
    'black': '#000000',
    'blue': '#0000FF',
    'red': '#FF0000',
}


def parse_args():
    parser = argparse.ArgumentParser(
        description='Output a diagram of a splicing graph')
    parser.add_argument(
        '--gtf',
        required=True,
        help='The path to a .gtf file with isoform definitions')
    parser.add_argument('--out-file',
                        required=True,
                        help='The file for the output plot')
    parser.add_argument(
        '--gene-id', help='Only plot isoforms with this "gene_id" attribute')
    parser.add_argument('--asm-id',
                        help='Only plot isoforms with this "asm_id" attribute')
    parser.add_argument(
        '--transcript-ids',
        help=('A comma separated list of isoform IDs. Only plot isoforms with'
              ' one of these "transcript_id" attributes.'))
    parser.add_argument(
        '--intron-scaling',
        type=float,
        default=1,
        help=('The factor to use to reduce intron length in the plot.'
              ' A value of 2 would reduce introns to 1/2 of the'
              ' original plot length. (default %(default)s)'))
    parser.add_argument(
        '--equal-spacing',
        action='store_true',
        help=('Plot each splice site a fixed distance from the previous'
              ' splice site'))
    parser.add_argument(
        '--junction-counts',
        help=('A tab separated file with headers:'
              ' ["chr", "start", "end", "count"].'
              ' Any junction plotted between "start" and "end" will have'
              ' "count" displayed.'))
    parser.add_argument(
        '--color-by-isoform',
        action='store_true',
        help=('Use a color for each isoform and apply that color to'
              ' exons and junctions only used by that isoform'))
    parser.add_argument('--color-first-transcript-only',
                        action='store_true',
                        help=('Color all exons and junctions for the'
                              ' 1st value from --transcript-ids'))
    parser.add_argument(
        '--transcript-colors',
        help=('A tab separated file with headers:'
              ' ["transcript_id", "color"].'
              'The color column is an RGB hex string (#RRGGBB)'))
    parser.add_argument(
        '--ids-from-color-file',
        action='store_true',
        help='set --transcripts-ids from the transcript_id column of'
        ' --transcript-colors')
    parser.add_argument('--show-isoform-endpoint-symbols',
                        action='store_true',
                        help=("Add a symbol at each isoform's start"
                              " and a symbol at each isoform's end"))
    parser.add_argument(
        '--font-size',
        type=int,
        default=8,
        help='The font size for the plot (default %(default)s)')
    parser.add_argument('--plot-width',
                        type=float,
                        default=12,
                        help='The plot width in inches (default %(default)s)')
    parser.add_argument('--plot-height',
                        type=float,
                        default=6,
                        help='The plot height in inches (default %(default)s)')
    parser.add_argument(
        '--plot-dpi',
        type=float,
        default=300,
        help='The plot resolution in dots per inch (default %(default)s)')

    args = parser.parse_args()
    if args.transcript_ids:
        args.transcript_ids = args.transcript_ids.split(',')
    else:
        args.transcript_ids = list()

    return args


def validate_color_string(string):
    expected_chars = '0123456789abcdef'
    string = string.lower()
    if len(string) != 7:
        return None

    if string[0] != '#':
        return None

    for char in string[1:]:
        if char not in expected_chars:
            return None

    return string


def load_transcript_colors(transcript_colors_path):
    supplied_color_by_transcript = dict()
    transcripts = list()
    result = {
        'by_transcript': supplied_color_by_transcript,
        'transcripts': transcripts
    }
    if not transcript_colors_path:
        return result

    with open(transcript_colors_path, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            transcript_id = row['transcript_id']
            color = row['color']
            validated = validate_color_string(color)
            if not validated:
                raise Exception('Invalid color string: {}'.format(color))

            supplied_color_by_transcript[transcript_id] = validated
            transcripts.append(transcript_id)

    return result


def get_text_sizes(font_size, all_coords, isoform_defs):
    plot_width = 1
    plot_height = 1
    fig = plt.figure(figsize=(plot_width, plot_height))
    data_to_fig = fig.transFigure.inverted()

    text_x = 0
    text_y = 0
    text_by_category = dict()

    max_coord = max(all_coords)
    formatted_coord = format_coord(max_coord)
    text = fig.text(text_x, text_y, formatted_coord, fontsize=font_size)
    text_by_category['coord'] = [text]

    five_prime = five_prime_string()
    text = fig.text(text_x, text_y, five_prime, fontsize=font_size)
    text_by_category['five_prime'] = [text]

    isoform_texts = list()
    text_by_category['isoform_id'] = isoform_texts
    for isoform_id in isoform_defs.keys():
        text = fig.text(text_x, text_y, isoform_id, fontsize=font_size)
        isoform_texts.append(text)

    fig.draw_without_rendering()
    sizes = dict()
    for category, texts in text_by_category.items():
        max_width = 0
        max_height = 0
        for text in texts:
            bbox = text.get_window_extent()
            points = bbox.get_points()
            x0 = points[0][0]
            y0 = points[0][1]
            x1 = points[1][0]
            y1 = points[1][1]
            width = x1 - x0
            height = y1 - y0
            transformed = data_to_fig.transform([width, height])
            width = transformed[0]
            height = transformed[1]
            max_width = max(max_width, width)
            max_height = max(max_height, height)

        sizes[category] = {'width': max_width, 'height': max_height}

    plt.close(fig)
    points_per_inch = 72
    font_inches = font_size / points_per_inch
    sizes['font'] = {'width': font_inches, 'height': font_inches}
    return sizes


def sort_loaded_exons(by_isoform):
    for isoform_details in by_isoform.values():
        isoform_details['exons'].sort()


def check_that_transcripts_were_found(by_isoform, transcript_ids):
    for transcript_id in transcript_ids:
        if transcript_id not in by_isoform:
            print('WARNING: transcript ID {} requested but not found'.format(
                transcript_id))


def load_isoform_definitions(gtf, gene_id, asm_id, transcript_ids):
    by_isoform = dict()
    used_chr = None
    with open(gtf, 'rt') as handle:
        for line in handle:
            parsed = rmats_long_utils.parse_gtf_line(line)
            if not parsed:
                continue

            feature = parsed['feature']
            if feature != 'exon':
                continue

            chr_name = parsed['chr']
            start = parsed['start']
            end = parsed['end']
            strand = parsed['strand']
            attributes = parsed['attributes']
            found_gene_id = attributes.get('gene_id')
            if gene_id and (found_gene_id != gene_id):
                continue

            found_asm_id = attributes.get('asm_id')
            if asm_id and (found_asm_id != asm_id):
                continue

            found_transcript_id = attributes.get('transcript_id')
            if transcript_ids and (found_transcript_id not in transcript_ids):
                continue

            if used_chr is None:
                used_chr = chr_name
            elif chr_name != used_chr:
                raise Exception(
                    'Cannot plot for multiple chr values: {} {}'.format(
                        used_chr, chr_name))

            isoform_details = by_isoform.get(found_transcript_id)
            if not isoform_details:
                isoform_details = {
                    'exons': list(),
                    'strand': strand,
                }
                by_isoform[found_transcript_id] = isoform_details

            isoform_details['exons'].append((start, end))

    sort_loaded_exons(by_isoform)
    check_that_transcripts_were_found(by_isoform, transcript_ids)

    return {
        'defs': by_isoform,
        'chr': used_chr,
    }


def record_isoform_start_and_end(isoform, strand, exons, isoforms_by_start,
                                 isoforms_by_end):
    start = exons[0][0]
    end = exons[-1][-1]
    if strand == '-':
        start, end = end, start

    start_isoforms = rmats_long_utils.try_get_or_set_default(
        isoforms_by_start, start, set())
    start_isoforms.add(isoform)
    end_isoforms = rmats_long_utils.try_get_or_set_default(
        isoforms_by_end, end, set())
    end_isoforms.add(isoform)


def get_region_details(isoform_defs):
    junctions = set()
    isoforms_by_junction = dict()
    isoforms_by_start = dict()
    isoforms_by_end = dict()
    all_coords = set()
    strand = '+'
    for isoform, isoform_details in isoform_defs.items():
        exons = isoform_details['exons']
        strand = isoform_details['strand']
        prev_end = None
        for start, end in exons:
            all_coords.add(start)
            all_coords.add(end)
            if prev_end is not None:
                junction = (prev_end, start)
                junctions.add(junction)
                isoforms = rmats_long_utils.try_get_or_set_default(
                    isoforms_by_junction, junction, set())
                isoforms.add(isoform)

            prev_end = end

        record_isoform_start_and_end(isoform, strand, exons, isoforms_by_start,
                                     isoforms_by_end)

    junctions = sorted(junctions)
    all_coords = sorted(all_coords)
    coord_to_i = dict()
    is_exon_by_i = list()
    isoforms_by_i = list()
    for coord_i, coord in enumerate(all_coords):
        coord_to_i[coord] = coord_i
        is_exon_by_i.append(False)
        isoforms_by_i.append(set())

    for isoform, isoform_details in isoform_defs.items():
        exons = isoform_details['exons']
        for start, end in exons:
            start_i = coord_to_i[start]
            end_i = coord_to_i[end]
            for exon_i in range(start_i, end_i):
                is_exon_by_i[exon_i] = True
                isoforms_by_i[exon_i].add(isoform)

    return {
        'all_coords': all_coords,
        'is_exon_by_i': is_exon_by_i,
        'isoforms_by_i': isoforms_by_i,
        'junctions': junctions,
        'isoforms_by_junction': isoforms_by_junction,
        'isoforms_by_start': isoforms_by_start,
        'isoforms_by_end': isoforms_by_end,
        'strand': strand,
    }


def load_junction_counts(junction_counts_path, chr_name, junctions):
    loaded = dict()
    if not junction_counts_path:
        return loaded

    with open(junction_counts_path, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            found_chr = row['chr']
            if found_chr != chr_name:
                continue

            start = int(row['start'])
            end = int(row['end'])
            pair = [start, end]
            pair.sort()
            pair = tuple(pair)
            if pair not in junctions:
                continue

            count = row['count']
            loaded[pair] = count

    return loaded


def update_plot_coords_for_minus_strand(coord_to_plot_x):
    updated = dict()
    min_x = min(coord_to_plot_x.values())
    max_x = max(coord_to_plot_x.values())
    for coord, plot_x in coord_to_plot_x.items():
        from_end = max_x - plot_x
        new_x = min_x + from_end
        updated[coord] = new_x

    return updated


def determine_coord_spacing(all_coords, is_exon_by_i, strand, intron_scaling,
                            equal_spacing):
    coord_to_plot_x = dict()
    if equal_spacing:
        for coord_i, coord in enumerate(all_coords):
            coord_to_plot_x[coord] = coord_i
    else:
        scale_mult = 1 / intron_scaling
        prev_plot_x = None
        last_coord_i = len(all_coords) - 1
        for coord_i, coord in enumerate(all_coords):
            if coord_i == 0:
                coord_to_plot_x[coord] = coord
                prev_plot_x = coord

            if coord_i == last_coord_i:
                break

            next_coord = all_coords[coord_i + 1]
            to_next = next_coord - coord
            is_exon = is_exon_by_i[coord_i]
            if is_exon:
                plot_x_interval = to_next
            else:
                scaled = scale_mult * to_next
                plot_x_interval = max(1, scaled)

            plot_x = prev_plot_x + plot_x_interval
            coord_to_plot_x[next_coord] = plot_x
            prev_plot_x = plot_x

    if strand == '-':
        coord_to_plot_x = update_plot_coords_for_minus_strand(coord_to_plot_x)

    return coord_to_plot_x


def get_isoform_color(supplied_color_by_transcript, isoform_id, isoform_i,
                      num_isoforms, color_by_isoform):
    if not color_by_isoform:
        return COLORS['black']

    if supplied_color_by_transcript:
        color = supplied_color_by_transcript.get(isoform_id)
        if color is None:
            return COLORS['black']

        return color

    for_two = [COLORS['red'], COLORS['blue']]
    if num_isoforms <= 2:
        return for_two[isoform_i]

    dark2 = matplotlib.color_sequences['Dark2']
    if num_isoforms <= len(dark2):
        return dark2[isoform_i]

    # Use equally spaced grey colors (not white)
    percent_black = (isoform_i + 1) / (num_isoforms + 1)
    return str(percent_black)


def get_isoform_colors(isoform_defs, supplied_color_by_transcript,
                       transcript_ids, color_by_isoform,
                       color_first_transcript_only):
    isoform_colors = dict()
    isoform_ids = list()
    for transcript_id in transcript_ids:
        if transcript_id in isoform_defs:
            isoform_ids.append(transcript_id)

    sorted_keys = sorted(isoform_defs.keys())
    for isoform in sorted_keys:
        if isoform not in isoform_ids:
            isoform_ids.append(isoform)

    if color_first_transcript_only:
        first_id = isoform_ids[0]
        color = get_isoform_color(supplied_color_by_transcript, first_id, 0, 1,
                                  True)
        isoform_colors[first_id] = color
        return {'colors': isoform_colors, 'order': isoform_ids}

    num_isoforms = len(isoform_ids)
    for isoform_i, isoform in enumerate(isoform_ids):
        color = get_isoform_color(supplied_color_by_transcript, isoform,
                                  isoform_i, num_isoforms, color_by_isoform)
        isoform_colors[isoform] = color

    return {'colors': isoform_colors, 'order': isoform_ids}


def get_exon_color(coord_i, isoforms_by_i, isoform_colors,
                   color_first_transcript_only):
    coord_isoforms = isoforms_by_i[coord_i]
    if color_first_transcript_only:
        for isoform in coord_isoforms:
            color = isoform_colors.get(isoform)
            if color:
                return color

        return COLORS['black']

    if len(coord_isoforms) == 1:
        isoform = list(coord_isoforms)[0]
        return isoform_colors[isoform]

    return COLORS['black']


def get_junction_color(junction, isoforms_by_junction, isoform_colors,
                       color_first_transcript_only):
    junction_isoforms = isoforms_by_junction[junction]
    if color_first_transcript_only:
        for isoform in junction_isoforms:
            color = isoform_colors.get(isoform)
            if color:
                return color

        return COLORS['black']

    if len(junction_isoforms) == 1:
        isoform = list(junction_isoforms)[0]
        return isoform_colors[isoform]

    return COLORS['black']


def draw_rectangle(x1, y1, x2, y2, line_width, ax, face_color):
    anchor = (x1, y1)
    width = x2 - x1
    height = y2 - y1
    edge_color = COLORS['black']
    fill = face_color != edge_color
    rect = matplotlib.patches.Rectangle(anchor,
                                        width,
                                        height,
                                        linewidth=line_width,
                                        fill=fill,
                                        edgecolor=edge_color,
                                        facecolor=face_color)
    ax.add_patch(rect)


def draw_line(x1, y1, x2, y2, line_width, color):
    plt.plot([x1, x2], [y1, y2], linewidth=line_width, color=color)


def draw_vee(x1, y1, x2, y2, x3, y3, line_width, color):
    plt.plot([x1, x2, x3], [y1, y2, y3], linewidth=line_width, color=color)


def group_overlapping_coords(coords_and_fig_xs, coord_x_size):
    details_by_coord = dict()
    coord_groups = list()
    prev_x = None
    for coord, fig_x in coords_and_fig_xs:
        details_by_coord[coord] = {
            'orig_x': fig_x,
            'fig_x': fig_x,
            'alpha': 1,
        }

        if prev_x is None:
            new_group = True
        else:
            diff_x = fig_x - prev_x
            new_group = diff_x >= coord_x_size

        if new_group:
            coord_groups.append([coord])
        else:
            coord_groups[-1].append(coord)

        prev_x = fig_x

    return {'details_by_coord': details_by_coord, 'coord_groups': coord_groups}


def update_x_and_alpha_for_groups(coord_groups, coord_x_size,
                                  details_by_coord):
    half_coord_x_size = coord_x_size / 2
    for group in coord_groups:
        num_coords = len(group)
        if num_coords == 1:
            continue

        xs = list()
        for coord in group:
            details = details_by_coord[coord]
            orig_x = details['orig_x']
            xs.append(orig_x)

        min_x = min(xs)
        max_x = max(xs)
        x_len = max_x - min_x
        required_len = coord_x_size * (num_coords - 1)
        # Try to fit with an extra half width on each side
        extra_width = half_coord_x_size
        with_extra_width = x_len + (extra_width * 2)
        if with_extra_width < required_len:
            # Extend up to an extra width on each side
            extra_width = coord_x_size
            with_extra_width = x_len + (extra_width * 2)

        new_x = min_x - extra_width
        interval = with_extra_width / (num_coords - 1)
        density = required_len / with_extra_width
        alpha = 1 / math.ceil(density)
        for coord in group:
            details_by_coord[coord]['fig_x'] = new_x
            details_by_coord[coord]['alpha'] = alpha
            new_x += interval


def draw_coordinate_annotations(coord_to_plot_x, sg_fig_pos, isoforms_by_start,
                                isoforms_by_end, isoform_order, isoform_colors,
                                plot_width, plot_height,
                                show_isoform_endpoint_symbols, text_sizes,
                                font_size, fig, ax):
    coord_x_size = text_sizes['coord']['height'] / plot_width
    font_width = text_sizes['font']['width']
    sg_margin = font_width / plot_height
    coord_top_y = sg_fig_pos['bottom'] - sg_margin
    coord_y_size = text_sizes['coord']['width']
    coord_y_fig = coord_y_size / plot_height
    coord_bot_y = coord_top_y - coord_y_fig
    coords_and_fig_xs = list()
    for coord, plot_x in coord_to_plot_x.items():
        fig_x = convert_data_to_fig([(plot_x, 0)], fig, ax)[0][0]
        coords_and_fig_xs.append((coord, fig_x))
        if show_isoform_endpoint_symbols:
            add_isoform_start_and_end_symbols(coord, isoforms_by_start,
                                              isoforms_by_end, isoform_order,
                                              isoform_colors, plot_width,
                                              plot_height, text_sizes, fig_x,
                                              coord_bot_y, fig)

    coords_and_fig_xs.sort(key=lambda pair: pair[1])
    group_result = group_overlapping_coords(coords_and_fig_xs, coord_x_size)
    details_by_coord = group_result['details_by_coord']
    coord_groups = group_result['coord_groups']
    update_x_and_alpha_for_groups(coord_groups, coord_x_size, details_by_coord)
    for coord, details in details_by_coord.items():
        fig_x = details['fig_x']
        alpha = details['alpha']
        formatted_coord = format_coord(coord)
        fig.text(fig_x,
                 coord_top_y,
                 formatted_coord,
                 fontsize=font_size,
                 horizontalalignment='right',
                 verticalalignment='center',
                 rotation=90,
                 rotation_mode='anchor',
                 alpha=alpha)


def format_coord(coord):
    comma_separated = '{:,}'.format(coord)
    return comma_separated


def add_isoform_symbol(is_end, color, plot_width, plot_height, text_sizes,
                       fig_x, fig_y, fig):
    half_font = text_sizes['font']['width'] / 2
    left_x = fig_x
    top_y = fig_y
    symbol_x_size = half_font
    symbol_y_size = half_font
    width = symbol_x_size / plot_width
    height = symbol_y_size / plot_height
    if is_end:
        other_x = left_x - width
    else:
        other_x = left_x + width

    bot_y = top_y - height
    mid_y = (top_y + bot_y) / 2
    points = list()
    points.append((left_x, top_y))
    points.append((left_x, bot_y))
    points.append((other_x, mid_y))
    rect = matplotlib.patches.Polygon(points,
                                      fill=True,
                                      edgecolor=None,
                                      facecolor=color)
    fig.add_artist(rect)
    return bot_y


def add_isoform_start_symbol(color, plot_width, plot_height, text_sizes, fig_x,
                             fig_y, fig):
    is_end = False
    return add_isoform_symbol(is_end, color, plot_width, plot_height,
                              text_sizes, fig_x, fig_y, fig)


def add_isoform_end_symbol(color, plot_width, plot_height, font_size, fig_x,
                           fig_y, fig):
    is_end = True
    return add_isoform_symbol(is_end, color, plot_width, plot_height,
                              font_size, fig_x, fig_y, fig)


def write_count(x, y, text, font_size):
    if not text:
        return

    plt.text(x,
             y,
             text,
             fontsize=font_size,
             horizontalalignment='center',
             verticalalignment='bottom')


def five_prime_string():
    return "5'"


def write_5_prime(sg_fig_pos, plot_width, text_sizes, font_size, fig):
    text = five_prime_string()
    text_x_size = text_sizes['five_prime']['width']
    text_x_fig = text_x_size / plot_width
    font_width = text_sizes['font']['width']
    sg_margin = font_width / plot_width
    fig_x = sg_fig_pos['left'] - (text_x_fig + sg_margin)
    fig_y = (sg_fig_pos['bottom'] + sg_fig_pos['top']) / 2
    fig.text(fig_x, fig_y, text, fontsize=font_size)


def prepare_legend(isoform_colors, isoform_order, text_sizes):
    handles = list()
    for isoform in isoform_order:
        color = isoform_colors.get(isoform)
        if color is None:
            color = COLORS['black']

        handle = matplotlib.patches.Patch(color=color, label=isoform)
        handles.append(handle)

    if not handles:
        return {
            'handles': handles,
            'x_size': 0,
            'y_size': 0,
        }

    num_isoforms = len(handles)
    num_spaces = num_isoforms + 1
    line_height = text_sizes['isoform_id']['height']
    legend_y_size = (num_isoforms + (num_spaces / 2)) * line_height
    font_width = text_sizes['font']['width']
    icon_x_size = 3 * font_width
    x_spacing = 2 * font_width
    isoform_x_size = text_sizes['isoform_id']['width']
    legend_x_size = x_spacing + icon_x_size + isoform_x_size
    return {
        'handles': handles,
        'x_size': legend_x_size,
        'y_size': legend_y_size,
    }


def convert_nested_array_to_list_of_tuple(arrays):
    result = list()
    for nested in arrays:
        result.append(tuple(nested))

    return result


def convert_fig_to_ax_coords(coords, fig, ax):
    from_fig = fig.transFigure.transform
    to_ax = ax.transAxes.inverted().transform
    new_coords = to_ax(from_fig(coords))
    result = convert_nested_array_to_list_of_tuple(new_coords)
    return result


def convert_data_to_fig(coords, fig, ax):
    from_data = ax.transData.transform
    to_fig = fig.transFigure.inverted().transform
    new_coords = to_fig(from_data(coords))
    result = convert_nested_array_to_list_of_tuple(new_coords)
    return result


def convert_points_to_data(points, plot_width, plot_height, fig, ax):
    points_per_inch = 72
    inches = points / points_per_inch
    fig_x = inches / plot_width
    fig_y = inches / plot_height
    from_fig = fig.transFigure.transform
    to_data = ax.transData.inverted().transform
    coords = [(0, 0), (fig_x, fig_y)]
    new_coords = to_data(from_fig(coords))
    converted = convert_nested_array_to_list_of_tuple(new_coords)
    data_width = converted[1][0] - converted[0][0]
    data_height = converted[1][1] - converted[0][1]
    return (data_width, data_height)


def add_legend(legend_handles, legend_y_size, sg_fig_pos, plot_width,
               plot_height, font_size, text_sizes, fig, ax):
    font_width = text_sizes['font']['width']
    sg_margin = font_width / plot_width
    legend_left = sg_fig_pos['right'] + sg_margin
    graph_bottom = sg_fig_pos['bottom']
    graph_top = sg_fig_pos['top']
    graph_mid = (graph_bottom + graph_top) / 2
    legend_y_size_fig = legend_y_size / plot_height
    legend_bottom = graph_mid - (legend_y_size_fig / 2)
    loc_fig = (legend_left, legend_bottom)
    loc_ax = convert_fig_to_ax_coords([loc_fig], fig, ax)[0]
    ax.legend(handles=legend_handles, loc=loc_ax, fontsize=font_size)


def save_plot(plot_dpi, out_file):
    plt.savefig(out_file, dpi=plot_dpi)


# The splice graph is one part of the figure.
# The space needed for the other parts depends on
# the font size and the isoforms.
# Use the rest of the space for the splice graph.
#
#    |----------------| |----------|
# 5' |--splice graph--| |--legend--|
#    |----------------| |----------|
#        coordinates
def set_splice_graph_space_within_fig(plot_width, plot_height,
                                      show_isoform_endpoint_symbols,
                                      legend_x_size, text_sizes):
    font_width = text_sizes['font']['width']
    five_prime_x_size = text_sizes['five_prime']['width']
    sg_margin = font_width
    outer_margin = font_width
    coordinate_y_size = text_sizes['coord']['width']
    symbols_y_size = 0
    if show_isoform_endpoint_symbols:
        symbols_y_size = font_width * 4

    left_size = sg_margin + outer_margin + five_prime_x_size
    after_right_size = sg_margin + outer_margin + legend_x_size
    bottom_size = sg_margin + outer_margin + coordinate_y_size + symbols_y_size
    after_top_size = sg_margin + outer_margin
    left_fig = left_size / plot_width
    after_right_fig = after_right_size / plot_width
    bottom_fig = bottom_size / plot_height
    after_top_fig = after_top_size / plot_height
    right_fig = 1 - after_right_fig
    top_fig = 1 - after_top_fig
    plt.subplots_adjust(left=left_fig,
                        right=right_fig,
                        bottom=bottom_fig,
                        top=top_fig)
    return {
        'left': left_fig,
        'right': right_fig,
        'bottom': bottom_fig,
        'top': top_fig,
    }


def set_axis_limits(coord_to_plot_x, exon_bottom, vee_top, exon_line_width,
                    junction_line_width, plot_width, plot_height, fig, ax):
    min_x = min(coord_to_plot_x.values())
    max_x = max(coord_to_plot_x.values())
    min_y = exon_bottom
    max_y = vee_top
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    max_line_width = max(exon_line_width, junction_line_width)
    data_line_sizes = convert_points_to_data(max_line_width, plot_width,
                                             plot_height, fig, ax)
    extra_x = data_line_sizes[0]
    extra_y = data_line_sizes[1]
    ax.set_xlim(min_x - extra_x, max_x + extra_x)
    ax.set_ylim(min_y - extra_y, max_y + extra_y)


def add_isoform_start_and_end_symbols(coord, isoforms_by_start,
                                      isoforms_by_end, isoform_order,
                                      isoform_colors, plot_width, plot_height,
                                      text_sizes, fig_x, fig_y, fig):
    font_y_fig = text_sizes['font']['width'] / plot_height
    symbol_y = fig_y - font_y_fig
    start_isoforms = isoforms_by_start.get(coord, list())
    end_isoforms = isoforms_by_end.get(coord, list())
    for isoform in isoform_order:
        color = isoform_colors.get(isoform)
        if color is None:
            color = COLORS['black']

        if isoform in start_isoforms:
            symbol_y = add_isoform_start_symbol(color, plot_width, plot_height,
                                                text_sizes, fig_x, symbol_y,
                                                fig)
        if isoform in end_isoforms:
            symbol_y = add_isoform_end_symbol(color, plot_width, plot_height,
                                              text_sizes, fig_x, symbol_y, fig)


def make_plot(all_coords, is_exon_by_i, isoforms_by_i, junctions,
              isoforms_by_junction, isoforms_by_start, isoforms_by_end,
              isoform_defs, strand, counts, supplied_color_by_transcript,
              transcript_ids, intron_scaling, equal_spacing, color_by_isoform,
              color_first_transcript_only, show_isoform_endpoint_symbols,
              font_size, plot_width, plot_height, plot_dpi, out_file):
    text_sizes = get_text_sizes(font_size, all_coords, isoform_defs)
    exon_mid = 0
    exon_height = 1
    vee_height = 1
    exon_top = exon_mid + exon_height
    exon_bottom = exon_mid - exon_height
    vee_top = exon_top + vee_height
    exon_line_width = 1
    junction_line_width = 1

    fig = plt.figure(figsize=(plot_width, plot_height))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_axis_off()
    isoform_colors_result = get_isoform_colors(isoform_defs,
                                               supplied_color_by_transcript,
                                               transcript_ids,
                                               color_by_isoform,
                                               color_first_transcript_only)
    isoform_colors = isoform_colors_result['colors']
    isoform_order = isoform_colors_result['order']
    legend_details = prepare_legend(isoform_colors, isoform_order, text_sizes)
    legend_x_size = legend_details['x_size']
    legend_y_size = legend_details['y_size']
    legend_handles = legend_details['handles']
    sg_fig_pos = set_splice_graph_space_within_fig(
        plot_width, plot_height, show_isoform_endpoint_symbols, legend_x_size,
        text_sizes)
    add_legend(legend_handles, legend_y_size, sg_fig_pos, plot_width,
               plot_height, font_size, text_sizes, fig, ax)
    write_5_prime(sg_fig_pos, plot_width, text_sizes, font_size, fig)

    coord_to_plot_x = determine_coord_spacing(all_coords, is_exon_by_i, strand,
                                              intron_scaling, equal_spacing)
    flat_junctions = set()
    last_coord_i = len(all_coords) - 1
    for coord_i, coord in enumerate(all_coords):
        plot_x = coord_to_plot_x[coord]
        if coord_i == last_coord_i:
            break

        next_coord = all_coords[coord_i + 1]
        next_plot_x = coord_to_plot_x[next_coord]
        is_exon = is_exon_by_i[coord_i]
        next_is_exon = is_exon_by_i[coord_i + 1]
        if is_exon:
            exon_color = get_exon_color(coord_i, isoforms_by_i, isoform_colors,
                                        color_first_transcript_only)
            draw_rectangle(plot_x, exon_bottom, next_plot_x, exon_top,
                           exon_line_width, ax, exon_color)
        elif next_is_exon:
            flat_junctions.add((coord, next_coord))

    for junction in junctions:
        start, end = junction
        start_x = coord_to_plot_x[start]
        end_x = coord_to_plot_x[end]
        mid_x = (start_x + end_x) / 2
        count = counts.get(junction, '')
        junction_color = get_junction_color(junction, isoforms_by_junction,
                                            isoform_colors,
                                            color_first_transcript_only)
        if junction in flat_junctions:
            draw_line(start_x, exon_mid, end_x, exon_mid, junction_line_width,
                      junction_color)
            write_count(mid_x, exon_mid, count, font_size)
        else:
            draw_vee(start_x, exon_top, mid_x, vee_top, end_x, exon_top,
                     junction_line_width, junction_color)
            write_count(mid_x, vee_top, count, font_size)

    # Coordinates are shown outside the ax, but rely on the x-axis
    # being finalized to convert the position from ax to fig.
    set_axis_limits(coord_to_plot_x, exon_bottom, vee_top, exon_line_width,
                    junction_line_width, plot_width, plot_height, fig, ax)
    draw_coordinate_annotations(coord_to_plot_x, sg_fig_pos, isoforms_by_start,
                                isoforms_by_end, isoform_order, isoform_colors,
                                plot_width, plot_height,
                                show_isoform_endpoint_symbols, text_sizes,
                                font_size, fig, ax)
    # TODO try 2 or 3 vee heights to reduce intersections:
    #    * plot junctions in sorted order starting with shortest distance.
    #    * for each junction try all heights and pick lowest with least intersections
    #    * https://www.geeksforgeeks.org/dsa/check-if-two-given-line-segments-intersect/
    save_plot(plot_dpi, out_file)
    plt.close(fig)


def main():
    # Select a non-interactive backend
    matplotlib.use('Agg')

    args = parse_args()
    color_result = load_transcript_colors(args.transcript_colors)
    supplied_color_by_transcript = color_result['by_transcript']
    if args.ids_from_color_file:
        args.transcript_ids = color_result['transcripts']

    if args.color_first_transcript_only and (not args.transcript_ids):
        raise Exception(
            '--color-first-transcript-only requires --transcript-ids')

    load_defs_result = load_isoform_definitions(args.gtf, args.gene_id,
                                                args.asm_id,
                                                args.transcript_ids)
    isoform_defs = load_defs_result['defs']
    chr_name = load_defs_result['chr']
    region_details = get_region_details(isoform_defs)
    all_coords = region_details['all_coords']
    is_exon_by_i = region_details['is_exon_by_i']
    isoforms_by_i = region_details['isoforms_by_i']
    junctions = region_details['junctions']
    isoforms_by_junction = region_details['isoforms_by_junction']
    isoforms_by_start = region_details['isoforms_by_start']
    isoforms_by_end = region_details['isoforms_by_end']
    strand = region_details['strand']
    counts = load_junction_counts(args.junction_counts, chr_name, junctions)
    make_plot(all_coords, is_exon_by_i, isoforms_by_i, junctions,
              isoforms_by_junction, isoforms_by_start, isoforms_by_end,
              isoform_defs, strand, counts, supplied_color_by_transcript,
              args.transcript_ids, args.intron_scaling, args.equal_spacing,
              args.color_by_isoform, args.color_first_transcript_only,
              args.show_isoform_endpoint_symbols, args.font_size,
              args.plot_width, args.plot_height, args.plot_dpi, args.out_file)

    print('plot_simple_splice_graph.py finished')


if __name__ == '__main__':
    main()
