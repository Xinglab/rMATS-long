import argparse
import os
import os.path

import pydot

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Plot the splice graph for a gene or ASM')
    parser.add_argument(
        '--event-dir',
        required=True,
        help='The output directory from detect_splicing_events.py')
    parser.add_argument('--chr',
                        required=False,
                        help='The chromosome name with the gene or ASM')
    parser.add_argument('--gene-id',
                        required=False,
                        help='The gene ID to plot')
    parser.add_argument('--asm-id', required=False, help='The ASM ID to plot')
    parser.add_argument('--out-file',
                        required=True,
                        help='The file for the output plot')
    parser.add_argument(
        '--edge-label-lines',
        action='store_true',
        help='draw a line connecting each edge label to its edge')
    parser.add_argument('--min-edge-weight',
                        type=int,
                        help='only draw edges with at least this weight')
    parser.add_argument(
        '--show-extra-nodes',
        action='store_true',
        help='include nodes other than the main annotated coordinates')

    args = parser.parse_args()
    if (args.gene_id is None) and (args.asm_id is None):
        parser.error('At least one of --gene-id or --asm-id is required')

    return args


def find_asm_id_in_event_file(asm_id, path):
    gene_id = None
    with open(path, 'rt') as handle:
        for row in rmats_long_utils.row_iterator_for_tsv_with_header(handle):
            if row['gene_id']:
                gene_id = row['gene_id']

            if row['asm_id'] == asm_id:
                row['gene_id'] = gene_id
                isoforms = rmats_long_utils.parse_isoforms_str(row['isoforms'])
                row['isoforms'] = isoforms
                return row

    return None


def find_info_for_asm_id(chr_id, asm_id, event_dir):
    if chr_id:
        event_path = rmats_long_utils.get_chr_file_path(event_dir, chr_id)
        asm_info = find_asm_id_in_event_file(asm_id, event_path)
        if not asm_info:
            raise Exception('asm_id {} not found in {}'.format(
                asm_id, event_path))

        asm_info['chr_id'] = chr_id
        return asm_info

    names = os.listdir(event_dir)
    prefix = 'chr_id_'
    suffix = '.tsv'
    for name in names:
        has_prefix = name.startswith(prefix)
        has_suffix = name.endswith(suffix)
        if not (has_prefix and has_suffix):
            continue

        chr_id = name[len(prefix):-len(suffix)]
        path = os.path.join(event_dir, name)
        asm_info = find_asm_id_in_event_file(asm_id, path)
        if asm_info:
            asm_info['chr_id'] = chr_id
            return asm_info

    raise Exception('asm_id {} not found in {}'.format(asm_id, event_dir))


def find_gene_graph(chr_id, gene_id, event_dir):
    if chr_id:
        graph_path = rmats_long_utils.get_graph_file_path(event_dir, chr_id)
        graph_lines = rmats_long_utils.find_gene_graph_in_graph_file(
            gene_id, graph_path)
        if not graph_lines:
            raise Exception('{} not found in {}'.format(gene_id, graph_path))

        return graph_lines

    names = os.listdir(event_dir)
    prefix = 'graph_'
    suffix = '.txt'
    for name in names:
        has_prefix = name.startswith(prefix)
        has_suffix = name.endswith(suffix)
        if not (has_prefix and has_suffix):
            continue

        graph_path = os.path.join(event_dir, name)
        graph_lines = rmats_long_utils.find_gene_graph_in_graph_file(
            gene_id, graph_path)
        if graph_lines:
            return graph_lines

    raise Exception('gene_id {} not found in {}'.format(gene_id, event_dir))


def parse_attributes(parts):
    attributes = dict()
    key = None
    value_parts = list()
    for part in parts:
        if part.startswith('['):
            part = part[1:]
        if part.endswith('];'):
            part = part[:-2]

        if part == '=':
            continue

        if key is None:
            key = part
            continue

        if part.endswith(','):
            part = part[:-1]
            value_parts.append(part)
            attributes[key] = ' '.join(value_parts)
            key = None
            value_parts = list()
            continue

        value_parts.append(part)

    if value_parts:
        attributes[key] = ' '.join(value_parts)

    return attributes


def parse_graph_line(parts):
    if parts[0] == '}':
        return None

    if parts[0] == 'strict':
        return {'type': parts[1], 'name': parts[2]}

    if parts[0] in ['graph', 'edge']:
        type_str = parts[0]
        if type_str == 'edge':
            type_str = 'edge_attrs'

        attribute_parts = parts[1:]
        attributes = parse_attributes(attribute_parts)
        return {'type': type_str, 'attributes': attributes}

    if parts[:4] == ['{', 'rank', '=', 'same;']:
        raw_nodes = parts[4:-1]
        nodes = list()
        for node in raw_nodes:
            nodes.append(node.rstrip(';'))

        return {'type': 'rank', 'nodes': nodes}

    if parts[1] == '->':
        nodes = list()
        attributes = list()
        attributes_string = ''
        for part_i, part in enumerate(parts):
            if part == '->':
                continue

            if part.startswith('['):
                attribute_parts = parts[part_i:]
                attributes = parse_attributes(attribute_parts)
                attributes_string = ' '.join(attribute_parts)
                break

            nodes.append(part)

        return {
            'type': 'edge',
            'nodes': nodes,
            'attributes': attributes,
            'attributes_string': attributes_string
        }

    node = parts[0]
    attribute_parts = parts[1:]
    attributes = parse_attributes(attribute_parts)
    attributes_string = ' '.join(attribute_parts)
    return {
        'type': 'node',
        'nodes': [node],
        'attributes': attributes,
        'attributes_string': attributes_string
    }


def get_asm_node_ids(isoforms, strand):
    nodes = set()
    for isoform in isoforms:
        if strand == '-':
            isoform = reversed(isoform)

        is_start = True
        for coord in isoform:
            if is_start:
                node = 's{}'.format(coord)
                is_start = False
            else:
                node = 'e{}'.format(coord)
                is_start = True

            # The start or end node ID
            nodes.add(node)
            # The coordinate ID
            nodes.add(str(coord))

    return nodes


def collapse_nodes(pending_collapse, added_nodes, used_nodes,
                   collapse_target_by_node):
    target_coord = None
    for details in pending_collapse:
        if details['is_mid'] or (target_coord is None):
            target_coord = details['coord']

    target_node = 'combined_{}'.format(target_coord)
    mid_node = 'mid_{}'.format(target_coord)
    added_nodes.add(target_node)
    used_nodes.add(mid_node)
    for details in pending_collapse:
        coord = details['coord']
        is_mid = details['is_mid']
        if is_mid:
            templates = ['ms{}', 'me{}']
        else:
            templates = ['s{}', 'e{}']

        for template in templates:
            collapse_target_by_node[template.format(coord)] = target_node


def get_node_order_and_edges_from_graph_lines(graph_lines):
    node_order = None
    edges = list()
    for line in graph_lines:
        parts = line.split()
        parsed = parse_graph_line(parts)
        if (parsed is None) or (parsed['type'] != 'edge'):
            continue

        nodes = parsed['nodes']
        attrs = parsed['attributes']
        if attrs.get('style') == '"invis"':
            node_order = nodes
            continue

        edges.append(nodes)

    if node_order is None:
        raise Exception('Missing node order line')

    return {'order': node_order, 'edges': edges}


def determine_used_nodes_and_collapse_others(graph_lines, asm_nodes,
                                             show_extra_nodes):
    order_and_edges = get_node_order_and_edges_from_graph_lines(graph_lines)
    node_order = order_and_edges['order']
    edges = order_and_edges['edges']
    start_end_counts_by_node = dict()

    used_nodes = set()
    collapse_target_by_node = dict()
    added_nodes = set()
    pending_collapse = list()
    for coord_name in node_order:
        # the names in the node order are one of:
        #   a coordinate, mid_{coordinate}, {source,sink}_order
        if coord_name.endswith('_order'):
            continue

        name_split = coord_name.split('mid_')
        is_mid = len(name_split) == 2
        if is_mid:
            coord = name_split[1]
        else:
            coord = coord_name

        any_used = False
        start_node = 's{}'.format(coord)
        end_node = 'e{}'.format(coord)
        if (asm_nodes is None) or (start_node in asm_nodes):
            used_nodes.add(start_node)
            start_end_counts_by_node[start_node] = 0
            any_used = True
        if (asm_nodes is None) or (end_node in asm_nodes):
            used_nodes.add(end_node)
            start_end_counts_by_node[end_node] = 0
            any_used = True

        if is_mid:
            if show_extra_nodes:
                pending_collapse.append({'coord': coord, 'is_mid': is_mid})

            continue

        if any_used:
            used_nodes.add(coord)
            if pending_collapse:
                collapse_nodes(pending_collapse, added_nodes, used_nodes,
                               collapse_target_by_node)

            pending_collapse = list()
        elif show_extra_nodes:
            pending_collapse.append({'coord': coord, 'is_mid': is_mid})

    if pending_collapse:
        collapse_nodes(pending_collapse, added_nodes, used_nodes,
                       collapse_target_by_node)

    for nodes in edges:
        any_used = False
        for node in nodes:
            if node in used_nodes:
                any_used = True
                old_count = start_end_counts_by_node.get(node)
                if old_count is not None:
                    start_end_counts_by_node[node] = old_count + 1

        if not any_used:
            continue

        if 'source' in nodes:
            used_nodes.add('source')
            used_nodes.add('source_order')

        if 'sink' in nodes:
            used_nodes.add('sink')
            used_nodes.add('sink_order')

    start_or_end_to_swap = dict()
    for node in start_end_counts_by_node.keys():
        coord = node[1:]
        start_node = 's{}'.format(coord)
        end_node = 'e{}'.format(coord)
        start_count = start_end_counts_by_node.get(start_node)
        end_count = start_end_counts_by_node.get(end_node)
        if start_count is None:
            start_or_end_to_swap[start_node] = end_node
        elif end_count is None:
            start_or_end_to_swap[end_node] = start_node
        elif start_count >= end_count:
            start_or_end_to_swap[end_node] = start_node
        else:
            start_or_end_to_swap[start_node] = end_node

    return {
        'used_nodes': used_nodes,
        'collapse_target_by_node': collapse_target_by_node,
        'added_nodes': added_nodes,
        'start_or_end_to_swap': start_or_end_to_swap
    }


def collapsed_html_table_string(has_left=True, has_right=True):
    height = '20'
    width = '20'
    point_size = '40'
    if has_left:
        left = '.'
    else:
        left = ' '

    if has_right:
        right = '.'
    else:
        right = ' '

    return ('<<TABLE BORDER="0">'
            '<TR>'
            '<TD PORT="left" HEIGHT="{h}" WIDTH="{w}">'
            '<FONT POINT-SIZE="{point}">{left}</FONT><BR/> </TD>'
            '<TD PORT="right" HEIGHT="{h}" WIDTH="{w}">'
            '<FONT POINT-SIZE="{point}">{right}</FONT><BR/> </TD>'
            '</TR>'
            '</TABLE>>'.format(h=height,
                               w=width,
                               point=point_size,
                               left=left,
                               right=right))


def format_attributes(attrs):
    new_parts = list()
    for key, value in attrs.items():
        new_parts.append('{} = {}'.format(key, value))

    return '[{}];'.format(', '.join(new_parts))


def filter_node_order_to_used(parsed, used_nodes):
    new_nodes = list()
    for node in parsed['nodes']:
        if node in used_nodes:
            new_nodes.append(node)

    new_parts = list()
    for node in new_nodes:
        if new_parts:
            new_parts.append('->')

        new_parts.append(node)

    new_parts.append(parsed['attributes_string'])
    new_line = ' '.join(new_parts)
    return new_line


def filter_edge_attrs_line(line, parsed, edge_label_lines, filtered):
    if not edge_label_lines:
        filtered.append(line)
        return

    attrs = parsed['attributes']
    attrs['decorate'] = '"true"'
    new_attr_string = format_attributes(attrs)
    filtered.append('edge {}'.format(new_attr_string))


def filter_rank_line(parsed, used_nodes, start_or_end_to_swap, filtered):
    if not used_nodes.intersection(set(parsed['nodes'])):
        return

    new_rank_line = None
    new_nodes = list()
    for node in parsed['nodes']:
        if node.startswith('mid_'):
            coord = node[len('mid_'):]
            added = 'combined_{}'.format(coord)
            new_rank_line = '{{ rank = same; {}; {}}};'.format(node, added)
            break

        if node not in start_or_end_to_swap:
            new_nodes.append(node)

    if not new_rank_line:
        new_rank_line = '{{ rank = same; {}}};'.format('; '.join(new_nodes))

    filtered.append(new_rank_line)


def filter_edge_line(line, parsed, used_nodes, used_except_source_sink,
                     collapse_targets, start_or_end_to_swap, show_extra_nodes,
                     filtered, edges_by_node_pair):
    # There is an edge line with the arrow lines disabled that is
    # used to show the coordinates.
    # Filter that edge line to just the coordinates used in the ASM.
    if parsed['attributes'].get('style') == '"invis"':
        new_line = filter_node_order_to_used(parsed, used_nodes)
        filtered.append(new_line)
        return

    if len(parsed['nodes']) != 2:
        raise Exception('Expected 2 nodes for edge: {}'.format(line))

    node_a = parsed['nodes'][0]
    node_b = parsed['nodes'][1]
    swapped_a = start_or_end_to_swap.get(node_a, node_a)
    swapped_b = start_or_end_to_swap.get(node_b, node_b)
    if not {swapped_a, swapped_b}.intersection(used_except_source_sink):
        return

    a_collapsed = swapped_a not in used_nodes
    b_collapsed = swapped_b not in used_nodes
    if (not show_extra_nodes) and (a_collapsed or b_collapsed):
        return

    if a_collapsed:
        swapped_a = collapse_targets[swapped_a]
    if b_collapsed:
        swapped_b = collapse_targets[swapped_b]

    if swapped_a == swapped_b:
        return

    new_nodes = (swapped_a, swapped_b)
    for_pair = edges_by_node_pair.get(new_nodes)
    if not for_pair:
        for_pair = list()
        edges_by_node_pair[new_nodes] = for_pair

    for_pair.append({
        'parsed': parsed,
        'a_collapsed': a_collapsed,
        'b_collapsed': b_collapsed
    })


def get_edges_by_pair_and_filter_other_lines(graph_lines, used_nodes,
                                             collapse_targets,
                                             start_or_end_to_swap,
                                             edge_label_lines,
                                             show_extra_nodes, filtered):
    used_except_source_sink = used_nodes.difference({'source', 'sink'})
    edges_by_node_pair = dict()
    pending_node_details = list()
    seen_end = False
    for line in graph_lines:
        if seen_end:
            raise Exception('graph lines contined after ending "}"')

        parts = line.split()
        parsed = parse_graph_line(parts)
        if parsed is None:
            seen_end = True
            continue

        if parsed['type'] in ['digraph', 'graph']:
            filtered.append(line)
            continue

        if parsed['type'] == 'edge_attrs':
            filter_edge_attrs_line(line, parsed, edge_label_lines, filtered)
            continue

        if parsed['type'] == 'node':
            # node lines may not be written if no edges are found
            node = parsed['nodes'][0]
            if (node in used_nodes) and (node not in start_or_end_to_swap):
                node_details = {'parsed': parsed, 'node': node, 'line': line}
                pending_node_details.append(node_details)

            continue

        if parsed['type'] == 'rank':
            filter_rank_line(parsed, used_nodes, start_or_end_to_swap,
                             filtered)
            continue

        if parsed['type'] == 'edge':
            filter_edge_line(line, parsed, used_nodes, used_except_source_sink,
                             collapse_targets, start_or_end_to_swap,
                             show_extra_nodes, filtered, edges_by_node_pair)
            continue

        raise Exception('unrecognized graph line {}'.format(line))

    if not seen_end:
        raise Exception('graph lines not properly delimited')

    return {
        'edges_by_node_pair': edges_by_node_pair,
        'pending_nodes': pending_node_details
    }


def filter_and_adjust_graph(asm_info, edge_label_lines, min_edge_weight,
                            show_extra_nodes, graph_lines):
    filtered = list()
    if asm_info:
        asm_nodes = get_asm_node_ids(asm_info['isoforms'], asm_info['strand'])
    else:
        asm_nodes = None

    node_results = determine_used_nodes_and_collapse_others(
        graph_lines, asm_nodes, show_extra_nodes)
    used_nodes = node_results['used_nodes']
    added_nodes = node_results['added_nodes']
    collapse_targets = node_results['collapse_target_by_node']
    start_or_end_to_swap = node_results['start_or_end_to_swap']

    get_edges_result = get_edges_by_pair_and_filter_other_lines(
        graph_lines, used_nodes, collapse_targets, start_or_end_to_swap,
        edge_label_lines, show_extra_nodes, filtered)
    edges_by_node_pair = get_edges_result['edges_by_node_pair']
    pending_node_details = get_edges_result['pending_nodes']

    nodes_with_out_edge = set()
    nodes_with_in_edge = set()
    for pair, pair_edge_details in edges_by_node_pair.items():
        total = 0
        by_style = dict()
        for edge_details in pair_edge_details:
            parsed = edge_details['parsed']
            a_collapsed = edge_details['a_collapsed']
            b_collapsed = edge_details['b_collapsed']
            is_collapsed = a_collapsed or b_collapsed
            attrs = parsed['attributes']
            label = attrs.get('label')
            without_quotes = label.strip('"')
            if without_quotes == '':
                value = 0
            else:
                value = int(without_quotes)

            if is_collapsed and (value == 0):
                continue

            total += value
            style = attrs.get('style')
            # keep one per style to use as a template
            by_style[style] = edge_details

        if (((not by_style) or ((min_edge_weight is not None) and
                                (total < min_edge_weight)))):
            continue

        # prefer solid (None), then dashed then dotted
        style = list(by_style.keys())[0]
        if None in by_style:
            style = None
        elif '"dashed"' in by_style:
            style = '"dashed"'
        elif '"dotted"' in by_style:
            style = '"dotted"'

        edge_details = by_style[style]
        attrs = edge_details['parsed']['attributes']
        a_collapsed = edge_details['a_collapsed']
        b_collapsed = edge_details['b_collapsed']
        if a_collapsed:
            attrs['tailport'] = '"right:c"'

        if b_collapsed:
            attrs['headport'] = '"left:c"'

        if 'label' in attrs:
            if edge_label_lines:
                # leading spaces improve the line connecting the label
                attrs['label'] = '"  {}"'.format(total)
            else:
                attrs['label'] = '"{}"'.format(total)

        new_attr_string = format_attributes(attrs)
        new_edge_line = '{} -> {} {}'.format(pair[0], pair[1], new_attr_string)
        nodes_with_out_edge.add(pair[0])
        nodes_with_in_edge.add(pair[1])
        filtered.append(new_edge_line)

    for node_details in pending_node_details:
        node = node_details['node']
        node_line = node_details['line']
        parsed = node_details['parsed']
        attrs = parsed['attributes']
        is_coord_node = attrs['label'].isdigit()
        # Keep nodes with edges and the nodes showing the coordinate order
        if (((node in nodes_with_in_edge) or (node in nodes_with_out_edge)
             or is_coord_node)):
            filtered.append(node_line)
            continue

        # Hide other nodes
        attrs['label'] = '""'
        attrs['shape'] = '"none"'
        new_attr_string = format_attributes(attrs)
        new_node_line = '{} {}'.format(node, new_attr_string)
        filtered.append(new_node_line)
        continue

    for node in added_nodes:
        has_left = node in nodes_with_in_edge
        has_right = node in nodes_with_out_edge
        table_label = collapsed_html_table_string(has_left=has_left,
                                                  has_right=has_right)
        node_line = '{} [label = {}, shape = "none"];'.format(
            node, table_label)
        filtered.append(node_line)

    filtered.append('}')
    return filtered


def plot_graph(graph_lines, out_file):
    graph_string = '\n'.join(graph_lines)
    graphs = pydot.graph_from_dot_data(graph_string)
    graph = graphs[0]
    _, extension = os.path.splitext(out_file)
    extension = extension.lower()
    if extension == '.png':
        graph.write_png(out_file)
    elif extension == '.pdf':
        graph.write_pdf(out_file)
    elif extension == '.svg':
        graph.write_svg(out_file)
    else:
        print('Unrecognized format for {}. Writing .svg content'.format(
            out_file))
        graph.write_svg(out_file)


def plot_splice_graph(chrom, gene_id, asm_id, event_dir, edge_label_lines,
                      min_edge_weight, show_extra_nodes, out_file):
    asm_info = None
    chr_id_mapping = rmats_long_utils.parse_chr_mapping(event_dir)
    chr_to_id = chr_id_mapping['chr_to_id']
    chr_id = None
    if chrom:
        chr_id = chr_to_id.get(chrom)
        if chr_id is None:
            raise Exception(
                'chr_id not found for {} in mapping from {}'.format(
                    chrom, event_dir))

    if asm_id:
        asm_info = find_info_for_asm_id(chr_id, asm_id, event_dir)
        chr_id = asm_info['chr_id']
        if gene_id:
            if asm_info['gene_id'] != gene_id:
                raise Exception(
                    'asm did not match provided gene_id: {} {}'.format(
                        gene_id, asm_info))
        else:
            gene_id = asm_info['gene_id']

    graph_lines = find_gene_graph(chr_id, gene_id, event_dir)
    graph_lines = filter_and_adjust_graph(asm_info, edge_label_lines,
                                          min_edge_weight, show_extra_nodes,
                                          graph_lines)
    plot_graph(graph_lines, out_file)


def main():
    args = parse_args()
    plot_splice_graph(args.chr, args.gene_id, args.asm_id, args.event_dir,
                      args.edge_label_lines, args.min_edge_weight,
                      args.show_extra_nodes, args.out_file)
    print('plot_splice_graph.py finished')


if __name__ == '__main__':
    main()
