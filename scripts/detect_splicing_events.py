import argparse
import bisect
import multiprocessing
import os
import os.path

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Detect alternative splicing events from'
                     ' a set of isoforms'))
    parser.add_argument(
        '--gtf-dir',
        required=True,
        help='The output directory from organize_gene_info_by_chr.py')
    parser.add_argument('--align-dir',
                        required=False,
                        help=('The output directory from'
                              ' organize_alignment_info_by_gene_and_chr.py'))
    parser.add_argument(
        '--out-dir',
        required=True,
        help='The directory to create and where chr files will be written')
    parser.add_argument(
        '--max-nodes-in-event',
        default=50,
        type=int,
        help=('Only look for events between nodes (splice sites) in the'
              ' splice graph that are at most --max-nodes-in-event apart.'
              ' Default: %(default)s'))
    parser.add_argument(
        '--max-paths-in-event',
        default=50,
        type=int,
        help=('Only output events with at most --max-paths-in-event.'
              ' Default: %(default)s'))
    parser.add_argument('--num-threads',
                        default=1,
                        type=int,
                        help='how many threads to use. Default: %(default)s')
    parser.add_argument(
        '--min-reads-per-edge',
        default=5,
        type=int,
        help=('Only include an edge in the splice graph if there are at least'
              ' this many supporting reads. Default: %(default)s'))
    parser.add_argument('--output-full-gene-asm',
                        action='store_true',
                        help=('Output only ASMs that cover the entire gene'
                              ' (from transcript start to end)'))
    parser.add_argument('--simplify-gene-isoform-endpoints',
                        action='store_true',
                        help=('Combine gene isoforms where the only difference'
                              ' is the transcripts start and/or end'))
    parser.add_argument(
        '--filter-gene-isoforms-by-edge',
        action='store_true',
        help=('With --output-full-gene-asm, require each isoform to have'
              ' --min-reads-per-edge for each splice junction'))

    args = parser.parse_args()
    # Without the align_dir the edges cannot be filtered by read count
    if not args.align_dir:
        args.min_reads_per_edge = 0

    return args


def parse_gtf_line(line):
    parsed = dict()
    columns = rmats_long_utils.read_tsv_line(line)
    num_columns = len(columns)
    if num_columns == 4:
        gene_id = columns[2]
        parsed['gene'] = gene_id
        return parsed

    start_str = columns[0]
    transcript_start = int(start_str)
    end_str = columns[1]
    transcript_end = int(end_str)
    sjs_string = columns[2]
    sjs = rmats_long_utils.parse_sjs_string(sjs_string)
    strand = columns[3]
    transcript_id = columns[4]
    parsed['transcript_id'] = transcript_id
    parsed['strand'] = strand
    exons = list()
    exons.append([transcript_start, None])
    for sj_start, sj_end in sjs:
        exons[-1][1] = sj_start
        exons.append([sj_end, None])

    exons[-1][1] = transcript_end
    parsed['exons'] = exons
    return parsed


class Node:
    '''
    A Node in the splice graph is one of:
    * supersource (^): All transcripts start here
    * supersink ($): All transcripts end here
    * exon start: The start coordinate of an exon
    * exon end: The end coordinate of an exon

    Example:
    * Transcripts
      + (^)->(a,b)->(c,d)->(e,f)-> (g,h)->(i,j)->($)
      + (^)->(a,b)->       (e,f)-> (g,h)->(i,j)->($)
      + (^)->(a,b)->(c,d)->(e,f)->(g*,h)->(i,j)->($)
      + (g*,h) is like exon (g,h) but with an earlier start coordinate
      + a, c, e, g, i are exon start nodes
      + b, d, f, h, j are exon end nodes
    * Final splice graph:
      ^ -> a -> b ------> e -> f -> g*
                 \       /      \    \
                  c -> d         g -> h -> i -> j -> $

    The implementations of hash and equality do not include self.edges.
    This allows keeping the edges within the Node as a dictionary.
    This node has a directed edge to node_b if node_b in self.edges
    The value in the dictionary is the edge weight which defaults to
    zero and can be assigned later.
    '''
    def __init__(self, coordinate, node_type):
        self.edges = dict()
        self.from_nodes = set()
        self.coordinate = coordinate
        self.node_type = node_type
        self.removed = dict()
        self.num_edges_from_source = 0
        self.num_edges_from_non_source = 0

    def _compare_key(self):
        return (self.coordinate, self.node_type)

    def __hash__(self):
        return hash(self._compare_key())

    def __eq__(self, other):
        return self._compare_key() == other._compare_key()

    def __lt__(self, other):
        if self == other:
            return False

        if (self.node_type == 'source') or (other.node_type == 'sink'):
            return True

        if (self.node_type == 'sink') or (other.node_type == 'source'):
            return False

        if self.coordinate is None and other.coordinate is None:
            return self.node_type < other.node_type

        if self.coordinate is None:
            return True

        if other.coordinate is None:
            return False

        return self._compare_key() < other._compare_key()

    def __str__(self):
        return str(self._compare_key())

    def __repr__(self):
        return repr({
            'edges': [str(x) for x in self.edges],
            'coordinate': self.coordinate,
            'node_type': self.node_type
        })


def reverse_minus_strand_exons(exons):
    reversed_exons = list()
    for exon in reversed(exons):
        end, start = exon
        reversed_exons.append((start, end))

    return reversed_exons


def add_node_edge(from_node, to_node):
    from_node.edges[to_node] = 0
    if from_node.node_type == 'source':
        to_node.num_edges_from_source += 1
    else:
        to_node.num_edges_from_non_source += 1

    to_node.from_nodes.add(from_node)


def remove_node_edge(from_node, to_node):
    del from_node.edges[to_node]
    if from_node.node_type == 'source':
        to_node.num_edges_from_source -= 1
    else:
        to_node.num_edges_from_non_source -= 1

    to_node.from_nodes.remove(from_node)


def update_if_closer(container, key_node, new_node):
    old_node = container.get(key_node, None)
    if old_node is None or old_node.coordinate is None:
        container[key_node] = new_node
        return

    if new_node.coordinate is None or key_node.coordinate is None:
        return

    old_diff = abs(old_node.coordinate - key_node.coordinate)
    new_diff = abs(new_node.coordinate - key_node.coordinate)
    if new_diff < old_diff:
        container[key_node] = new_node


# build_splice_graph returns:
#  splice_graph: the source node of the graph
#  sink: the sink node of the graph
#  closest_downstream_by_node: stores one outgoing edge for each node
#  closest_upstream_by_node: stores one incoming edge for each node
# The extra closest_downstream_by_node and closest_upstream_by_node are used
# to provide missing exon start/end coordinates for splicing events
# that don't depend on the full exon
def build_splice_graph(transcripts, strand):
    all_nodes = dict()
    closest_downstream_by_node = dict()
    closest_upstream_by_node = dict()
    splice_graph = Node(None, 'source')
    all_nodes[splice_graph] = splice_graph
    sink_node = Node(None, 'sink')
    all_nodes[sink_node] = sink_node
    for transcript_details in transcripts:
        exons = transcript_details['exons']
        if strand == '-':
            exons = reverse_minus_strand_exons(exons)

        node = splice_graph
        for exon in exons:
            start, end = exon
            start_node = Node(start, 'start')
            existing_start = all_nodes.get(start_node)
            if existing_start:
                start_node = existing_start
            else:
                all_nodes[start_node] = start_node

            add_node_edge(node, start_node)
            update_if_closer(closest_downstream_by_node, node, start_node)
            update_if_closer(closest_upstream_by_node, start_node, node)
            # move to exon start node
            node = start_node

            end_node = Node(end, 'end')
            existing_end = all_nodes.get(end_node)
            if existing_end:
                end_node = existing_end
            else:
                all_nodes[end_node] = end_node

            add_node_edge(node, end_node)
            update_if_closer(closest_downstream_by_node, node, end_node)
            update_if_closer(closest_upstream_by_node, end_node, node)
            # move to exon end node
            node = end_node

        # add sink node at end of every transcript
        add_node_edge(node, sink_node)
        update_if_closer(closest_downstream_by_node, node, sink_node)
        update_if_closer(closest_upstream_by_node, sink_node, node)

    results = dict()
    results['graph'] = splice_graph
    results['sink'] = sink_node
    results['closest_down'] = closest_downstream_by_node
    results['closest_up'] = closest_upstream_by_node
    return results


def filter_running_paths(paths, cutoff_node, strand):
    if cutoff_node is None or cutoff_node.coordinate is None:
        return paths

    cutoff = cutoff_node.coordinate
    new_paths_set = set()
    new_paths = list()
    for path in paths:
        for node_i, node in enumerate(path):
            if (node.coordinate is not None
                    and ((strand == '+') and (node.coordinate >= cutoff) or
                         (strand == '-') and (node.coordinate <= cutoff))):
                new_path = path[node_i:]
                new_path_tuple = tuple(new_path)
                if new_path_tuple not in new_paths_set:
                    new_paths_set.add(new_path_tuple)
                    new_paths.append(new_path)

                break

    return new_paths


def detect_bubbles_visit_node(node, nodes_to_visit, running_paths_by_node,
                              cutoff_node, strand, max_paths_in_event,
                              bubbles):
    if node.node_type == 'source':
        running_path = [node]
        running_paths_by_node[node] = [running_path]

    paths_by_start = dict()
    # Stop tracking paths that end at this node
    running_paths = running_paths_by_node.pop(node)
    running_paths = filter_running_paths(running_paths, cutoff_node, strand)
    for running_path in running_paths:
        # Add all start-to-node paths that could be part of a bubble.
        # If a start node does not have multiple outgoing edges then it
        # cannot be the start of a bubble.
        for start_node_i in range(len(running_path) - 1):
            start_node = running_path[start_node_i]
            if len(start_node.edges) < 2:
                continue

            paths_for_start_end = paths_by_start.get(start_node)
            if not paths_for_start_end:
                paths_for_start_end = set()
                paths_by_start[start_node] = paths_for_start_end

            paths_for_start_end.add(tuple(running_path[start_node_i:]))

    detect_bubbles_from_paths(paths_by_start, max_paths_in_event, bubbles)
    if len(running_paths_by_node) == 0:
        # All currently tracked paths end at this node.
        # Reset to a single path with this node.
        running_path = [node]
        running_paths = [running_path]

    for other in node.edges:
        nodes_to_visit.add(other)
        paths_to_other = running_paths_by_node.get(other)
        if not paths_to_other:
            paths_to_other = list()
            running_paths_by_node[other] = paths_to_other

        for running_path in running_paths:
            extended = running_path + [other]
            paths_to_other.append(extended)


def choose_next_node(nodes_to_visit, strand):
    next_node = None
    for node in nodes_to_visit:
        if next_node is None:
            next_node = node
            continue

        if node.coordinate is None:
            continue

        if next_node.coordinate is None:
            next_node = node
            continue

        # It's possible for a coordinate to be both a start and end of an exon
        if (((node.coordinate == next_node.coordinate)
             and (node.node_type == 'start'))):
            next_node = node
            continue

        if (strand == '+') and (node.coordinate < next_node.coordinate):
            next_node = node
            continue

        if (strand == '-') and (node.coordinate > next_node.coordinate):
            next_node = node
            continue

    nodes_to_visit.remove(next_node)
    return next_node


def detect_bubbles_from_paths(by_start, max_paths_in_event, bubbles):
    for paths in by_start.values():
        num_paths = len(paths)
        if (num_paths < 2) or (num_paths > max_paths_in_event):
            continue

        shared_nodes = None
        for path in paths:
            if shared_nodes is None:
                shared_nodes = set(path)
            else:
                shared_nodes.intersection_update(set(path))

        # A bubble is multiple paths between two nodes.
        # If all the paths share a node other than the
        # start and end nodes then there is a smaller bubble to find
        if len(shared_nodes) == 2:
            bubble = {'paths': paths, 'path_to_id': dict()}
            bubbles.append(bubble)


def detect_bubbles(splice_graph, strand, max_nodes_in_event,
                   max_paths_in_event):
    # Start at the source node and traverse the graph.
    # When choosing the next node prefer nodes with lower coordinates, but
    # if the graph is for the minus strand then prefer higher coordinates.
    # That way, when a node is visited all nodes that link to it will already
    # have been visited.
    # This assumes no loops and that strand information is correct.
    # While doing the traversal keep track of "running paths".
    # If at any time all running paths share a node then stop tracking
    # those paths. This is because any bubble would end at the shared node.
    # To limit the search and ensure reasonable running time and memory usage,
    # store the order of nodes_visited. If a running path contains a node that
    # is more than max_nodes_in_event back in nodes_visited then remove that
    # node from the running path and simplify the remaining running paths.
    nodes_to_visit = {splice_graph}
    running_paths_by_node = dict()
    nodes_visited = list()
    bubbles = list()
    while nodes_to_visit:
        next_node = choose_next_node(nodes_to_visit, strand)
        nodes_visited.append(next_node)
        cutoff_node = None
        if len(nodes_visited) > max_nodes_in_event:
            cutoff_node = nodes_visited[-max_nodes_in_event]

        detect_bubbles_visit_node(next_node, nodes_to_visit,
                                  running_paths_by_node, cutoff_node, strand,
                                  max_paths_in_event, bubbles)

    return bubbles


def get_isoform_exons(path, closest_downstream_by_node,
                      closest_upstream_by_node, output_full_gene_asm):
    exons = list()
    for node in path:
        if node.node_type in ['source', 'sink']:
            continue

        if node.node_type == 'start':
            if output_full_gene_asm:
                exons.append([node.coordinate, None])
            else:
                default_end = closest_downstream_by_node[node]
                exons.append([node.coordinate, default_end.coordinate])
        elif node.node_type == 'end':
            if exons:
                exons[-1][-1] = node.coordinate
            else:
                default_start = closest_upstream_by_node[node]
                exons.append([default_start.coordinate, node.coordinate])

    return [tuple(exon) for exon in exons]


def is_exon_skipping(paths):
    if len(paths) != 2:
        return False

    include_path = None
    skip_path = None
    for path in paths:
        if len(path) == 2:
            skip_path = path
        if len(path) == 4:
            include_path = path

    if include_path is None or skip_path is None:
        return False

    return ((include_path[0] == skip_path[0])
            and (include_path[3] == skip_path[1])
            and (include_path[0].node_type == 'end')
            and (include_path[1].node_type == 'start')
            and (include_path[2].node_type == 'end')
            and (include_path[3].node_type == 'start'))


def is_alt_3_splice_site(paths):
    if len(paths) != 2:
        return False

    path_list = list()
    for path in paths:
        if len(path) != 3:
            return False

        path_list.append(path)

    path_a = path_list[0]
    path_b = path_list[1]
    return ((path_a[0] == path_b[0]) and (path_a[2] == path_b[2])
            and (path_a[1].node_type == path_b[1].node_type)
            and (path_a[1].coordinate != path_b[1].coordinate)
            and (path_a[0].node_type == 'end')
            and (path_a[1].node_type == 'start')
            and (path_a[2].node_type == 'end'))


def is_alt_5_splice_site(paths):
    if len(paths) != 2:
        return False

    path_list = list()
    for path in paths:
        if len(path) != 3:
            return False

        path_list.append(path)

    path_a = path_list[0]
    path_b = path_list[1]
    return ((path_a[0] == path_b[0]) and (path_a[2] == path_b[2])
            and (path_a[1].node_type == path_b[1].node_type)
            and (path_a[1].coordinate != path_b[1].coordinate)
            and (path_a[0].node_type == 'start')
            and (path_a[1].node_type == 'end')
            and (path_a[2].node_type == 'start'))


def is_alt_first_exon(paths):
    if len(paths) != 2:
        return False

    path_list = list()
    for path in paths:
        if len(path) != 4:
            return False

        path_list.append(path)

    path_a = path_list[0]
    path_b = path_list[1]
    return ((path_a[0] == path_b[0]) and (path_a[3] == path_b[3])
            and (path_a[1].node_type == path_b[1].node_type)
            and (path_a[1].coordinate != path_b[1].coordinate)
            and (path_a[2].node_type == path_b[2].node_type)
            and (path_a[2].coordinate != path_b[2].coordinate)
            and (path_a[0].node_type == 'source')
            and (path_a[1].node_type == 'start')
            and (path_a[2].node_type == 'end')
            and (path_a[3].node_type == 'start'))


def is_alt_last_exon(paths):
    if len(paths) != 2:
        return False

    path_list = list()
    for path in paths:
        if len(path) != 4:
            return False

        path_list.append(path)

    path_a = path_list[0]
    path_b = path_list[1]
    return ((path_a[0] == path_b[0]) and (path_a[3] == path_b[3])
            and (path_a[1].node_type == path_b[1].node_type)
            and (path_a[1].coordinate != path_b[1].coordinate)
            and (path_a[2].node_type == path_b[2].node_type)
            and (path_a[2].coordinate != path_b[2].coordinate)
            and (path_a[0].node_type == 'end')
            and (path_a[1].node_type == 'start')
            and (path_a[2].node_type == 'end')
            and (path_a[3].node_type == 'sink'))


def is_intron_retention(paths):
    if len(paths) != 2:
        return False

    retained_path = None
    spliced_path = None
    for path in paths:
        if len(path) == 4:
            spliced_path = path
        elif len(path) == 2:
            retained_path = path
        else:
            return False

    return ((retained_path is not None) and (spliced_path is not None)
            and (spliced_path[0] == retained_path[0])
            and (spliced_path[3] == retained_path[1])
            and (spliced_path[0].node_type == 'start')
            and (spliced_path[1].node_type == 'end')
            and (spliced_path[2].node_type == 'start')
            and (spliced_path[3].node_type == 'end'))


def is_mutually_exclusive_exons(paths):
    if len(paths) != 2:
        return False

    path_list = list()
    for path in paths:
        if len(path) != 4:
            return False

        path_list.append(path)

    path_a = path_list[0]
    path_b = path_list[1]
    return ((path_a[0] == path_b[0]) and (path_a[3] == path_b[3])
            and (path_a[1].node_type == path_b[1].node_type)
            and (path_a[1].coordinate != path_b[1].coordinate)
            and (path_a[2].node_type == path_b[2].node_type)
            and (path_a[2].coordinate != path_b[2].coordinate)
            and (path_a[0].node_type == 'end')
            and (path_a[1].node_type == 'start')
            and (path_a[2].node_type == 'end')
            and (path_a[3].node_type == 'start'))


def group_path_info_by_internal_nodes(paths, path_to_id):
    paths_by_internal = dict()
    result = {'any_alt': False, 'paths_by_internal': paths_by_internal}
    for path_i, path in enumerate(paths):
        # All paths share the first and last node. Just check on 1st path
        if path_i == 0:
            has_start = path[0].node_type == 'source'
            has_end = path[-1].node_type == 'sink'
            # Only simplify if the paths include the source or sink node.
            if not (has_start or has_end):
                return result

            internal_start_offset = 1
            internal_end_offset = 1
            if has_start:
                internal_start_offset += 1
            if has_end:
                internal_end_offset += 1

        # The ends of the paths are always the same for paths in an event.
        # Also transcript start and end sites are not "internal"
        internal_list = list()
        for i in range(internal_start_offset, len(path) - internal_end_offset):
            internal_list.append(path[i].coordinate)

        internal = tuple(internal_list)
        details = paths_by_internal.get(internal)
        if not details:
            details = {
                'paths': list(),
                'starts': set(),
                'ends': set(),
                'starts_by_coord': dict(),
                'ends_by_coord': dict(),
                'ids_by_start': dict(),
                'ids_by_end': dict()
            }
            paths_by_internal[internal] = details
        else:
            result['any_alt'] = True

        # convert from a tuple to a list to allow modifying
        details['paths'].append(list(path))
        path_id = path_to_id.get(path)
        if has_start:
            start_coord = path[1].coordinate
            details['starts'].add(start_coord)
            details['starts_by_coord'][start_coord] = path[1]
            ids_for_start = rmats_long_utils.try_get_or_set_default(
                details['ids_by_start'], start_coord, set())
            ids_for_start.add(path_id)

        if has_end:
            end_coord = path[-2].coordinate
            details['ends'].add(end_coord)
            details['ends_by_coord'][end_coord] = path[-2]
            ids_for_end = rmats_long_utils.try_get_or_set_default(
                details['ids_by_end'], end_coord, set())
            ids_for_end.add(path_id)

    return result


def check_for_shared_starts_and_ends(paths_by_internal):
    shared_starts = None
    shared_ends = None
    for details in paths_by_internal.values():
        if shared_starts is None:
            shared_starts = details['starts']
        else:
            shared_starts = shared_starts.intersection(details['starts'])

        if shared_ends is None:
            shared_ends = details['ends']
        else:
            shared_ends = shared_ends.intersection(details['ends'])

    return {'starts': shared_starts, 'ends': shared_ends}


def get_adjusted_start_for_strand(strand, details):
    if strand == '+':
        return min(details['starts'])

    return max(details['starts'])


def get_adjusted_end_for_strand(strand, details):
    if strand == '+':
        return max(details['ends'])

    return min(details['ends'])


def update_path_id_from_isoforms_after_simplification(start_ids, end_ids,
                                                      tuple_path, path_to_id):
    if start_ids and end_ids:
        intersection = start_ids.intersection(end_ids)
        if intersection:
            first_id = sorted(intersection)[0]
            path_to_id[tuple_path] = first_id
        else:
            first_start_id = sorted(start_ids)[0]
            first_end_id = sorted(end_ids)[0]
            sorted_ids = sorted([first_start_id, first_end_id])
            combined_id = '_'.join(sorted_ids)
            path_to_id[tuple_path] = combined_id
    elif start_ids:
        first_id = sorted(start_ids)[0]
        path_to_id[tuple_path] = first_id
    elif end_ids:
        first_id = sorted(end_ids)[0]
        path_to_id[tuple_path] = first_id


# If there is a pair of paths where the only difference is the
# transcript start or the transcript end then remove one of the paths
def simplify_event_paths(paths, strand, output_full_gene_asm,
                         simplify_gene_isoform_endpoints, path_to_id):
    if output_full_gene_asm and (not simplify_gene_isoform_endpoints):
        return paths

    group_result = group_path_info_by_internal_nodes(paths, path_to_id)
    any_alt = group_result['any_alt']
    paths_by_internal = group_result['paths_by_internal']
    if not any_alt:
        return paths

    shared_starts_and_ends = check_for_shared_starts_and_ends(
        paths_by_internal)
    shared_starts = shared_starts_and_ends['starts']
    shared_ends = shared_starts_and_ends['ends']
    # If there is a transcript start or end that is available for all internal
    # paths then there should be a smaller event using that shared node.
    if (shared_starts or shared_ends) and (not output_full_gene_asm):
        return list()

    new_paths = list()
    for details in paths_by_internal.values():
        if len(details['paths']) == 1:
            new_paths.append(tuple(details['paths'][0]))
            continue

        start_ids = set()
        end_ids = set()
        adjusted_path = details['paths'][0][:]
        if details['starts']:
            adjusted_start = get_adjusted_start_for_strand(strand, details)
            adjusted_path[1] = details['starts_by_coord'][adjusted_start]
            start_ids = details['ids_by_start'][adjusted_start]

        if details['ends']:
            adjusted_end = get_adjusted_end_for_strand(strand, details)
            adjusted_path[-2] = details['ends_by_coord'][adjusted_end]
            end_ids = details['ids_by_end'][adjusted_end]

        tuple_path = tuple(adjusted_path)
        new_paths.append(tuple_path)
        if simplify_gene_isoform_endpoints:
            update_path_id_from_isoforms_after_simplification(
                start_ids, end_ids, tuple_path, path_to_id)

    return new_paths


def get_start_and_end_from_simplified_paths(paths):
    first_path = list(paths)[0]
    start_node = first_path[0]
    end_node = first_path[-1]
    if start_node.node_type == 'source':
        start = 'source'
    else:
        start = start_node.coordinate

    if end_node.node_type == 'sink':
        end = 'sink'
    else:
        end = end_node.coordinate

    return {
        'start_node': start_node,
        'end_node': end_node,
        'start': start,
        'end': end
    }


def get_splice_site_always_and_never(start_node, end_node):
    end_always_ss = True
    end_never_ss = True
    for after_end_node in end_node.edges:
        if after_end_node.node_type == 'sink':
            end_always_ss = False
        else:
            end_never_ss = False

    if end_node.node_type == 'sink':
        end_always_ss = False
        end_never_ss = False

    start_always_ss = True
    start_never_ss = True
    if start_node.num_edges_from_source > 0:
        start_always_ss = False
    if start_node.num_edges_from_non_source > 0:
        start_never_ss = False
    if start_node.node_type == 'source':
        start_always_ss = False
        start_never_ss = False

    return {
        'start_always_ss': start_always_ss,
        'start_never_ss': start_never_ss,
        'end_always_ss': end_always_ss,
        'end_never_ss': end_never_ss
    }


def edge_starts_at_or_after_end_node(edge_start, end_node, strand):
    if end_node.node_type == 'sink':
        return False

    if edge_start.node_type == 'source':
        return False

    if strand == '+':
        return edge_start.coordinate >= end_node.coordinate

    return edge_start.coordinate <= end_node.coordinate


def edge_ends_at_or_before_start_node(edge_end, start_node, strand):
    if start_node.node_type == 'source':
        return False

    if edge_end.node_type == 'sink':
        return False

    if strand == '+':
        return edge_end.coordinate <= start_node.coordinate

    return edge_end.coordinate >= start_node.coordinate


# An ASM meets the strict definition if there are no edges that
# overlap the ASM which are not part of the ASM
def check_if_strict_asm(paths, all_edges, strand, start_node, end_node):
    asm_edges = set()
    prev_node = None
    for path in paths:
        seen_start = False
        for node in path:
            if not seen_start:
                if node == start_node:
                    seen_start = True
                    prev_node = node

                continue

            edge = (prev_node, node)
            asm_edges.add(edge)
            prev_node = node
            if node == end_node:
                break

    non_asm_edges = all_edges.difference(asm_edges)
    for edge_start, edge_end in non_asm_edges:
        if ((edge_starts_at_or_after_end_node(edge_start, end_node, strand) or
             edge_ends_at_or_before_start_node(edge_end, start_node, strand))):
            continue

        return False

    return True


def classify_bubble(bubble, closest_downstream_by_node,
                    closest_upstream_by_node, all_edges, strand,
                    output_full_gene_asm, simplify_gene_isoform_endpoints):
    paths = bubble['paths']
    path_to_id = bubble['path_to_id']
    isoforms_and_ids = list()
    simplified_paths = simplify_event_paths(paths, strand,
                                            output_full_gene_asm,
                                            simplify_gene_isoform_endpoints,
                                            path_to_id)
    if len(simplified_paths) < 2:
        return None

    for path in simplified_paths:
        isoform_id = path_to_id.get(path)
        isoform_exons = get_isoform_exons(path, closest_downstream_by_node,
                                          closest_upstream_by_node,
                                          output_full_gene_asm)
        if strand == '-':
            isoform_exons = reverse_minus_strand_exons(isoform_exons)

        isoforms_and_ids.append((isoform_exons, isoform_id))

    start_end_details = get_start_and_end_from_simplified_paths(
        simplified_paths)
    start_node = start_end_details['start_node']
    end_node = start_end_details['end_node']
    start = start_end_details['start']
    end = start_end_details['end']
    ss_details = get_splice_site_always_and_never(start_node, end_node)
    is_strict = check_if_strict_asm(paths, all_edges, strand, start_node,
                                    end_node)
    isoforms_and_ids.sort()
    isoforms = list()
    isoform_ids = list()
    for isoform, isoform_id in isoforms_and_ids:
        isoforms.append(isoform)
        isoform_ids.append(isoform_id)

    event = {
        'strand': strand,
        'isoforms': isoforms,
        'isoform_ids': isoform_ids,
        'start': start,
        'end': end,
        'start_always_ss': ss_details['start_always_ss'],
        'start_never_ss': ss_details['start_never_ss'],
        'end_always_ss': ss_details['end_always_ss'],
        'end_never_ss': ss_details['end_never_ss'],
        'is_strict': is_strict
    }
    if is_exon_skipping(simplified_paths):
        event['type'] = 'exon skipping'
    elif is_alt_3_splice_site(simplified_paths):
        event['type'] = "alternative 3'-splice site usage"
    elif is_alt_5_splice_site(simplified_paths):
        event['type'] = "alternative 5'-splice site usage"
    elif is_alt_first_exon(simplified_paths):
        event['type'] = 'alternative first exon usage'
    elif is_alt_last_exon(simplified_paths):
        event['type'] = 'alternative last exon usage'
    elif is_intron_retention(simplified_paths):
        event['type'] = 'intron retention'
    elif is_mutually_exclusive_exons(simplified_paths):
        event['type'] = 'mutually exclusive exons'
    else:
        event['type'] = 'complex splicing'

    return event


def get_edge_set_from_graph(splice_graph):
    edges = set()
    to_visit = set()
    visited = set()
    to_visit.add(splice_graph)
    while to_visit:
        node = to_visit.pop()
        visited.add(node)
        for other in node.edges:
            edge = (node, other)
            edges.add(edge)
            if other not in visited:
                to_visit.add(other)

    return edges


def classify_bubbles(bubbles, closest_downstream_by_node,
                     closest_upstream_by_node, splice_graph, strand,
                     output_full_gene_asm, simplify_gene_isoform_endpoints):
    events = list()
    all_edges = get_edge_set_from_graph(splice_graph)
    for bubble in bubbles:
        event = classify_bubble(bubble, closest_downstream_by_node,
                                closest_upstream_by_node, all_edges, strand,
                                output_full_gene_asm,
                                simplify_gene_isoform_endpoints)
        if event:
            events.append(event)

    return events


def create_bubble_from_transcripts(transcripts, strand):
    if not transcripts:
        return list()

    bubble_paths = set()
    path_to_id = dict()
    for transcript in transcripts:
        path = list()
        path.append(Node(None, 'source'))
        exons = transcript['exons']
        if strand == '-':
            exons = reverse_minus_strand_exons(exons)

        for start, end in exons:
            path.append(Node(start, 'start'))
            path.append(Node(end, 'end'))

        path.append(Node(None, 'sink'))
        path_tuple = tuple(path)
        # Transcripts could have the same exon coords:
        # (ENST00000602634.5, ENST00000644515.1)
        if path_tuple in bubble_paths:
            continue

        bubble_paths.add(path_tuple)
        path_to_id[path_tuple] = transcript['transcript_id']

    bubble = {'paths': bubble_paths, 'path_to_id': path_to_id}
    return [bubble]


def separate_plus_and_minus_strand_transcripts(transcripts):
    plus_transcripts = list()
    minus_transcripts = list()
    for transcript in transcripts:
        strand = transcript['strand']
        if strand == '+':
            plus_transcripts.append(transcript)
        elif strand == '-':
            minus_transcripts.append(transcript)
        else:
            print('skipping {} with unknown strand {}'.format(
                transcript['transcript_id'], strand))

    return {'plus': plus_transcripts, 'minus': minus_transcripts}


def reverse_sjs(sjs):
    new_sjs = list()
    for sj_start, sj_end in sjs:
        new_sjs.append((sj_end, sj_start))

    return reversed(new_sjs)


def lookup_node_details(node_details, node_type, coordinate=None):
    if node_type == 'source':
        return node_details['source']

    if node_type == 'sink':
        return node_details['sink']

    orig_nodes = node_details['orig_nodes']
    orig_node_key = Node(coordinate, node_type)
    node = orig_nodes.get(orig_node_key)
    if node:
        return node

    if node_type == 'start':
        mid_type = 'mid_start'
    else:
        mid_type = 'mid_end'

    extra_nodes = node_details['extra_nodes']
    # For the minus strand, negated_coords are in ascending order and
    # can be used with bisect.
    coords = node_details['coords']
    neg_coords = node_details['negated_coords']
    if neg_coords is not None:
        index = bisect.bisect_left(neg_coords, -coordinate)
    else:
        index = bisect.bisect_left(coords, coordinate)

    # index is the offset for the first coord that is >= coordinate.
    # If no such value then index is len(coords).
    if index == len(coords):
        mid_node_key = Node(coords[-1], mid_type)
        return extra_nodes[mid_node_key]

    matched_coord = coords[index]
    if matched_coord == coordinate:
        extra_key = Node(coordinate, node_type)
        return extra_nodes[extra_key]

    if index == 0:
        after_source_key = Node(None, mid_type)
        return extra_nodes[after_source_key]

    lower_coord = coords[index - 1]
    mid_key = Node(lower_coord, mid_type)
    return extra_nodes[mid_key]


def increment_edge_weight(node_details, from_node, to_node):
    from_is_extra = from_node not in node_details['orig_nodes']
    old_count = from_node.edges.get(to_node)
    if old_count is not None:
        from_node.edges[to_node] = old_count + 1
        return

    if from_is_extra:
        from_node.edges[to_node] = 1
        return

    # Separately track the weights from original nodes for edges that
    # are not in the original splice graph.
    extra_edges = node_details['extra_edges']
    extra_for_from = rmats_long_utils.try_get_or_set_default(
        extra_edges, from_node, dict())
    old_count = extra_for_from.get(to_node, 0)
    extra_for_from[to_node] = old_count + 1


def add_read_count_to_matching_edges(start, end, sjs, node_details):
    if not node_details['coords']:
        # No transcripts for this strand
        return

    source = lookup_node_details(node_details, 'source')
    sink = lookup_node_details(node_details, 'sink')
    read_start_node = lookup_node_details(node_details,
                                          'start',
                                          coordinate=start)
    read_end_node = lookup_node_details(node_details, 'end', coordinate=end)
    increment_edge_weight(node_details, source, read_start_node)
    increment_edge_weight(node_details, read_end_node, sink)

    prev_exon_start_node = read_start_node
    for sj in sjs:
        sj_start, sj_end = sj
        # end of exon
        end_node = lookup_node_details(node_details,
                                       'end',
                                       coordinate=sj_start)
        # start of exon
        start_node = lookup_node_details(node_details,
                                         'start',
                                         coordinate=sj_end)
        increment_edge_weight(node_details, prev_exon_start_node, end_node)
        increment_edge_weight(node_details, end_node, start_node)
        prev_exon_start_node = start_node

    increment_edge_weight(node_details, prev_exon_start_node, read_end_node)


def remove_edges_based_on_read_counts(min_reads_per_edge, splice_graph):
    to_visit = set()
    visited = set()
    to_visit.add(splice_graph)
    while to_visit:
        node = to_visit.pop()
        visited.add(node)

        node_is_end = node.node_type == 'end'
        to_remove = list()
        for other, count in node.edges.items():
            if other not in visited:
                to_visit.add(other)

            if ((node_is_end and (other.node_type == 'start')
                 and (count < min_reads_per_edge))):
                to_remove.append(other)
                node.removed[other] = count

        for other in to_remove:
            remove_node_edge(node, other)


def initialize_extra_node_details(graph_results, extra_details, strand):
    orig_nodes = dict()
    extra_details['orig_nodes'] = orig_nodes
    source = graph_results['graph']
    orig_nodes[source] = source
    extra_details['source'] = source
    sink = graph_results['sink']
    orig_nodes[sink] = sink
    extra_details['sink'] = sink
    extra_details['extra_edges'] = dict()
    extra_nodes = dict()
    extra_details['extra_nodes'] = extra_nodes
    coordinates = list()
    extra_details['coords'] = coordinates
    for node in graph_results['closest_down']:
        if node in [source, sink]:
            continue

        coordinates.append(node.coordinate)
        orig_nodes[node] = node
        if node.node_type == 'start':
            other_type = 'end'
        else:
            other_type = 'start'

        other_node = Node(node.coordinate, other_type)
        extra_nodes[other_node] = other_node

    # negated_coords are used with bisect to determine
    # "mid" nodes for extra coordinates.
    # Bisect requires the list to be in ascending order.
    if strand == '-':
        coordinates.sort(reverse=True)
        extra_details['negated_coords'] = [-x for x in coordinates]
    else:
        coordinates.sort()
        extra_details['negated_coords'] = None

    mid_start_after_source = Node(None, 'mid_start')
    mid_end_after_source = Node(None, 'mid_end')
    extra_nodes[mid_start_after_source] = mid_start_after_source
    extra_nodes[mid_end_after_source] = mid_end_after_source
    for coord in coordinates:
        mid_start_after_coord = Node(coord, 'mid_start')
        mid_end_after_coord = Node(coord, 'mid_end')
        extra_nodes[mid_start_after_coord] = mid_start_after_coord
        extra_nodes[mid_end_after_coord] = mid_end_after_coord


def filter_splice_graphs(gene_i, min_reads_per_edge, align_handle_and_line,
                         splice_graph_results_plus,
                         splice_graph_results_minus):
    extra_details_plus = dict()
    extra_details_minus = dict()
    result = {'plus': extra_details_plus, 'minus': extra_details_minus}
    initialize_extra_node_details(splice_graph_results_plus,
                                  extra_details_plus, '+')
    initialize_extra_node_details(splice_graph_results_minus,
                                  extra_details_minus, '-')
    if not align_handle_and_line:
        return result

    line = align_handle_and_line.read_line()
    while line:
        columns = rmats_long_utils.read_tsv_line(line)
        gene_i_str = columns[0]
        found_gene_i = int(gene_i_str)
        if found_gene_i > gene_i:
            break

        line = align_handle_and_line.clear_and_read_line()
        if found_gene_i != gene_i:
            continue

        # sample = columns[2]
        start_str = columns[3]
        start = int(start_str)
        end_str = columns[4]
        end = int(end_str)
        sjs_str = columns[5]
        sjs = rmats_long_utils.parse_sjs_string(sjs_str)
        # The alignment strand is not used. The sequencing is assumed
        # to be unstranded and each alignment will be checked against
        # the splice graphs for both '+' and '-' strand transcripts.
        #
        # align_strand = columns[6]
        add_read_count_to_matching_edges(start, end, sjs, extra_details_plus)
        reversed_sjs = reverse_sjs(sjs)
        # start and end are swapped for the minus strand
        add_read_count_to_matching_edges(end, start, reversed_sjs,
                                         extra_details_minus)

    remove_edges_based_on_read_counts(min_reads_per_edge,
                                      splice_graph_results_plus['graph'])
    remove_edges_based_on_read_counts(min_reads_per_edge,
                                      splice_graph_results_minus['graph'])
    return result


def html_table_string(top=True,
                      bottom=True,
                      left=True,
                      right=True,
                      height='20',
                      width='20'):
    sides = ''
    if top:
        sides += 'T'
    if bottom:
        sides += 'B'
    if left:
        sides += 'L'
    if right:
        sides += 'R'

    return ('<<TABLE BORDER="0" CELLBORDER="2">'
            '<TR>'
            '<TD PORT="label" SIDES="{}" HEIGHT="{}" WIDTH="{}"></TD>'
            '</TR>'
            '</TABLE>>'.format(sides, height, width))


def write_line(handle, line):
    handle.write(line)
    handle.write('\n')


def add_double_quotes(string):
    return '"{}"'.format(string)


def get_node_name_and_label(node):
    if node.node_type in ['source', 'sink']:
        node_name = node.node_type
        node_label = add_double_quotes(node_name)
    elif node.node_type == 'start':
        node_name = 's{}'.format(node.coordinate)
        node_label = html_table_string(right=False)
    elif node.node_type == 'end':
        node_name = 'e{}'.format(node.coordinate)
        node_label = html_table_string(left=False)
    elif node.node_type == 'mid_start':
        node_name = 'ms{}'.format(node.coordinate)
        node_label = add_double_quotes('|')
    elif node.node_type == 'mid_end':
        node_name = 'me{}'.format(node.coordinate)
        node_label = add_double_quotes('|')
    else:
        raise Exception('unrecognized node_type: {}'.format(node.node_type))

    return node_name, node_label


def get_node_line(node):
    name, label = get_node_name_and_label(node)
    node_line = '{} [label = {}, shape = "none"];'.format(name, label)
    return node_line


def get_node_line_for_coord(coord, is_mid=False):
    if is_mid:
        label = '""'
    else:
        label = coord

    node_line = '{} [label = {}, shape = "none"];'.format(coord, label)
    return node_line


def append_edge_lines_for_node(node, edges, removed, extra_for_node,
                               show_counts, is_extra, node_lines):
    for next_node, weight in edges.items():
        edge_line = get_edge_line(node,
                                  next_node,
                                  weight,
                                  show_counts,
                                  is_removed=False,
                                  is_extra=is_extra)
        node_lines.append(edge_line)

    for next_node, weight in extra_for_node.items():
        edge_line = get_edge_line(node,
                                  next_node,
                                  weight,
                                  show_counts,
                                  is_removed=False,
                                  is_extra=True)
        node_lines.append(edge_line)

    for next_node, weight in removed.items():
        edge_line = get_edge_line(node,
                                  next_node,
                                  weight,
                                  show_counts,
                                  is_removed=True,
                                  is_extra=is_extra)
        node_lines.append(edge_line)


def get_from_source_edge_line(to_node, to_node_name, weight, is_removed,
                              is_extra):
    if to_node.node_type not in ['start', 'mid_start']:
        raise Exception(
            'source connects to a non-start node: {}'.format(to_node))

    style = ''
    if is_removed or is_extra:
        style = 'style = "dotted", '

    return ('source -> {}'
            ' [label = "{}", {}headport = "label:w",'
            ' headclip = "false"];'.format(to_node_name, weight, style))


def get_from_start_edge_line(from_node, from_node_name, to_node, to_node_name,
                             weight, is_removed, is_extra):
    if to_node.node_type not in ['end', 'mid_end']:
        raise Exception('a start node connects to a non-end node {} {}'.format(
            from_node, to_node))

    style = 'dashed'
    if is_removed or is_extra:
        style = 'dotted'

    return ('{} -> {} [label = "{}", style = "{}",'
            ' arrowhead = "none", headport = "label:e",'
            ' tailport = "label:w", headclip = "false",'
            ' tailclip = "false"];'.format(from_node_name, to_node_name,
                                           weight, style))


def get_from_end_edge_line(from_node, from_node_name, to_node, to_node_name,
                           weight, is_removed, is_extra):
    if to_node.node_type in ['end', 'mid_end']:
        raise Exception('end node connects to an end node {} {}'.format(
            from_node, to_node))

    style = ''
    if is_removed or is_extra:
        style = 'style = "dotted", '

    if to_node.node_type == 'sink':
        return ('{} -> sink'
                ' [label = "{}", {}tailport = "label:e",'
                ' tailclip = "false"];'.format(from_node_name, weight, style))

    return ('{} -> {}'
            ' [label = "{}", {}headport = "label:w",'
            ' tailport = "label:e", headclip = "false",'
            ' tailclip = "false"];'.format(from_node_name, to_node_name,
                                           weight, style))


def get_edge_line(from_node,
                  to_node,
                  weight,
                  show_counts,
                  is_removed=False,
                  is_extra=False):
    if not show_counts:
        weight = ''

    from_node_name, _ = get_node_name_and_label(from_node)
    to_node_name, _ = get_node_name_and_label(to_node)
    if from_node.node_type == 'source':
        return get_from_source_edge_line(to_node, to_node_name, weight,
                                         is_removed, is_extra)

    if from_node.node_type in ['start', 'mid_start']:
        return get_from_start_edge_line(from_node, from_node_name, to_node,
                                        to_node_name, weight, is_removed,
                                        is_extra)

    if from_node.node_type in ['end', 'mid_end']:
        return get_from_end_edge_line(from_node, from_node_name, to_node,
                                      to_node_name, weight, is_removed,
                                      is_extra)

    raise Exception('Unexpected from node type: {}'.format(
        from_node.node_type))


def get_rank_line(coords):
    coords_string = '; '.join(coords)
    return '{{ rank = same; {} }};'.format(coords_string)


def get_lines_based_on_coordinates_and_nodes_used(coordinates, lines_by_node):
    node_lines = list()
    rank_lines = list()
    coord_edge_parts = list()
    # None is used as the coordinate for the mid nodes next to the source.
    for coord in [None] + coordinates:
        mid_coord = 'mid_{}'.format(coord)
        regular_rank_parts = [str(coord)]
        mid_rank_parts = [mid_coord]
        regular_nodes = [Node(coord, 'start'), Node(coord, 'end')]
        mid_nodes = [Node(coord, 'mid_start'), Node(coord, 'mid_end')]

        regular_used = False
        for node in regular_nodes:
            if node not in lines_by_node:
                continue

            regular_used = True
            name, _ = get_node_name_and_label(node)
            regular_rank_parts.append(name)

        mid_used = False
        for node in mid_nodes:
            if node not in lines_by_node:
                continue

            mid_used = True
            name, _ = get_node_name_and_label(node)
            mid_rank_parts.append(name)

        if regular_used:
            coord_edge_parts.append(str(coord))
            node_line = get_node_line_for_coord(str(coord))
            node_lines.append(node_line)
            rank_line = get_rank_line(regular_rank_parts)
            rank_lines.append(rank_line)

        if mid_used:
            coord_edge_parts.append(mid_coord)
            node_line = get_node_line_for_coord(mid_coord, is_mid=True)
            node_lines.append(node_line)
            rank_line = get_rank_line(mid_rank_parts)
            rank_lines.append(rank_line)

    coord_order_line = None
    if coord_edge_parts:
        source_order = 'source_order'
        sink_order = 'sink_order'
        coord_edge_parts = [source_order] + coord_edge_parts + [sink_order]
        edge_details = ' -> '.join(coord_edge_parts)
        coord_order_line = '{} [style = "invis"];'.format(edge_details)
        source_order_line = get_node_line_for_coord(source_order, is_mid=True)
        sink_order_line = get_node_line_for_coord(sink_order, is_mid=True)
        node_lines.extend([source_order_line, sink_order_line])
        source_rank_line = get_rank_line(['source', source_order])
        sink_rank_line = get_rank_line(['sink', sink_order])
        rank_lines.extend([source_rank_line, sink_rank_line])

    return {
        'order_line': coord_order_line,
        'node_lines': node_lines,
        'rank_lines': rank_lines
    }


def output_graph_description(gene, extra_node_details, show_counts,
                             out_handle):
    coordinates = extra_node_details['coords']
    if not coordinates:
        return

    write_line(out_handle, 'strict digraph {} {{'.format(gene))
    write_line(out_handle, 'graph [rankdir = "LR", splines = "line"];')
    write_line(out_handle, 'edge [arrowhead = "empty", fontsize = "10"];')

    orig_nodes = extra_node_details['orig_nodes']
    extra_edges = extra_node_details['extra_edges']
    extra_nodes = extra_node_details['extra_nodes']
    sink_node = extra_node_details['sink']
    lines_by_node = dict()
    nodes_to_process = list(set(orig_nodes).union(set(extra_nodes)))
    for node in nodes_to_process:
        extra_for_node = extra_edges.get(node, dict())
        if not (node.edges or node.removed or extra_for_node):
            # sink node should not have any outgoing edges
            if node == sink_node:
                lines_by_node[node] = [get_node_line(node)]

            continue

        node_lines = list()
        lines_by_node[node] = node_lines
        node_line = get_node_line(node)
        node_lines.append(node_line)
        is_extra = node in extra_nodes
        append_edge_lines_for_node(node, node.edges, node.removed,
                                   extra_for_node, show_counts, is_extra,
                                   node_lines)

    coord_results = get_lines_based_on_coordinates_and_nodes_used(
        coordinates, lines_by_node)
    coord_results_lines = ([coord_results['order_line']]
                           + coord_results['node_lines']
                           + coord_results['rank_lines'])

    nodes_to_write = sorted(lines_by_node.keys())
    for node_to_write in nodes_to_write:
        lines = lines_by_node[node_to_write]
        for line in lines:
            write_line(out_handle, line)

    for line in coord_results_lines:
        if line is None:
            continue

        write_line(out_handle, line)

    write_line(out_handle, '}')


def find_transcripts_compatible_with_splice_graph(transcripts, strand, graph,
                                                  orig_nodes):
    new_transcripts = list()
    for transcript in transcripts:
        exons = transcript['exons']
        if strand == '-':
            exons = reverse_minus_strand_exons(exons)

        any_mismatch = False
        node = graph
        for exon in exons:
            start, end = exon
            start_node = Node(start, 'start')
            end_node = Node(end, 'end')
            start_node = orig_nodes.get(start_node)
            end_node = orig_nodes.get(end_node)
            if not (start_node and end_node):
                any_mismatch = True
                break

            if start_node not in node.edges:
                any_mismatch = True
                break

            if end_node not in start_node.edges:
                any_mismatch = True
                break

            node = end_node

        if not any_mismatch:
            new_transcripts.append(transcript)

    return new_transcripts


def detect_splicing_events_for_gene(gene, gene_i, all_transcripts,
                                    max_nodes_in_event, max_paths_in_event,
                                    min_reads_per_edge, output_full_gene_asm,
                                    simplify_gene_isoform_endpoints,
                                    filter_gene_isoforms_by_edge, show_counts,
                                    chr_id, chr_event_i, align_handle_and_line,
                                    graph_out_handle, out_handle):
    separated = separate_plus_and_minus_strand_transcripts(all_transcripts)
    events = list()

    splice_graph_results_plus = build_splice_graph(separated['plus'], '+')
    splice_graph_results_minus = build_splice_graph(separated['minus'], '-')
    extra_node_details = filter_splice_graphs(gene_i, min_reads_per_edge,
                                              align_handle_and_line,
                                              splice_graph_results_plus,
                                              splice_graph_results_minus)
    extra_details_plus = extra_node_details['plus']
    extra_details_minus = extra_node_details['minus']
    output_graph_description(gene, extra_details_plus, show_counts,
                             graph_out_handle)
    output_graph_description(gene, extra_details_minus, show_counts,
                             graph_out_handle)
    for strand in ['+', '-']:
        if strand == '+':
            results = splice_graph_results_plus
            transcripts = separated['plus']
            orig_nodes = extra_details_plus['orig_nodes']
        else:
            results = splice_graph_results_minus
            transcripts = separated['minus']
            orig_nodes = extra_details_minus['orig_nodes']

        if output_full_gene_asm:
            if filter_gene_isoforms_by_edge:
                transcripts = find_transcripts_compatible_with_splice_graph(
                    transcripts, strand, results['graph'], orig_nodes)

            bubbles = create_bubble_from_transcripts(transcripts, strand)
        else:
            bubbles = detect_bubbles(results['graph'], strand,
                                     max_nodes_in_event, max_paths_in_event)

        new_events = classify_bubbles(bubbles, results['closest_down'],
                                      results['closest_up'], results['graph'],
                                      strand, output_full_gene_asm,
                                      simplify_gene_isoform_endpoints)
        events.extend(new_events)

    events.sort(key=lambda event: event['type'])
    gene_id_value = gene
    for event_i, event in enumerate(events):
        asm_id = '{}_{}'.format(chr_id, chr_event_i)
        chr_event_i += 1
        columns = [
            asm_id, gene_i, gene_id_value, event['type'], event['strand'],
            event['start'], event['end']
        ]
        if event_i == 0:
            gene_id_value = ''

        isoforms = event['isoforms']
        isoform_strs = isoform_strings_from_isoforms(isoforms)
        columns.append(';'.join(isoform_strs))
        isoform_id_strs = list()
        isoform_ids = event['isoform_ids']
        for isoform_i, isoform_id in enumerate(isoform_ids):
            if isoform_id is None:
                isoform_id = '{}_{}'.format(asm_id, isoform_i)

            isoform_id_strs.append(isoform_id)

        columns.append(';'.join(isoform_id_strs))
        columns.extend([
            event['start_always_ss'], event['start_never_ss'],
            event['end_always_ss'], event['end_never_ss'], event['is_strict']
        ])
        rmats_long_utils.write_tsv_line(out_handle, [str(x) for x in columns])

    return chr_event_i


def isoform_strings_from_isoforms(isoforms):
    isoform_strs = list()
    for isoform in isoforms:
        coord_strs = list()
        for exon in isoform:
            start, end = exon
            coord_strs.append(str(start))
            coord_strs.append(str(end))

        isoform_str = ','.join(coord_strs)
        isoform_strs.append(isoform_str)

    return isoform_strs


def detect_and_write_events_thread_with_handles(
        max_nodes_in_event, max_paths_in_event, min_reads_per_edge,
        output_full_gene_asm, simplify_gene_isoform_endpoints,
        filter_gene_isoforms_by_edge, show_counts, chr_id, gtf_handle,
        align_handle, graph_out_handle, out_handle):
    headers = [
        'asm_id', 'gene_i', 'gene_id', 'event_type', 'strand', 'start', 'end',
        'isoforms', 'isoform_ids', 'start_always_ss', 'start_never_ss',
        'end_always_ss', 'end_never_ss', 'is_strict'
    ]
    rmats_long_utils.write_tsv_line(out_handle, headers)

    chr_event_i = 0
    gene = None
    gene_i = 0
    if align_handle:
        align_handle_and_line = rmats_long_utils.HandleAndNextLine(
            align_handle)
    else:
        align_handle_and_line = None

    transcripts = list()
    for line in gtf_handle:
        parsed = parse_gtf_line(line)
        if parsed.get('transcript_id'):
            details = dict()
            details['strand'] = parsed['strand']
            details['exons'] = parsed['exons']
            transcripts.append(parsed)
            continue

        if gene is not None:
            chr_event_i = detect_splicing_events_for_gene(
                gene, gene_i, transcripts, max_nodes_in_event,
                max_paths_in_event, min_reads_per_edge, output_full_gene_asm,
                simplify_gene_isoform_endpoints, filter_gene_isoforms_by_edge,
                show_counts, chr_id, chr_event_i, align_handle_and_line,
                graph_out_handle, out_handle)
            gene_i += 1

        gene = parsed['gene']
        transcripts = list()

    if gene is not None:
        detect_splicing_events_for_gene(
            gene, gene_i, transcripts, max_nodes_in_event, max_paths_in_event,
            min_reads_per_edge, output_full_gene_asm,
            simplify_gene_isoform_endpoints, filter_gene_isoforms_by_edge,
            show_counts, chr_id, chr_event_i, align_handle_and_line,
            graph_out_handle, out_handle)


def detect_and_write_events_thread(input_queue, max_nodes_in_event,
                                   max_paths_in_event, min_reads_per_edge,
                                   output_full_gene_asm,
                                   simplify_gene_isoform_endpoints,
                                   filter_gene_isoforms_by_edge, show_counts):
    while True:
        arguments = input_queue.get()
        if arguments is None:
            return

        gtf_path = arguments['gtf_path']
        align_path = arguments['align_path']
        graph_out_path = arguments['graph_out_path']
        out_path = arguments['out_path']
        chr_id = arguments['chr_id']
        with open(graph_out_path, 'wt') as graph_out_handle:
            with open(out_path, 'wt') as out_handle:
                with open(gtf_path, 'rt') as gtf_handle:
                    if align_path:
                        with open(align_path, 'rt') as align_handle:
                            detect_and_write_events_thread_with_handles(
                                max_nodes_in_event, max_paths_in_event,
                                min_reads_per_edge, output_full_gene_asm,
                                simplify_gene_isoform_endpoints,
                                filter_gene_isoforms_by_edge, show_counts,
                                chr_id, gtf_handle, align_handle,
                                graph_out_handle, out_handle)
                    else:
                        align_handle = None
                        detect_and_write_events_thread_with_handles(
                            max_nodes_in_event, max_paths_in_event,
                            min_reads_per_edge, output_full_gene_asm,
                            simplify_gene_isoform_endpoints,
                            filter_gene_isoforms_by_edge, show_counts, chr_id,
                            gtf_handle, align_handle, graph_out_handle,
                            out_handle)


def detect_splicing_events(max_nodes_in_event, max_paths_in_event, num_threads,
                           min_reads_per_edge, output_full_gene_asm,
                           simplify_gene_isoform_endpoints,
                           filter_gene_isoforms_by_edge, gtf_dir, align_dir,
                           out_dir):
    # TODO only apply max_nodes_in_event after reaching a timeout on the full search
    out_dir = os.path.abspath(out_dir)
    rmats_long_utils.create_output_dir(out_dir, check_empty=True)
    rmats_long_utils.copy_chr_name_mapping(gtf_dir, out_dir)
    show_counts = align_dir is not None
    threads = list()
    thread_inputs = multiprocessing.Queue(num_threads)
    for thread_i in range(num_threads):
        thread = multiprocessing.Process(
            target=detect_and_write_events_thread,
            args=(thread_inputs, max_nodes_in_event, max_paths_in_event,
                  min_reads_per_edge, output_full_gene_asm,
                  simplify_gene_isoform_endpoints,
                  filter_gene_isoforms_by_edge, show_counts))
        threads.append(thread)
        thread.start()

    gtf_file_names = os.listdir(gtf_dir)
    for file_name in gtf_file_names:
        if not file_name.startswith('chr_id_'):
            continue

        chr_id = rmats_long_utils.get_chr_id_from_path(file_name)
        gtf_path = os.path.join(gtf_dir, file_name)
        align_path = None
        if align_dir:
            align_path = os.path.join(align_dir, file_name)
            if not os.path.exists(align_path):
                align_path = None
                if min_reads_per_edge > 0:
                    # No events will pass the read cutoff
                    continue

        graph_out_name = 'graph_{}.txt'.format(chr_id)
        graph_out_path = os.path.join(out_dir, graph_out_name)
        out_path = os.path.join(out_dir, file_name)
        arguments = {
            'gtf_path': gtf_path,
            'align_path': align_path,
            'graph_out_path': graph_out_path,
            'out_path': out_path,
            'chr_id': chr_id
        }
        thread_inputs.put(arguments)

    # signal threads to stop
    for thread in threads:
        thread_inputs.put(None)

    for thread in threads:
        thread.join()
        if thread.exitcode != 0:
            raise Exception('thread exited with value: {}'.format(
                thread.exitcode))

    thread_inputs.close()


def main():
    args = parse_args()
    gtf_dir = os.path.abspath(args.gtf_dir)
    out_dir = os.path.abspath(args.out_dir)
    align_dir = None
    if args.align_dir:
        align_dir = os.path.abspath(args.align_dir)

    detect_splicing_events(args.max_nodes_in_event, args.max_paths_in_event,
                           args.num_threads, args.min_reads_per_edge,
                           args.output_full_gene_asm,
                           args.simplify_gene_isoform_endpoints,
                           args.filter_gene_isoforms_by_edge, gtf_dir,
                           align_dir, out_dir)
    print('detect_splicing_events.py finished')


if __name__ == '__main__':
    main()
