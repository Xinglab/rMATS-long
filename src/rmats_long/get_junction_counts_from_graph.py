import argparse
import os
import os.path

from rmats_long import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Create a .tsv with junction counts based on'
                     ' a graph_{id}.txt file'))
    parser.add_argument(
        '--graph-file',
        required=True,
        help='A graph_{id}.txt file from detect_splicing_events.py')
    parser.add_argument('--chr-name-id-mapping',
                        required=True,
                        help='The chr_name_id_mapping.tsv file')
    parser.add_argument('--out-tsv',
                        required=True,
                        help='Where to write junction counts')

    args = parser.parse_args()
    return args


def get_id_from_graph_file(graph_file):
    name = os.path.basename(graph_file)
    prefix = 'graph_'
    suffix = '.txt'
    if not (name.startswith(prefix) and name.endswith(suffix)):
        raise Exception('Unexpected format for {}'.format(graph_file))

    chr_id_str = name[len(prefix):-len(suffix)]
    chr_id = int(chr_id_str)
    return chr_id


# e123 -> s456 [label = "789", ...];
def parse_junction_line(line):
    tokens = line.split()
    if len(tokens) < 6:
        return None

    if [tokens[1], tokens[3], tokens[4]] != ['->', '[label', '=']:
        return None

    if not tokens[0].startswith('e'):
        return None

    if (not tokens[2].startswith('s')) or (tokens[2] == 'sink'):
        return None

    start_coord_str = tokens[0][1:]
    start_coord = int(start_coord_str)

    end_coord_str = tokens[2][1:]
    end_coord = int(end_coord_str)

    count_str = tokens[5].strip('",')
    count = int(count_str)
    return {'start': start_coord, 'end': end_coord, 'count': count}


def read_junction_counts(graph_file):
    count_by_junction = dict()
    with open(graph_file, 'rt') as in_handle:
        for line in in_handle:
            parsed_junction = parse_junction_line(line)
            if not parsed_junction:
                continue

            start = parsed_junction['start']
            end = parsed_junction['end']
            count = parsed_junction['count']
            if count == 0:
                continue

            junction = (start, end)
            old_count = count_by_junction.get(junction, 0)
            if count > old_count:
                count_by_junction[junction] = count

    return count_by_junction


def get_junction_counts_from_graph(graph_file, chr_name_id_mapping, out_tsv):
    chr_id = get_id_from_graph_file(graph_file)
    mapping_dir = os.path.dirname(chr_name_id_mapping)
    parsed_mapping = rmats_long_utils.parse_chr_mapping(mapping_dir)
    chr_name = parsed_mapping['id_to_chr'][chr_id]
    count_by_junction = read_junction_counts(graph_file)
    sorted_junctions = sorted(count_by_junction.keys())
    headers = ['chr', 'start', 'end', 'count']
    with open(out_tsv, 'wt') as out_handle:
        rmats_long_utils.write_tsv_line(out_handle, headers)
        for junction in sorted_junctions:
            start, end = junction
            count = count_by_junction[junction]
            columns = [chr_name, str(start), str(end), str(count)]
            rmats_long_utils.write_tsv_line(out_handle, columns)


def main():
    args = parse_args()
    get_junction_counts_from_graph(args.graph_file, args.chr_name_id_mapping,
                                   args.out_tsv)
    print('get_junction_counts_from_graph.py finished')


if __name__ == '__main__':
    main()
