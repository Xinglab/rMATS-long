import argparse


def parse_args():
    parser = argparse.ArgumentParser(description=(
        'Extract the sequence for an isoform from a GTF and FASTA'))
    parser.add_argument('--gtf',
                        required=True,
                        help='.gtf file with coordinates')
    parser.add_argument('--fasta',
                        required=True,
                        help='.fa file with sequence')
    parser.add_argument('--transcript-id',
                        required=True,
                        help='the ID used to lookup the isoform in the GTF')

    return parser.parse_args()


def trim_quotes(string):
    return string.strip('"')


def read_coordinates_from_gtf(path, transcript_id):
    coordinates = list()
    strand = None
    chr_name = None
    with open(path, 'rt') as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('#'):
                continue

            columns = line.split('\t')
            feature = columns[2]
            if feature != 'exon':
                continue

            attribute_string = columns[8]
            attrs = attribute_string.split(';')
            for attr in attrs:
                attr = attr.strip()
                if not attr:
                    continue

                parts = attr.split()
                if parts[0] == 'transcript_id':
                    found_id = trim_quotes(parts[1])

            if found_id != transcript_id:
                continue

            row_chr_name = columns[0]
            if chr_name is None:
                chr_name = row_chr_name
            elif row_chr_name != chr_name:
                raise Exception('found conflicting chr info: {} != {}'.format(
                    row_chr_name, chr_name))
            start_str = columns[3]
            start = int(start_str)
            end_str = columns[4]
            end = int(end_str)
            row_strand = columns[6]
            if strand is None:
                strand = row_strand
            elif row_strand != strand:
                raise Exception(
                    'found conflicting strand info: {} != {}'.format(
                        row_strand, strand))

            # gtf could have multiple entries for the same transcript
            exon = [start, end]
            if exon not in coordinates:
                coordinates.append(exon)

    return {'coords': coordinates, 'strand': strand, 'chr': chr_name}


def read_sequence_from_fasta(coord_details, path):
    coordinates = sorted(coord_details['coords'])
    strand = coord_details['strand']
    chr_name = coord_details['chr']
    min_coord = coordinates[0][0]
    max_coord = coordinates[-1][1]
    target_len = (max_coord - min_coord) + 1
    sequence_list_from_min = list()
    len_copied = 0
    found_chr = False
    coords_to_skip = min_coord - 1
    with open(path, 'rt') as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not found_chr:
                if not line.startswith('>'):
                    continue

                parts = line.split()
                found_chr = parts[0] == '>{}'.format(chr_name)
                continue

            line_len = len(line)
            if line_len < coords_to_skip:
                coords_to_skip -= line_len
                continue

            if coords_to_skip:
                line_seq = line[coords_to_skip:]
                coords_to_skip = 0
            else:
                line_seq = line

            sequence_list_from_min.append(line_seq)
            len_copied += len(line_seq)
            if len_copied >= target_len:
                break

    if strand == '-':
        exon_coords = sorted(coordinates, reverse=True)
    else:
        exon_coords = coordinates

    sequence_from_min = ''.join(sequence_list_from_min)
    sequence_details = list()
    for i, exon_coords in enumerate(exon_coords):
        name = 'exon_{}'.format(i + 1)
        start, end = exon_coords
        adjusted_start = start - min_coord
        adjusted_end = end - min_coord
        sequence = sequence_from_min[adjusted_start:adjusted_end + 1]
        if strand == '-':
            sequence = reverse_complement(sequence)

        sequence_details.append({'name': name, 'sequence': sequence})

    return sequence_details


def reverse_complement(sequence):
    new_seq = list()
    for char in sequence:
        char = char.upper()
        if char == 'A':
            new_seq.append('T')
        elif char == 'C':
            new_seq.append('G')
        elif char == 'G':
            new_seq.append('C')
        elif char == 'T':
            new_seq.append('A')
        elif char == 'N':
            new_seq.append('N')
        else:
            raise Exception('unrecognized sequence: {} from {}'.format(
                char, sequence))

    return ''.join(reversed(new_seq))


def print_sequence(details):
    for detail in details:
        print('{}:'.format(detail['name']))
        print(detail['sequence'])
        print()


def main():
    args = parse_args()
    gtf_path = args.gtf
    fasta_path = args.fasta
    transcript_id = args.transcript_id
    coord_details = read_coordinates_from_gtf(gtf_path, transcript_id)
    sequence_details = read_sequence_from_fasta(coord_details, fasta_path)
    print_sequence(sequence_details)


if __name__ == '__main__':
    main()
