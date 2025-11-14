import argparse
import multiprocessing
import os
import os.path
import subprocess
import tempfile

import rmats_long_utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Determine read counts for isoforms in ASMs')
    parser.add_argument(
        '--event-dir',
        required=True,
        help='The output directory from detect_splicing_events.py')
    parser.add_argument('--align-dir',
                        required=True,
                        help=('The output directory from'
                              ' organize_alignment_info_by_gene_and_chr.py'))
    parser.add_argument(
        '--gtf-dir',
        required=True,
        help='The output directory from organize_gene_info_by_chr.py')
    parser.add_argument('--out-tsv',
                        required=True,
                        help='The output file for ASM read counts')
    parser.add_argument('--num-threads',
                        default=1,
                        type=int,
                        help='how many threads to use')
    parser.add_argument('--sort-buffer-size',
                        default='2G',
                        help=('Used for the --buffer-size argument of sort.'
                              ' Default: %(default)s'))

    return parser.parse_args()


def read_next_gene_id(gtf_handle):
    gene_id = None
    line = gtf_handle.readline()
    while line:
        columns = rmats_long_utils.read_tsv_line(line)
        if len(columns) != 4:
            line = gtf_handle.readline()
            continue

        gene_id = columns[2]
        break

    return gene_id


def load_asms_for_gene_id(gene_id, event_handle_and_line):
    asms = list()
    line = event_handle_and_line.read_line()
    expected_headers = [
        'asm_id', 'gene_i', 'gene_id', 'event_type', 'strand', 'start', 'end',
        'isoforms', 'isoform_ids', 'start_always_ss', 'start_never_ss',
        'end_always_ss', 'end_never_ss', 'is_strict'
    ]
    matched_the_gene = False
    while line:
        columns = rmats_long_utils.read_tsv_line(line)
        if columns == expected_headers:
            line = event_handle_and_line.clear_and_read_line()
            continue

        found_gene_id = columns[2]
        if not matched_the_gene:
            if found_gene_id == gene_id:
                matched_the_gene = True
            else:
                break
        elif found_gene_id != '':
            break

        line = event_handle_and_line.clear_and_read_line()
        asm_id = columns[0]
        strand = columns[4]
        start = columns[5]
        end = columns[6]
        isoforms_str = columns[7]
        isoforms = rmats_long_utils.parse_isoforms_str(isoforms_str)
        isoform_ids_str = columns[8]
        isoform_ids = rmats_long_utils.parse_isoform_ids_str(isoform_ids_str)
        start_always_ss = rmats_long_utils.parse_bool(columns[9])
        start_never_ss = rmats_long_utils.parse_bool(columns[10])
        end_always_ss = rmats_long_utils.parse_bool(columns[11])
        end_never_ss = rmats_long_utils.parse_bool(columns[12])
        isoform_details = create_isoform_details(isoforms, start, end, strand,
                                                 start_always_ss,
                                                 start_never_ss, end_always_ss,
                                                 end_never_ss)
        asm = {
            'asm_id': asm_id,
            'strand': strand,
            'isoforms': isoform_details,
            'isoform_ids': isoform_ids,
        }
        asms.append(asm)

    return asms


def determine_start_type(isoforms, start, start_is_source, start_always_ss,
                         start_never_ss):
    if start_is_source:
        return 'isoform_start'

    start = int(start)
    first_isoform_start = isoforms[0][0]
    if first_isoform_start != start:
        return 'any_exon_start'
    if start_always_ss:
        return 'exon_start'
    if start_never_ss:
        return 'isoform_start'

    return 'exon_or_isoform_start'


def determine_end_type(isoforms, end, end_is_sink, end_always_ss,
                       end_never_ss):
    if end_is_sink:
        return 'isoform_end'

    end = int(end)
    last_isoform_end = isoforms[-1][-1]
    if last_isoform_end != end:
        return 'any_exon_end'
    if end_always_ss:
        return 'exon_end'
    if end_never_ss:
        return 'isoform_end'

    return 'exon_or_isoform_end'


# The isoform details for a minus strand isoform will be in ascending
# coordinate order like a plus strand isoform. The exon start types
# will actually be at exon end coordinates. The exon end types will be
# at actual exon start coords. The start and end columns from the ASM
# definition need to be interpreted based on the strand.
#
# Having minus strand isoforms treated this way simplifies the
# compatibility checks.
def create_isoform_details(isoforms, start, end, strand, start_always_ss,
                           start_never_ss, end_always_ss, end_never_ss):
    if strand == '-':
        start_is_source = end == 'sink'
        start_type = determine_start_type(isoforms, end, start_is_source,
                                          end_always_ss, end_never_ss)
        end_is_sink = start == 'source'
        end_type = determine_end_type(isoforms, start, end_is_sink,
                                      start_always_ss, start_never_ss)
    else:
        start_is_source = start == 'source'
        start_type = determine_start_type(isoforms, start, start_is_source,
                                          start_always_ss, start_never_ss)
        end_is_sink = end == 'sink'
        end_type = determine_end_type(isoforms, end, end_is_sink,
                                      end_always_ss, end_never_ss)

    isoform_details = list()
    for isoform in isoforms:
        coords = list()
        details = {'coords': coords}
        isoform_details.append(details)
        is_exon_start = True
        for coordinate in isoform:
            if is_exon_start:
                coord_type = 'exon_start'
                is_exon_start = False
            else:
                coord_type = 'exon_end'
                is_exon_start = True

            coords.append({'value': coordinate, 'type': coord_type})

        coords[0]['type'] = start_type
        coords[-1]['type'] = end_type

    return isoform_details


def isoform_coord_type_is_start(coord_type):
    return coord_type in [
        'exon_start', 'isoform_start', 'exon_or_isoform_start',
        'any_exon_start'
    ]


def isoform_coord_type_is_end(coord_type):
    return coord_type in [
        'exon_end', 'isoform_end', 'exon_or_isoform_end', 'any_exon_end'
    ]


class ReadPositionDetails():
    def __init__(self, start, end, sjs):
        self.start = start
        self.end = end
        self.sjs = sjs
        self.pos = start
        self.prev_pos = None
        self.num_sjs = len(sjs)
        self.sj_i = -1  # -1 is read start and num_sjs is read_end
        self.i_within_sj = None

    def is_at_start(self):
        return self.sj_i == -1

    def is_at_end(self):
        return self.sj_i == self.num_sjs

    def is_at_exon_start(self):
        if self.is_at_start():
            return True

        if self.is_at_end():
            return False

        return self.i_within_sj == 1

    def is_at_exon_end(self):
        return not self.is_at_exon_start()

    def is_at_first_exon_end(self):
        return self.sj_i == 0 and self.i_within_sj == 0

    def _set_pos_from_sj_offsets(self):
        if self.is_at_start():
            self.pos = self.start
            return

        if self.is_at_end():
            self.pos = self.end
            return

        self.pos = self.sjs[self.sj_i][self.i_within_sj]

    def advance(self):
        if self.is_at_end():
            raise Exception('Cannot advance past end')

        self.prev_pos = self.pos
        if self.is_at_start():
            self.sj_i = 0
            self.i_within_sj = 0
        elif self.i_within_sj == 0:
            self.i_within_sj = 1
        else:
            self.sj_i += 1
            self.i_within_sj = 0

        self._set_pos_from_sj_offsets()

    def __str__(self):
        return ('start: {}, end: {}, pos: {}, prev_pos: {}, num_sjs: {},'
                ' sj_i: {}, i_within_sj: {}, sjs: {}'.format(
                    self.start, self.end, self.pos, self.prev_pos,
                    self.num_sjs, self.sj_i, self.i_within_sj, self.sjs))


class IsoformPositionDetails():
    def __init__(self, isoform):
        self.coords = isoform['coords']
        self.coord_i = 0
        self.coord = self.coords[self.coord_i]
        self.pos = self.coord['value']
        self.type = self.coord['type']
        self.num_coords = len(self.coords)

    def is_at_start(self):
        return self.coord_i == 0

    def is_at_end(self):
        return self.coord_i == (self.num_coords - 1)

    def prev_pos(self):
        if self.is_at_start():
            raise Exception('Nothing before start')

        return self.coords[self.coord_i - 1]['value']

    def prev_type(self):
        if self.is_at_start():
            raise Exception('Nothing before start')

        return self.coords[self.coord_i - 1]['type']

    def advance(self):
        if self.is_at_end():
            raise Exception('Cannot advance past end')

        self.coord_i += 1
        self.coord = self.coords[self.coord_i]
        self.pos = self.coord['value']
        self.type = self.coord['type']

    def __str__(self):
        return ('coord_i: {}, pos: {}, type: {}, num_coords: {}, coords: {}'
                .format(self.coord_i, self.pos, self.type, self.num_coords,
                        self.coords))


# After this is run, at least one of the read or the isoform will
# still be at start.  If they have the same start then both will still
# be at start.  Otherwise the one not at start will be advanced until
# it hits or crosses the other start.
def find_read_isoform_initial_overlap(read_pos_details, isoform_pos_details):
    if read_pos_details.pos == isoform_pos_details.pos:
        return True

    if read_pos_details.pos < isoform_pos_details.pos:
        while not read_pos_details.is_at_end():
            read_pos_details.advance()
            if read_pos_details.pos >= isoform_pos_details.pos:
                return True

        return False

    while not isoform_pos_details.is_at_end():
        isoform_pos_details.advance()
        if isoform_pos_details.pos >= read_pos_details.pos:
            return True

    return False


def check_start_boundary(read_pos_details, isoform_pos_details):
    if isoform_pos_details.is_at_start():
        if read_pos_details.pos == read_pos_details.start:
            return True

        if read_pos_details.is_at_exon_start():
            # read spliced to 1st isoform exon start
            if read_pos_details.pos == isoform_pos_details.pos:
                return isoform_pos_details.type != 'isoform_start'

            # read spliced into 1st isoform exon
            return isoform_pos_details.type == 'any_exon_start'

        # read crosses isoform_start
        if isoform_pos_details.type == 'exon_start':
            return False

        if isoform_pos_details.type == 'any_exon_start':
            return True

        if isoform_pos_details.type in [
                'isoform_start', 'exon_or_isoform_start'
        ]:
            # read 1st exon reads over isoform start
            return read_pos_details.is_at_first_exon_end()

        raise Exception('Unexpected isoform start type: {}'.format(
            isoform_pos_details.type))

    # read starts within isoform
    if isoform_coord_type_is_end(isoform_pos_details.type):
        # read starts within an exon
        return True

    if isoform_pos_details.pos == read_pos_details.pos:
        # read starts at exon start
        return True

    # read starts in an intron
    return False


def check_end_boundary(read_pos_details, isoform_pos_details):
    if read_pos_details.is_at_end():
        if isoform_pos_details.is_at_end():
            if read_pos_details.pos <= isoform_pos_details.pos:
                return True

            # read crossed the isoform end
            return isoform_pos_details.type != 'exon_end'

        # read ends within internal exon of isoform
        if read_pos_details.pos <= isoform_pos_details.pos:
            return True

        # read crossed internal exon boundary
        return False

    if isoform_pos_details.is_at_end():
        # isoform ends, but read continues to another exon
        if isoform_pos_details.type == 'any_exon_end':
            return True

        if isoform_pos_details.type in ['exon_end', 'exon_or_isoform_end']:
            return read_pos_details.pos == isoform_pos_details.pos

        return False

    read_at_exon_start = read_pos_details.is_at_exon_start()
    isoform_at_exon_start = isoform_coord_type_is_start(
        isoform_pos_details.type)
    if read_at_exon_start != isoform_at_exon_start:
        raise Exception(
            'read and isoform should be at the same boundary type'
            ' read_at_exon_start: {}, isoform_at_exon_start: {}'.format(
                read_at_exon_start, isoform_at_exon_start))

    # different internal exon boundary
    return False


def is_read_compatible_with_isoform(start, end, sjs, isoform):
    read_pos_details = ReadPositionDetails(start, end, sjs)
    isoform_pos_details = IsoformPositionDetails(isoform)
    had_overlap = find_read_isoform_initial_overlap(read_pos_details,
                                                    isoform_pos_details)
    if not had_overlap:
        return False

    start_boundary_is_ok = check_start_boundary(read_pos_details,
                                                isoform_pos_details)
    if not start_boundary_is_ok:
        return False

    # The read and isoform coordinates have been advanced until the first overlap.
    # If only one is at an exon start then advance it to the exon end.
    if isoform_coord_type_is_start(isoform_pos_details.type):
        if read_pos_details.is_at_exon_end():
            isoform_pos_details.advance()
    elif read_pos_details.is_at_exon_start():
        read_pos_details.advance()

    # The read and isoform are now both at an exon start or both at an exon end.
    # Advance along both the read and isoform while the coords match.
    while read_pos_details.pos == isoform_pos_details.pos:
        if read_pos_details.is_at_end() or isoform_pos_details.is_at_end():
            break

        read_pos_details.advance()
        isoform_pos_details.advance()

    return check_end_boundary(read_pos_details, isoform_pos_details)


def count_read_for_asm(gene_id, read_i, sample, start, end, sjs, asm,
                       out_handle):
    asm_id = asm['asm_id']
    asm_id_parts = asm_id.split('_')
    chr_id = asm_id_parts[0]
    read_id = '{}_{}'.format(chr_id, read_i)
    isoforms = asm['isoforms']
    isoform_ids = asm['isoform_ids']
    # The last column is the isoform ID
    out_columns = [asm_id, gene_id, read_id, sample, None]
    for isoform_i, isoform in enumerate(isoforms):
        is_match = is_read_compatible_with_isoform(start, end, sjs, isoform)
        if is_match:
            isoform_id = isoform_ids[isoform_i]
            out_columns[4] = isoform_id
            rmats_long_utils.write_tsv_line(out_handle, out_columns)


def count_reads_for_gene(gene_i, gene_id, read_i, asms, align_handle_and_line,
                         out_handle):
    line = align_handle_and_line.read_line()
    while line:
        columns = rmats_long_utils.read_tsv_line(line)
        gene_i_str = columns[0]
        found_gene_i = int(gene_i_str)
        if found_gene_i > gene_i:
            break

        line = align_handle_and_line.clear_and_read_line()
        read_i += 1
        # gene_id = columns[1]
        sample = columns[2]
        start_str = columns[3]
        start = int(start_str)
        end_str = columns[4]
        end = int(end_str)
        sjs_str = columns[5]
        sjs = rmats_long_utils.parse_sjs_string(sjs_str)
        # strand = columns[6]
        for asm in asms:
            count_read_for_asm(gene_id, read_i, sample, start, end, sjs, asm,
                               out_handle)

    return read_i


def count_reads_for_asms_thread_with_handles(gtf_handle, event_handle,
                                             align_handle, out_handle):
    gene_i = 0
    read_i = 0
    event_handle_and_line = rmats_long_utils.HandleAndNextLine(event_handle)
    align_handle_and_line = rmats_long_utils.HandleAndNextLine(align_handle)
    while True:
        gene_id = read_next_gene_id(gtf_handle)
        if gene_id is None:
            break

        asms = load_asms_for_gene_id(gene_id, event_handle_and_line)
        read_i = count_reads_for_gene(gene_i, gene_id, read_i, asms,
                                      align_handle_and_line, out_handle)
        gene_i += 1


def count_reads_for_asms_thread(input_queue):
    while True:
        arguments = input_queue.get()
        if arguments is None:
            return

        gtf_path = arguments['gtf_path']
        event_path = arguments['event_path']
        align_path = arguments['align_path']
        out_path = arguments['out_path']
        with open(out_path, 'wt') as out_handle:
            with open(gtf_path, 'rt') as gtf_handle:
                with open(event_path, 'rt') as event_handle:
                    with open(align_path, 'rt') as align_handle:
                        count_reads_for_asms_thread_with_handles(
                            gtf_handle, event_handle, align_handle, out_handle)


# A read_id from the original alignment/fastq can appear in multiple genes.
# The original read_id will not be preserved. Instead it will have a different
# '{chr_i}_{read_i}' for each gene.
def combine_out_files(buffer_size, out_paths, out_tsv):
    with open(out_tsv, 'wt') as out_handle:
        for path in out_paths:
            with open(path, 'rt') as in_handle:
                for line in in_handle:
                    out_handle.write(line)

    dir_path = os.path.dirname(os.path.abspath(out_tsv))
    with tempfile.TemporaryDirectory(suffix='_tmp',
                                     prefix='combine_out_files',
                                     dir=dir_path) as temp_dir:
        tmp_path = os.path.join(temp_dir, 'sorted.tmp')
        sort_combined_file(out_tsv, tmp_path, temp_dir, buffer_size)
        with open(out_tsv, 'wt') as out_handle:
            headers = [
                'asm_id', 'gene_id', 'read_id', 'sample_id', 'isoform_id'
            ]
            rmats_long_utils.write_tsv_line(out_handle, headers)
            with open(tmp_path, 'rt') as in_handle:
                prev_asm = None
                for line in in_handle:
                    columns = rmats_long_utils.read_tsv_line(line)
                    asm = columns[0]
                    if (prev_asm is not None) and (asm == prev_asm):
                        # only report the gene_id once per ASM
                        columns[1] = ''

                    prev_asm = asm
                    rmats_long_utils.write_tsv_line(out_handle, columns)


def sort_combined_file(path, tmp_path, temp_dir, buffer_size):
    asm_key_arg = '-k1,1'
    read_key_arg = '-k3,3'
    env = {'LC_ALL': 'C'}  # to ensure sort order
    command = [
        'sort', '--buffer-size', buffer_size, '--temporary-directory',
        temp_dir, asm_key_arg, read_key_arg, '--output', tmp_path, path
    ]
    subprocess.run(command, env=env, check=True)


def count_reads_for_asms(buffer_size, event_dir, align_dir, gtf_dir,
                         num_threads, out_tsv):
    out_dir = os.path.dirname(os.path.abspath(out_tsv))
    gtf_file_names = os.listdir(gtf_dir)
    out_paths = list()
    threads = list()
    thread_inputs = multiprocessing.Queue(len(gtf_file_names))
    for thread_i in range(num_threads):
        thread = multiprocessing.Process(target=count_reads_for_asms_thread,
                                         args=(thread_inputs, ))
        threads.append(thread)
        thread.start()

    with tempfile.TemporaryDirectory(suffix='_tmp',
                                     prefix='count_reads_for_asms',
                                     dir=out_dir) as temp_dir:
        for file_name in gtf_file_names:
            if not file_name.startswith('chr_id_'):
                continue

            gtf_path = os.path.join(gtf_dir, file_name)
            event_path = os.path.join(event_dir, file_name)
            align_path = os.path.join(align_dir, file_name)
            out_path = os.path.join(temp_dir, file_name)
            # Some chrs may not have events or alignments
            if not (os.path.exists(event_path) and os.path.exists(align_path)):
                continue

            out_paths.append(out_path)
            arguments = {
                'gtf_path': gtf_path,
                'event_path': event_path,
                'align_path': align_path,
                'out_path': out_path
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
        combine_out_files(buffer_size, out_paths, out_tsv)


def main():
    args = parse_args()
    gtf_dir = os.path.abspath(args.gtf_dir)
    align_dir = os.path.abspath(args.align_dir)
    event_dir = os.path.abspath(args.event_dir)
    count_reads_for_asms(args.sort_buffer_size, event_dir, align_dir, gtf_dir,
                         args.num_threads, args.out_tsv)
    print('count_reads_for_asms.py finished')


if __name__ == '__main__':
    main()
