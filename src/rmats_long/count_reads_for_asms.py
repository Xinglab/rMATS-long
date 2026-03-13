import argparse
import multiprocessing
import os
import os.path
import shutil
import subprocess
import tempfile
import time

import cython
from cython.cimports.libc.stdint import int8_t, int64_t
from cython.cimports.libc.stdlib import free, malloc

from rmats_long import rmats_long_utils


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
    parser.add_argument(
        '--out-dir',
        required=True,
        help='The directory to write ASM read counts in output files by chr')
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


Int64Pair = cython.typedef(int64_t[2])
SpliceJunction = cython.typedef(Int64Pair)

# sj_i: -1 is read start, and num_sjs is read_end
ReadPositionDetails = cython.struct(start=int64_t,
                                    end=int64_t,
                                    sjs=cython.pointer[SpliceJunction],
                                    pos=int64_t,
                                    num_sjs=cython.size_t,
                                    sj_i=int64_t,
                                    i_within_sj=int8_t)
p_ReadPositionDetails = cython.typedef(cython.pointer[ReadPositionDetails])


@cython.cfunc
@cython.exceptval(check=True)
def read_position_details_init(in_start: int, in_end: int,
                               in_sjs: list) -> ReadPositionDetails:
    num_sjs: cython.size_t
    malloc_size: cython.size_t
    sjs: cython.pointer[SpliceJunction]
    sj_start: int
    sj_end: int
    sj_i: int64_t
    pos: int64_t
    out_sj_i: int64_t
    i_within_sj: int8_t
    num_sjs = len(in_sjs)
    malloc_size = num_sjs * cython.sizeof(SpliceJunction)
    sjs = cython.cast(cython.pointer[SpliceJunction], malloc(malloc_size))
    if sjs is cython.NULL:
        raise Exception('read_position_details_init: null return from malloc')

    for sj_i in range(num_sjs):
        sj_start = in_sjs[sj_i][0]
        sj_end = in_sjs[sj_i][1]
        sjs[sj_i][0] = sj_start
        sjs[sj_i][1] = sj_end

    pos = in_start
    out_sj_i = -1
    i_within_sj = -1
    return ReadPositionDetails(in_start, in_end, sjs, pos, num_sjs, out_sj_i,
                               i_within_sj)


@cython.cfunc
@cython.exceptval(check=False)
def read_position_details_free(rpd: p_ReadPositionDetails) -> cython.void:
    free(rpd.sjs)


@cython.cfunc
@cython.exceptval(check=False)
def read_position_details_reset(rpd: p_ReadPositionDetails) -> cython.void:
    rpd.pos = rpd.start
    rpd.sj_i = -1
    rpd.i_within_sj = -1


@cython.cfunc
@cython.exceptval(check=False)
def read_position_details_is_at_start(
        rpd: p_ReadPositionDetails) -> cython.bint:
    return rpd.sj_i == -1


@cython.cfunc
@cython.exceptval(check=False)
def read_position_details_is_at_end(rpd: p_ReadPositionDetails) -> cython.bint:
    return rpd.sj_i == rpd.num_sjs


@cython.cfunc
@cython.exceptval(check=False)
def read_position_details_is_at_exon_start(
        rpd: p_ReadPositionDetails) -> cython.bint:
    if read_position_details_is_at_start(rpd):
        return True

    if read_position_details_is_at_end(rpd):
        return False

    return rpd.i_within_sj == 1


@cython.cfunc
@cython.exceptval(check=False)
def read_position_details_is_at_exon_end(
        rpd: p_ReadPositionDetails) -> cython.bint:
    return not read_position_details_is_at_exon_start(rpd)


@cython.cfunc
@cython.exceptval(check=False)
def read_position_details_is_at_first_exon_end(
        rpd: p_ReadPositionDetails) -> cython.bint:
    return rpd.sj_i == 0 and rpd.i_within_sj == 0


@cython.cfunc
@cython.exceptval(check=False)
def read_position_details_set_pos_from_sj_offsets(
        rpd: p_ReadPositionDetails) -> cython.void:
    if read_position_details_is_at_start(rpd):
        rpd.pos = rpd.start
        return

    if read_position_details_is_at_end(rpd):
        rpd.pos = rpd.end
        return

    rpd.pos = rpd.sjs[rpd.sj_i][rpd.i_within_sj]


@cython.cfunc
@cython.exceptval(check=False)
def read_position_details_advance(rpd: p_ReadPositionDetails) -> cython.void:
    if read_position_details_is_at_start(rpd):
        rpd.sj_i = 0
        rpd.i_within_sj = 0
    elif rpd.i_within_sj == 0:
        rpd.i_within_sj = 1
    else:
        rpd.sj_i += 1
        rpd.i_within_sj = 0

    read_position_details_set_pos_from_sj_offsets(rpd)


CoordTypeEnumValueType = cython.typedef(int8_t)
CoordTypeEnumType = cython.struct(exon_start=CoordTypeEnumValueType,
                                  exon_end=CoordTypeEnumValueType,
                                  isoform_start=CoordTypeEnumValueType,
                                  isoform_end=CoordTypeEnumValueType,
                                  exon_or_isoform_start=CoordTypeEnumValueType,
                                  exon_or_isoform_end=CoordTypeEnumValueType,
                                  any_exon_start=CoordTypeEnumValueType,
                                  any_exon_end=CoordTypeEnumValueType)


@cython.cfunc
@cython.exceptval(check=False)
def get_coord_type_enum() -> CoordTypeEnumType:
    coord_type_enum: CoordTypeEnumType
    coord_type_enum = CoordTypeEnumType(
        1,  # exon_start
        2,  # exon_end
        3,  # isoform_start
        4,  # isoform_end
        5,  # exon_or_isoform_start
        6,  # exon_or_isoform_end
        7,  # any_exon_start
        8  # any_exon_end
    )
    return coord_type_enum


@cython.cfunc
@cython.exceptval(check=True)
def convert_to_coord_type_enum(string: str) -> CoordTypeEnumValueType:
    enum: CoordTypeEnumType
    enum = get_coord_type_enum()
    if string == 'exon_start':
        return enum.exon_start
    if string == 'exon_end':
        return enum.exon_end
    if string == 'isoform_start':
        return enum.isoform_start
    if string == 'isoform_end':
        return enum.isoform_end
    if string == 'exon_or_isoform_start':
        return enum.exon_or_isoform_start
    if string == 'exon_or_isoform_end':
        return enum.exon_or_isoform_end
    if string == 'any_exon_start':
        return enum.any_exon_start
    if string == 'any_exon_end':
        return enum.any_exon_end

    raise Exception(
        'convert_to_coord_type_enum: {} not in enum'.format(string))


@cython.cfunc
@cython.exceptval(check=False)
def isoform_coord_type_is_start(
        coord_type: CoordTypeEnumValueType) -> cython.bint:
    enum: CoordTypeEnumType
    enum = get_coord_type_enum()
    return coord_type in [
        enum.exon_start, enum.isoform_start, enum.exon_or_isoform_start,
        enum.any_exon_start
    ]


@cython.cfunc
@cython.exceptval(check=False)
def isoform_coord_type_is_end(
        coord_type: CoordTypeEnumValueType) -> cython.bint:
    enum: CoordTypeEnumType
    enum = get_coord_type_enum()
    return coord_type in [
        enum.exon_end, enum.isoform_end, enum.exon_or_isoform_end,
        enum.any_exon_end
    ]


IsoformPositionDetails = cython.struct(
    coords=cython.pointer[int64_t],
    coord_types=cython.pointer[CoordTypeEnumValueType],
    coord_i=int64_t,
    pos=int64_t,
    coord_type=CoordTypeEnumValueType,
    num_coords=cython.size_t)
p_IsoformPositionDetails = cython.typedef(
    cython.pointer[IsoformPositionDetails])
pp_IsoformPositionDetails = cython.typedef(
    cython.pointer[p_IsoformPositionDetails])


@cython.cfunc
@cython.exceptval(check=True)
def isoform_position_details_init(isoform: dict) -> IsoformPositionDetails:
    isoform_coords: list
    isoform_coord: dict
    num_coords: cython.size_t
    malloc_size: cython.size_t
    coords: cython.pointer[int64_t]
    coord_types: cython.pointer[CoordTypeEnumValueType]
    coord_i: int64_t
    pos: int64_t
    coord_type: CoordTypeEnumValueType
    isoform_coords = isoform['coords']
    num_coords = len(isoform_coords)
    malloc_size = num_coords * cython.sizeof(int64_t)
    coords = cython.cast(cython.pointer[int64_t], malloc(malloc_size))
    if coords is cython.NULL:
        raise Exception(
            'isoform_position_details_init: null return from malloc coords')

    malloc_size = num_coords * cython.sizeof(CoordTypeEnumValueType)
    coord_types = cython.cast(cython.pointer[CoordTypeEnumValueType],
                              malloc(malloc_size))
    if coord_types is cython.NULL:
        raise Exception(
            'isoform_position_details_init: null return from malloc coord_types'
        )

    for coord_i in range(num_coords):
        isoform_coord = isoform_coords[coord_i]
        coords[coord_i] = isoform_coord['value']
        coord_types[coord_i] = convert_to_coord_type_enum(
            isoform_coord['type'])

    coord_i = 0
    pos = coords[coord_i]
    coord_type = coord_types[coord_i]
    return IsoformPositionDetails(coords, coord_types, coord_i, pos,
                                  coord_type, num_coords)


@cython.cfunc
@cython.exceptval(check=False)
def isoform_position_details_free(
        ipd: p_IsoformPositionDetails) -> cython.void:
    free(ipd.coords)
    free(ipd.coord_types)


@cython.cfunc
@cython.exceptval(check=False)
def isoform_position_details_reset(
        ipd: p_IsoformPositionDetails) -> cython.void:
    ipd.coord_i = 0
    ipd.pos = ipd.coords[ipd.coord_i]
    ipd.coord_type = ipd.coord_types[ipd.coord_i]


@cython.cfunc
@cython.exceptval(check=False)
def isoform_position_details_is_at_start(
        ipd: p_IsoformPositionDetails) -> cython.bint:
    return ipd.coord_i == 0


@cython.cfunc
@cython.exceptval(check=False)
def isoform_position_details_is_at_end(
        ipd: p_IsoformPositionDetails) -> cython.bint:
    return ipd.coord_i == (ipd.num_coords - 1)


@cython.cfunc
@cython.exceptval(check=False)
def isoform_position_details_advance(
        ipd: p_IsoformPositionDetails) -> cython.void:
    ipd.coord_i += 1
    ipd.pos = ipd.coords[ipd.coord_i]
    ipd.coord_type = ipd.coord_types[ipd.coord_i]


# After this is run, at least one of the read or the isoform will
# still be at start.  If they have the same start then both will still
# be at start.  Otherwise the one not at start will be advanced until
# it hits or crosses the other start.
@cython.cfunc
@cython.exceptval(check=False)
def find_read_isoform_initial_overlap(
        read_pos_details: p_ReadPositionDetails,
        isoform_pos_details: p_IsoformPositionDetails) -> cython.bint:
    if read_pos_details.pos == isoform_pos_details.pos:
        return True

    if read_pos_details.pos < isoform_pos_details.pos:
        while not read_position_details_is_at_end(read_pos_details):
            read_position_details_advance(read_pos_details)
            if read_pos_details.pos >= isoform_pos_details.pos:
                return True

        return False

    while not isoform_position_details_is_at_end(isoform_pos_details):
        isoform_position_details_advance(isoform_pos_details)
        if isoform_pos_details.pos >= read_pos_details.pos:
            return True

    return False


@cython.cfunc
@cython.exceptval(check=False)
def check_start_boundary(
        read_pos_details: p_ReadPositionDetails,
        isoform_pos_details: p_IsoformPositionDetails) -> cython.bint:
    enum: CoordTypeEnumType
    enum = get_coord_type_enum()
    if isoform_position_details_is_at_start(isoform_pos_details):
        if read_pos_details.pos == read_pos_details.start:
            return True

        if read_position_details_is_at_exon_start(read_pos_details):
            # read spliced to 1st isoform exon start
            if read_pos_details.pos == isoform_pos_details.pos:
                return isoform_pos_details.coord_type != enum.isoform_start

            # read spliced into 1st isoform exon
            return isoform_pos_details.coord_type == enum.any_exon_start

        # read crosses isoform_start
        if isoform_pos_details.coord_type == enum.exon_start:
            return False

        if isoform_pos_details.coord_type == enum.any_exon_start:
            return True

        # isoform_pos_details.coord_type is one of:
        #   'isoform_start', 'exon_or_isoform_start'
        #
        # read 1st exon reads over isoform start
        return read_position_details_is_at_first_exon_end(read_pos_details)

    # read starts within isoform
    if isoform_coord_type_is_end(isoform_pos_details.coord_type):
        # read starts within an exon
        return True

    if isoform_pos_details.pos == read_pos_details.pos:
        # read starts at exon start
        return True

    # read starts in an intron
    return False


@cython.cfunc
@cython.exceptval(check=False)
def check_end_boundary(
        read_pos_details: p_ReadPositionDetails,
        isoform_pos_details: p_IsoformPositionDetails) -> cython.bint:
    enum: CoordTypeEnumType
    enum = get_coord_type_enum()
    if read_position_details_is_at_end(read_pos_details):
        if isoform_position_details_is_at_end(isoform_pos_details):
            if read_pos_details.pos <= isoform_pos_details.pos:
                return True

            # read crossed the isoform end
            return isoform_pos_details.coord_type != enum.exon_end

        # read ends within internal exon of isoform
        if read_pos_details.pos <= isoform_pos_details.pos:
            return True

        # read crossed internal exon boundary
        return False

    if isoform_position_details_is_at_end(isoform_pos_details):
        # isoform ends, but read continues to another exon
        if isoform_pos_details.coord_type == enum.any_exon_end:
            return True

        if isoform_pos_details.coord_type in [
                enum.exon_end, enum.exon_or_isoform_end
        ]:
            return read_pos_details.pos == isoform_pos_details.pos

        return False

    # Read and isoform should be at the same boundary type.
    # Previously this was checked, but now it's assumed.
    #
    # read_at_exon_start = read_position_details_is_at_exon_start(read_pos_details)
    # isoform_at_exon_start = isoform_coord_type_is_start(
    #     isoform_pos_details.coord_type)

    # different internal exon boundary
    return False


# A run with cProfile showed almost all time being spent in this function.
# Cython types are used to decrease the running time.
@cython.cfunc
@cython.exceptval(check=True)
def is_read_compatible_with_isoform(
        p_read_pos_details: p_ReadPositionDetails,
        p_isoform_pos_details: p_IsoformPositionDetails) -> cython.bint:
    had_overlap: cython.bint
    start_boundary_is_ok: cython.bint
    is_boundary_ok: cython.bint
    read_position_details_reset(p_read_pos_details)
    isoform_position_details_reset(p_isoform_pos_details)
    had_overlap = find_read_isoform_initial_overlap(p_read_pos_details,
                                                    p_isoform_pos_details)
    if not had_overlap:
        return False

    start_boundary_is_ok = check_start_boundary(p_read_pos_details,
                                                p_isoform_pos_details)
    if not start_boundary_is_ok:
        return False

    # The read and isoform coordinates have been advanced until the first overlap.
    # If only one is at an exon start then advance it to the exon end.
    if isoform_coord_type_is_start(p_isoform_pos_details.coord_type):
        if read_position_details_is_at_exon_end(p_read_pos_details):
            isoform_position_details_advance(p_isoform_pos_details)
    elif read_position_details_is_at_exon_start(p_read_pos_details):
        read_position_details_advance(p_read_pos_details)

    # The read and isoform are now both at an exon start or both at an exon end.
    # Advance along both the read and isoform while the coords match.
    while p_read_pos_details.pos == p_isoform_pos_details.pos:
        if ((read_position_details_is_at_end(p_read_pos_details)
             or isoform_position_details_is_at_end(p_isoform_pos_details))):
            break

        read_position_details_advance(p_read_pos_details)
        isoform_position_details_advance(p_isoform_pos_details)

    is_boundary_ok = check_end_boundary(p_read_pos_details,
                                        p_isoform_pos_details)
    return is_boundary_ok


# A read_id from the original alignment/fastq can appear in multiple genes.
# The original read_id will not be preserved. Instead it will have a different
# {read_i} for each gene.
@cython.cfunc
@cython.exceptval(check=True)
def count_read_for_asm(gene_id, read_i, sample,
                       p_read_pos_details: p_ReadPositionDetails, asm,
                       isoform_pos_details_for_asm: p_IsoformPositionDetails,
                       out_handle) -> cython.void:
    p_isoform_pos_details: p_IsoformPositionDetails
    asm_id = asm['asm_id']
    parsed = rmats_long_utils.parse_asm_id(asm_id)
    asm_i = parsed['event_i']
    isoform_ids = asm['isoform_ids']
    # The last column is the isoform ID
    out_columns = [str(asm_i), str(read_i), gene_id, sample, None]
    for isoform_i, isoform_id in enumerate(isoform_ids):
        p_isoform_pos_details = cython.address(
            isoform_pos_details_for_asm[isoform_i])
        is_match = is_read_compatible_with_isoform(p_read_pos_details,
                                                   p_isoform_pos_details)
        if is_match:
            isoform_id = isoform_ids[isoform_i]
            out_columns[4] = isoform_id
            rmats_long_utils.write_tsv_line(out_handle, out_columns)


@cython.cfunc
@cython.exceptval(check=True)
def create_isoform_objects_for_asms(asms: list) -> pp_IsoformPositionDetails:
    malloc_size: cython.size_t
    num_asms = len(asms)
    malloc_size = num_asms * cython.sizeof(p_IsoformPositionDetails)
    isoform_pos_details_for_asms = cython.cast(pp_IsoformPositionDetails,
                                               malloc(malloc_size))
    if isoform_pos_details_for_asms is cython.NULL:
        raise Exception(
            'create_isoform_objects_for_asms: null return from malloc')

    for asm_i, asm in enumerate(asms):
        isoforms = asm['isoforms']
        num_isoforms = len(isoforms)
        malloc_size = num_isoforms * cython.sizeof(IsoformPositionDetails)
        isoform_pos_details_for_asm = cython.cast(p_IsoformPositionDetails,
                                                  malloc(malloc_size))
        if isoform_pos_details_for_asm is cython.NULL:
            raise Exception(
                'create_isoform_objects_for_asms: null return from malloc')

        isoform_pos_details_for_asms[asm_i] = isoform_pos_details_for_asm
        for isoform_i, isoform in enumerate(isoforms):
            isoform_pos_details = isoform_position_details_init(isoform)
            isoform_pos_details_for_asm[isoform_i] = isoform_pos_details

    return isoform_pos_details_for_asms


@cython.cfunc
@cython.exceptval(check=True)
def free_isoform_objects_for_asms(
        asms: list, ipd_for_asms: pp_IsoformPositionDetails) -> cython.void:
    isoform_pos_details_for_asms: pp_IsoformPositionDetails
    p_isoform_pos_details: p_IsoformPositionDetails
    isoform_pos_details_for_asms = ipd_for_asms
    for asm_i, asm in enumerate(asms):
        isoforms = asm['isoforms']
        for isoform_i in range(len(isoforms)):
            p_isoform_pos_details = cython.address(
                isoform_pos_details_for_asms[asm_i][isoform_i])
            isoform_position_details_free(p_isoform_pos_details)

        free(isoform_pos_details_for_asms[asm_i])

    free(isoform_pos_details_for_asms)


def count_reads_for_gene_with_handles(gene_i, gene_id, asms, align_handle,
                                      out_handle):
    # IsoformPositionDetails[][]
    isoform_pos_details_for_asms: pp_IsoformPositionDetails
    isoform_pos_details_for_asm: p_IsoformPositionDetails
    # Create the Cython object for each isoform here.
    # Then each read can reuse the isoform objects.
    isoform_pos_details_for_asms = create_isoform_objects_for_asms(asms)

    read_i = 0
    line = align_handle.readline()
    while line:
        columns = rmats_long_utils.read_tsv_line(line)
        gene_i_str = columns[0]
        found_gene_i = int(gene_i_str)
        if found_gene_i > gene_i:
            break

        line = align_handle.readline()
        read_i += 1
        # gene_id = columns[1]
        sample = columns[2]
        start_str = columns[3]
        start = int(start_str)
        end_str = columns[4]
        end = int(end_str)
        sjs_str = columns[5]
        sjs = rmats_long_utils.parse_sjs_string(sjs_str)
        read_pos_details = read_position_details_init(start, end, sjs)
        p_read_pos_details = cython.address(read_pos_details)
        # strand = columns[6]
        for asm_i, asm in enumerate(asms):
            isoform_pos_details_for_asm = isoform_pos_details_for_asms[asm_i]
            count_read_for_asm(gene_id, read_i, sample, p_read_pos_details,
                               asm, isoform_pos_details_for_asm, out_handle)

        read_position_details_free(p_read_pos_details)

    free_isoform_objects_for_asms(asms, isoform_pos_details_for_asms)


def get_out_files_by_basename(thread_dirs):
    by_basename = dict()
    for thread_dir in thread_dirs:
        file_names = os.listdir(thread_dir)
        for file_name in file_names:
            path = os.path.join(thread_dir, file_name)
            for_basename = by_basename.get(file_name)
            if for_basename is None:
                for_basename = list()
                by_basename[file_name] = for_basename

            for_basename.append(path)

    return by_basename


def get_asm_i_from_out_line(line):
    columns = rmats_long_utils.read_tsv_line(line)
    asm_i_str = columns[0]
    asm_i = int(asm_i_str)
    return asm_i


def update_next_asm_i_and_line(next_asm_i_and_line, handle):
    line = handle.readline()
    if not line:
        asm_i = None
        line = None
    else:
        asm_i = get_asm_i_from_out_line(line)

    next_asm_i_and_line['asm_i'] = asm_i
    next_asm_i_and_line['line'] = line


def initialize_next_asm_i_and_line(paths, in_handles,
                                   next_asm_i_and_line_by_handle):
    initial_asm_i = None
    for path in paths:
        handle = open(path, 'rt')
        in_handles.append(handle)
        next_asm_i_and_line = dict()
        next_asm_i_and_line_by_handle.append(next_asm_i_and_line)
        update_next_asm_i_and_line(next_asm_i_and_line, handle)
        asm_i = next_asm_i_and_line['asm_i']
        if initial_asm_i is None:
            initial_asm_i = asm_i
        elif (asm_i is not None) and (asm_i < initial_asm_i):
            initial_asm_i = asm_i

    return initial_asm_i


def output_updated_tsv_line(chr_i, line, previous_gene_id, out_handle):
    in_columns = rmats_long_utils.read_tsv_line(line)
    asm_i = in_columns[0]
    read_i = in_columns[1]
    gene_id = in_columns[2]
    sample = in_columns[3]
    isoform = in_columns[4]
    write_gene_id = (gene_id != previous_gene_id)
    previous_gene_id = gene_id
    if not write_gene_id:
        gene_id = ''

    asm_id = '{}_{}'.format(chr_i, asm_i)
    out_columns = [asm_id, gene_id, read_i, sample, isoform]
    rmats_long_utils.write_tsv_line(out_handle, out_columns)
    return previous_gene_id


def combine_out_files_for_chr(chr_i, paths, out_handle):
    in_handles = list()
    next_asm_i_and_line_by_handle = list()
    previous_asm_i = None
    previous_gene_id = None
    try:
        initial_asm_i = initialize_next_asm_i_and_line(
            paths, in_handles, next_asm_i_and_line_by_handle)
        if initial_asm_i is None:
            return

        previous_asm_i = initial_asm_i - 1
        any_progress = True
        while any_progress:
            any_progress = False
            lowest_next_asm_i = None
            for handle_i, handle in enumerate(in_handles):
                next_asm_i_and_line = next_asm_i_and_line_by_handle[handle_i]
                while True:
                    asm_i = next_asm_i_and_line['asm_i']
                    line = next_asm_i_and_line['line']
                    # Keep writing for previous_asm_i, or start the next asm_i
                    if (asm_i is None) or (asm_i > (previous_asm_i + 1)):
                        break

                    previous_asm_i = asm_i
                    lowest_next_asm_i = asm_i
                    previous_gene_id = output_updated_tsv_line(
                        chr_i, line, previous_gene_id, out_handle)
                    update_next_asm_i_and_line(next_asm_i_and_line, handle)

                handle_asm_i = next_asm_i_and_line['asm_i']
                if lowest_next_asm_i is None:
                    lowest_next_asm_i = handle_asm_i
                elif (handle_asm_i is not None) and (handle_asm_i
                                                     < lowest_next_asm_i):
                    lowest_next_asm_i = handle_asm_i

            if lowest_next_asm_i is not None:
                previous_asm_i = lowest_next_asm_i - 1
                any_progress = True

    finally:
        for handle in in_handles:
            handle.close()


def combine_out_files_thread(thread_inputs):
    headers = ['asm_id', 'gene_id', 'read_id', 'sample_id', 'isoform_id']
    while True:
        arguments = rmats_long_utils.try_get_from_queue_with_short_wait(
            thread_inputs)
        if arguments is None:
            return

        chr_i = arguments['chr_i']
        paths = arguments['paths']
        out_path = arguments['out_path']
        with open(out_path, 'wt') as out_handle:
            rmats_long_utils.write_tsv_line(out_handle, headers)
            combine_out_files_for_chr(chr_i, paths, out_handle)


def combine_out_files(writer_thread_infos, num_threads, out_dir):
    thread_dirs = list()
    for thread_info in writer_thread_infos:
        thread_dir = thread_info['work_dir']
        thread_dirs.append(thread_dir)

    out_files_by_basename = get_out_files_by_basename(thread_dirs)
    num_jobs = len(out_files_by_basename)
    threads = list()
    thread_inputs = multiprocessing.Queue(num_jobs)
    for basename, paths in out_files_by_basename.items():
        chr_i = rmats_long_utils.get_chr_id_from_path(basename)
        out_path = os.path.join(out_dir, basename)
        arguments = {'chr_i': chr_i, 'paths': paths, 'out_path': out_path}
        thread_inputs.put(arguments)

    for _ in range(num_threads):
        thread = multiprocessing.Process(target=combine_out_files_thread,
                                         args=(thread_inputs, ))
        threads.append(thread)
        thread.start()

    rmats_long_utils.join_threads_and_raise_if_error(threads)


def read_events_with_handles(align_path, next_queue_i, gtf_handle,
                             event_handle, queues_for_reader):
    num_queues = len(queues_for_reader)
    last_queue_i = num_queues - 1
    gene_i = 0
    event_handle_and_line = rmats_long_utils.HandleAndNextLine(event_handle)
    while True:
        this_gene_i = gene_i
        gene_i += 1
        gene_id = read_next_gene_id(gtf_handle)
        if gene_id is None:
            break

        asms = load_asms_for_gene_id(gene_id, event_handle_and_line)
        if not asms:
            continue

        arguments = {
            'asms': asms,
            'gene_id': gene_id,
            'gene_i': this_gene_i,
            'align_path': align_path,
        }
        while True:
            queue = queues_for_reader[next_queue_i]
            if next_queue_i == last_queue_i:
                next_queue_i = 0
            else:
                next_queue_i += 1

            was_put = rmats_long_utils.try_put_to_queue_with_short_wait(
                queue, arguments)
            if was_put:
                break

    return next_queue_i


def read_events_thread(reader_input_queue, queues_for_reader):
    next_queue_i = 0
    while True:
        arguments = reader_input_queue.get()
        if len(arguments) == 0:
            # pass signal to writer threads as empty dict()
            for writer_queue in queues_for_reader:
                writer_queue.put(dict())

            return

        gtf_path = arguments['gtf_path']
        event_path = arguments['event_path']
        align_path = arguments['align_path']
        with open(gtf_path, 'rt') as gtf_handle:
            with open(event_path, 'rt') as event_handle:
                next_queue_i = read_events_with_handles(
                    align_path, next_queue_i, gtf_handle, event_handle,
                    queues_for_reader)


def load_offsets_by_gene_from_align_index(align_index_path):
    offsets_by_gene = dict()
    with open(align_index_path, 'rt') as handle:
        for line in handle:
            columns = rmats_long_utils.read_tsv_line(line)
            gene_id = columns[0]
            offset = columns[1]
            offset = int(offset)
            offsets_by_gene[gene_id] = offset

    return offsets_by_gene


def get_offset_for_gene(gene_id, align_index_path, align_index_by_path):
    offsets_by_gene = align_index_by_path.get(align_index_path)
    if offsets_by_gene is None:
        offsets_by_gene = load_offsets_by_gene_from_align_index(
            align_index_path)
        align_index_by_path[align_index_path] = offsets_by_gene

    offset = offsets_by_gene.get(gene_id)
    return offset


def count_reads_for_gene(gene_i, gene_id, asms, align_path, align_index_path,
                         align_index_by_path, out_path, out_paths):
    offset = get_offset_for_gene(gene_id, align_index_path,
                                 align_index_by_path)
    if offset is None:
        return

    out_paths.add(out_path)
    with open(align_path, 'rt') as align_handle:
        align_handle.seek(offset)
        with open(out_path, 'at') as out_handle:
            count_reads_for_gene_with_handles(gene_i, gene_id, asms,
                                              align_handle, out_handle)


def sort_file_by_asm_and_read(sort_buffer_size, out_path):
    tmp_dir = os.path.dirname(out_path)
    tmp_path = '{}.tmp'.format(out_path)
    asm_key_arg = '-k1,1g'
    read_key_arg = '-k2,2g'
    env = {'LC_ALL': 'C'}  # to ensure sort order
    command = [
        'sort', '--buffer-size', sort_buffer_size, '--temporary-directory',
        tmp_dir, asm_key_arg, read_key_arg, '--output', tmp_path, out_path
    ]
    subprocess.run(command, env=env, check=True)

    shutil.move(tmp_path, out_path)


def process_and_write_results_thread_main(writer_queues, work_dir):
    expected_signals = len(writer_queues)
    out_paths = set()
    align_index_by_path = dict()
    while True:
        all_empty = True
        for queue in writer_queues:
            arguments = rmats_long_utils.try_get_from_queue_without_wait(queue)
            if arguments is None:
                continue

            if len(arguments) == 0:
                expected_signals -= 1
                if expected_signals == 0:
                    return out_paths

                continue

            all_empty = False
            asms = arguments['asms']
            gene_id = arguments['gene_id']
            gene_i = arguments['gene_i']
            align_path = arguments['align_path']
            align_index_path = '{}.index'.format(align_path)
            basename = os.path.basename(align_path)
            out_path = os.path.join(work_dir, basename)
            count_reads_for_gene(gene_i, gene_id, asms, align_path,
                                 align_index_path, align_index_by_path,
                                 out_path, out_paths)

        if all_empty:
            time.sleep(1)  # 1 second


def process_and_write_results_thread(sort_buffer_size, writer_queues,
                                     work_dir):
    os.makedirs(work_dir)
    out_paths = process_and_write_results_thread_main(writer_queues, work_dir)
    for out_path in out_paths:
        sort_file_by_asm_and_read(sort_buffer_size, out_path)


def get_file_info_for_chrs(event_dir, align_dir, gtf_dir):
    chr_file_infos = list()
    gtf_file_names = os.listdir(gtf_dir)
    for file_name in gtf_file_names:
        if not file_name.startswith('chr_id_'):
            continue

        gtf_path = os.path.join(gtf_dir, file_name)
        event_path = os.path.join(event_dir, file_name)
        align_path = os.path.join(align_dir, file_name)
        # Some chrs may not have events or alignments
        if not (os.path.exists(event_path) and os.path.exists(align_path)):
            continue

        size_in_bytes = rmats_long_utils.get_file_size(event_path)
        chr_info = {
            'gtf_path': gtf_path,
            'event_path': event_path,
            'align_path': align_path,
            'event_size': size_in_bytes,
        }
        chr_file_infos.append(chr_info)

    chr_file_infos.sort(key=lambda x: x['event_size'], reverse=True)
    return chr_file_infos


def create_reader_to_writer_queues(num_threads):
    queue_size = 1
    queues_by_reader = list()
    for _ in range(num_threads):
        queues_for_reader = list()
        queues_by_reader.append(queues_for_reader)
        for _ in range(num_threads):
            queue = multiprocessing.Queue(queue_size)
            queues_for_reader.append(queue)

    return queues_by_reader


def create_writer_threads(num_threads, sort_buffer_size,
                          reader_to_writer_queues, temp_dir):
    thread_infos = list()
    for writer_i in range(num_threads):
        writer_queues = list()
        for queues_for_reader in reader_to_writer_queues:
            writer_queue = queues_for_reader[writer_i]
            writer_queues.append(writer_queue)

        work_dir = os.path.join(temp_dir, str(writer_i))
        thread = multiprocessing.Process(
            target=process_and_write_results_thread,
            args=(sort_buffer_size, writer_queues, work_dir))
        info = {
            'thread': thread,
            'work_dir': work_dir,
        }
        thread_infos.append(info)
        thread.start()

    return thread_infos


def create_reader_threads(reader_input_queue, reader_to_writer_queues):
    threads = list()
    for queues_for_reader in reader_to_writer_queues:
        thread = multiprocessing.Process(target=read_events_thread,
                                         args=(reader_input_queue,
                                               queues_for_reader))
        threads.append(thread)
        thread.start()

    return threads


def cleanup_threads(reader_threads, writer_thread_infos, reader_input_queue,
                    reader_to_writer_queues):
    # signal threads to stop
    for _ in reader_threads:
        reader_input_queue.put(dict())

    all_threads = list()
    all_threads.extend(reader_threads)
    for thread_info in writer_thread_infos:
        all_threads.append(thread_info['thread'])

    # attempt to join all threads
    failed_threads = list()
    while all_threads and (not failed_threads):
        running_threads = list()
        for thread in all_threads:
            join_result = rmats_long_utils.try_join_thread(thread)
            if join_result == 'running':
                running_threads.append(thread)
                continue

            if join_result == 'joined':
                continue

            # join_result == 'error'
            failed_threads.append(thread)

        all_threads = running_threads

    rmats_long_utils.drain_queue(reader_input_queue)
    for queues_for_reader in reader_to_writer_queues:
        for queue in queues_for_reader:
            rmats_long_utils.drain_queue(queue)

    # Force threads to stop if there was an error
    messages = list()
    if all_threads:
        messages.append('{} thread(s) were not joined'.format(
            len(all_threads)))
        for thread in all_threads:
            thread.terminate()

    if failed_threads:
        messages.append('{} thread(s) had non-zero exit codes'.format(
            len(failed_threads)))

    if messages:
        raise Exception('. '.join(messages))


def count_reads_for_asms(event_dir, align_dir, gtf_dir, num_threads,
                         sort_buffer_size, out_dir):
    with tempfile.TemporaryDirectory(suffix='_tmp',
                                     prefix='count_reads_for_asms',
                                     dir=out_dir) as temp_dir:
        chr_file_infos = get_file_info_for_chrs(event_dir, align_dir, gtf_dir)
        # Enough space for each file and the stop signals
        reader_queue_size = len(chr_file_infos) + num_threads
        reader_input_queue = multiprocessing.Queue(reader_queue_size)
        # each writer thread gets an input queue from each reader thread
        reader_to_writer_queues = create_reader_to_writer_queues(num_threads)
        writer_thread_infos = create_writer_threads(num_threads,
                                                    sort_buffer_size,
                                                    reader_to_writer_queues,
                                                    temp_dir)
        reader_threads = create_reader_threads(reader_input_queue,
                                               reader_to_writer_queues)
        for chr_file_info in chr_file_infos:
            gtf_path = chr_file_info['gtf_path']
            event_path = chr_file_info['event_path']
            align_path = chr_file_info['align_path']
            arguments = {
                'gtf_path': gtf_path,
                'event_path': event_path,
                'align_path': align_path,
            }
            reader_input_queue.put(arguments)

        cleanup_threads(reader_threads, writer_thread_infos,
                        reader_input_queue, reader_to_writer_queues)

        combine_out_files(writer_thread_infos, num_threads, out_dir)


def main():
    args = parse_args()
    gtf_dir = os.path.abspath(args.gtf_dir)
    align_dir = os.path.abspath(args.align_dir)
    event_dir = os.path.abspath(args.event_dir)
    out_dir = os.path.abspath(args.out_dir)
    rmats_long_utils.create_output_dir(out_dir)
    count_reads_for_asms(event_dir, align_dir, gtf_dir, args.num_threads,
                         args.sort_buffer_size, out_dir)
    print('count_reads_for_asms.py finished')


if __name__ == '__main__':
    main()
