import os
import os.path
import shutil
import subprocess
import sys
import time
import unittest


def write_tsv_line(handle, columns):
    handle.write('{}\n'.format('\t'.join(columns)))


def read_tsv_line(line):
    return line.rstrip('\n').split('\t')


def remove_then_create_directories(dirs):
    for directory in dirs:
        if os.path.exists(directory):
            shutil.rmtree(directory)

        os.makedirs(directory)


def format_gtf_attributes(key_value_pairs):
    pair_strings = list()
    for key, value in key_value_pairs:
        string = '{} "{}";'.format(key, value)
        pair_strings.append(string)

    return ' '.join(pair_strings)


def write_gtf_gene_line(handle, chr_name, gene_id, start, end, strand):
    attributes = format_gtf_attributes([('gene_id', gene_id)])
    columns = [
        chr_name, 'base_test', 'gene',
        str(start),
        str(end), '.', strand, '.', attributes
    ]
    write_tsv_line(handle, columns)


def write_gtf_transcript_line(handle, chr_name, gene_id, transcript_id, start,
                              end, strand):
    attributes = format_gtf_attributes([('gene_id', gene_id),
                                        ('transcript_id', transcript_id)])
    columns = [
        chr_name, 'base_test', 'transcript',
        str(start),
        str(end), '.', strand, '.', attributes
    ]
    write_tsv_line(handle, columns)


def write_gtf_exon_line(handle, chr_name, gene_id, transcript_id, start, end,
                        strand):
    # exon_number is not needed
    attributes = format_gtf_attributes([('gene_id', gene_id),
                                        ('transcript_id', transcript_id)])
    columns = [
        chr_name, 'base_test', 'exon',
        str(start),
        str(end), '.', strand, '.', attributes
    ]
    write_tsv_line(handle, columns)


def write_gtf(path, chroms):
    with open(path, 'wt') as handle:
        for chrom in chroms:
            for gene in chrom.genes:
                if not gene.isoforms:
                    continue

                gene_info = gene.get_start_end_and_strand_from_isoforms()
                gene_start = gene_info[0] + 1
                gene_end = gene_info[1] + 1
                gene_strand = gene_info[2]
                write_gtf_gene_line(handle, chrom.name, gene.id, gene_start,
                                    gene_end, gene_strand)
                for isoform in gene.isoforms:
                    isoform_info = isoform.get_start_and_end_from_exons()
                    isoform_start = isoform_info[0] + 1
                    isoform_end = isoform_info[1] + 1
                    write_gtf_transcript_line(handle, chrom.name, gene.id,
                                              isoform.id, isoform_start,
                                              isoform_end, isoform.strand)
                    for exon in isoform.exons:
                        write_gtf_exon_line(handle, chrom.name, gene.id,
                                            isoform.id, exon.start + 1,
                                            exon.end + 1, isoform.strand)


def write_sam_chr_headers(handle, chroms):
    for chrom in chroms:
        length = 0
        for gene in chrom.genes:
            gene_info = gene.get_start_end_and_strand_from_isoforms()
            gene_end = gene_info[1]
            if gene_end > length:
                length = gene_end

        if length == 0:
            continue

        header_columns = [
            '@SQ', 'SN:{}'.format(chrom.name), 'LN:{}'.format(length + 1)
        ]
        write_tsv_line(handle, header_columns)


def write_sam_alignment_line(handle, align, next_i):
    if align.read_id is None:
        read_id = 'read-{}'.format(next_i)
    else:
        read_id = align.read_id

    if align.sequence is None:
        seq = '*'
    else:
        seq = align.sequence

    if align.sequence is None:
        qual = '*'
    else:
        qual = align.sequence

    ref_next = '*'
    pos_next = '0'
    template_length = '0'
    columns = [
        read_id,
        str(align.flag()), align.chr_name,
        str(align.start),
        str(align.mapq),
        str(align.cigar), ref_next, pos_next, template_length, seq, qual
    ]
    write_tsv_line(handle, columns)


def write_sam(path, chroms, aligns, next_i):
    with open(path, 'wt') as handle:
        write_sam_chr_headers(handle, chroms)
        for align in aligns:
            write_sam_alignment_line(handle, align, next_i)
            next_i += 1


def sam_flag_int(has_multiple_segments=False,
                 all_segments_aligned=True,
                 segment_unmapped=False,
                 next_segment_unmapped=False,
                 is_reversed=False,
                 next_segment_is_reversed=False,
                 is_first_segment=True,
                 is_last_segment=False,
                 is_secondary=False,
                 failed_filter=False,
                 is_duplicate=False,
                 is_supplementary=False):
    flags_int = 0
    if has_multiple_segments:
        flags_int += 1
    if all_segments_aligned:
        flags_int += 2
    if segment_unmapped:
        flags_int += 4
    if next_segment_unmapped:
        flags_int += 8
    if is_reversed:
        flags_int += 16
    if next_segment_is_reversed:
        flags_int += 32
    if is_first_segment:
        flags_int += 64
    if is_last_segment:
        flags_int += 128
    if is_secondary:
        flags_int += 256
    if failed_filter:
        flags_int += 512
    if is_duplicate:
        flags_int += 1024
    if is_supplementary:
        flags_int += 2048

    return flags_int


def append_copies(to_copy, num_copies, values):
    for _ in range(num_copies):
        copy = to_copy.copy()
        values.append(copy)


def run_command(command, log=None, out=None, check=True):
    if (log is not None) and (out is not None):
        with open(log, 'wb') as log_handle:
            with open(out, 'wb') as out_handle:
                process = subprocess.run(command,
                                         check=check,
                                         stdout=out_handle,
                                         stderr=log_handle)
    elif out is None:
        with open(log, 'wb') as log_handle:
            process = subprocess.run(command,
                                     check=check,
                                     stdout=log_handle,
                                     stderr=subprocess.STDOUT)
    elif log is None:
        with open(out, 'wb') as out_handle:
            process = subprocess.run(command, check=check, stdout=out_handle)
    else:
        process = subprocess.run(command, check=check)

    return process


def get_python_executable():
    exe = sys.executable
    if not exe:
        return 'python'

    return exe


def create_group_1_file(sample_names, dir_path):
    return create_group_file(sample_names, dir_path, '1')


def create_group_2_file(sample_names, dir_path):
    return create_group_file(sample_names, dir_path, '2')


def create_group_file(sample_names, dir_path, group_i):
    group_name = 'group_{}.txt'.format(group_i)
    group_path = os.path.join(dir_path, group_name)
    with open(group_path, 'wt') as handle:
        handle.write('{}\n'.format(','.join(sample_names)))

    return group_path


def start_timer():
    start_time = time.time()
    return start_time


class BaseTest(unittest.TestCase):
    def setUp(self):
        self._python_exe = get_python_executable()
        self._base_test_dir = os.path.dirname(__file__)
        self._base_parent_dir = os.path.dirname(self._base_test_dir)
        self._script_dir = os.path.join(self._base_parent_dir, 'scripts')
        self._rmats_long_py = os.path.join(self._script_dir, 'rmats_long.py')
        self._clean_espresso_gtf_py = os.path.join(self._script_dir,
                                                   'clean_espresso_gtf.py')
        self._organize_gene_info_by_chr_py = os.path.join(
            self._script_dir, 'organize_gene_info_by_chr.py')
        self._organize_alignment_info_by_gene_and_chr_py = os.path.join(
            self._script_dir, 'organize_alignment_info_by_gene_and_chr.py')
        self._simplify_alignment_info_py = os.path.join(
            self._script_dir, 'simplify_alignment_info.py')
        self._detect_splicing_events_py = os.path.join(
            self._script_dir, 'detect_splicing_events.py')
        self._count_reads_for_asms_py = os.path.join(
            self._script_dir, 'count_reads_for_asms.py')

    def assert_exists(self, path):
        if os.path.exists(path):
            return

        self.fail('{} does not exist'.format(path))

    def assert_seconds_at_most(self, start_time, seconds):
        end_time = time.time()
        total_seconds = end_time - start_time
        self.assertLess(total_seconds, seconds)

    def run_clean_gtf(self, in_gtf, out_gtf, log_dir):
        command = [
            self._python_exe, self._clean_espresso_gtf_py, '--in-gtf', in_gtf,
            '--out-gtf', out_gtf
        ]
        log_path = os.path.join(log_dir, 'clean_gtf.log')
        run_command(command, log=log_path)

    def run_organize_gtf(self, gtf, out_dir, log_dir):
        gtf_dir = os.path.join(out_dir, 'org_gtf')
        command = [
            self._python_exe, self._organize_gene_info_by_chr_py, '--gtf', gtf,
            '--out-dir', gtf_dir
        ]
        log_path = os.path.join(log_dir, 'organize_gtf.log')
        run_command(command, log=log_path)
        return gtf_dir

    def run_simplify_alignments(self, sams_by_sample, input_dir, log_dir):
        simplified_by_sample = dict()
        for sample, sams in sams_by_sample.items():
            for_sample = list()
            simplified_by_sample[sample] = for_sample
            for sam_i, sam in enumerate(sams):
                out_path = os.path.join(
                    input_dir, '{}_{}_simplified.tsv'.format(sample, sam_i))
                for_sample.append(out_path)
                command = [
                    self._python_exe, self._simplify_alignment_info_py,
                    '--in-file', sam, '--out-tsv', out_path
                ]
                log_path = os.path.join(
                    log_dir,
                    'simplify_alignments_{}_{}.log'.format(sample, sam_i))
                run_command(command, log=log_path)

        return simplified_by_sample

    def run_organize_alignments(self, gtf_dir, simplified_by_sample, input_dir,
                                out_dir, log_dir):
        align_dir = os.path.join(out_dir, 'org_align')
        samples_tsv = os.path.join(input_dir, 'samples.tsv')
        with open(samples_tsv, 'wt') as handle:
            for sample, tsvs in simplified_by_sample.items():
                for tsv in tsvs:
                    write_tsv_line(handle, [sample, tsv])

        command = [
            self._python_exe, self._organize_alignment_info_by_gene_and_chr_py,
            '--gtf-dir', gtf_dir, '--out-dir', align_dir, '--samples-tsv',
            samples_tsv
        ]
        log_path = os.path.join(log_dir, 'organize_alignments.log')
        run_command(command, log=log_path)
        return align_dir

    def run_detect_events(self,
                          gtf_dir,
                          out_dir,
                          log_dir,
                          align_dir=None,
                          min_reads=None,
                          max_nodes=None,
                          max_paths=None,
                          output_full_gene_asm=False,
                          simplify_gene_isoform_endpoints=False,
                          filter_gene_isoforms_by_edge=False):
        event_dir = os.path.join(out_dir, 'events')
        command = [
            self._python_exe, self._detect_splicing_events_py, '--gtf-dir',
            gtf_dir, '--out-dir', event_dir
        ]

        if align_dir:
            command.extend(['--align-dir', align_dir])

        if min_reads is not None:
            command.extend(['--min-reads-per-edge', str(min_reads)])

        if max_nodes is not None:
            command.extend(['--max-nodes-in-event', str(max_nodes)])

        if max_paths is not None:
            command.extend(['--max-paths-in-event', str(max_paths)])

        if output_full_gene_asm:
            command.append('--output-full-gene-asm')

        if simplify_gene_isoform_endpoints:
            command.append('--simplify-gene-isoform-endpoints')

        if filter_gene_isoforms_by_edge:
            command.append('--filter-gene-isoforms-by-edge')

        log_path = os.path.join(log_dir, 'detect_events.log')
        run_command(command, log=log_path)
        return event_dir

    def run_count_reads(self, gtf_dir, align_dir, event_dir, out_dir, log_dir):
        asm_counts = os.path.join(out_dir, 'read_counts.tsv')
        command = [
            self._python_exe, self._count_reads_for_asms_py, '--gtf-dir',
            gtf_dir, '--event-dir', event_dir, '--align-dir', align_dir,
            '--out-tsv', asm_counts
        ]
        log_path = os.path.join(log_dir, 'count_reads.log')
        run_command(command, log=log_path)
        return asm_counts

    def run_rmats_long(self,
                       group_1,
                       group_2,
                       align_dir,
                       event_dir,
                       asm_counts,
                       out_dir,
                       log_dir,
                       gtf=None,
                       average_cpm_per_group=None,
                       min_cpm_per_group=None,
                       asm_proportion_gene=None,
                       average_reads_per_group=None):
        rmats_long_out = os.path.join(out_dir, 'rmats_long')
        command = [
            self._python_exe, self._rmats_long_py, '--asm-counts', asm_counts,
            '--event-dir', event_dir, '--align-dir', align_dir, '--group-1',
            group_1, '--group-2', group_2, '--out-dir', rmats_long_out
        ]
        if gtf:
            command.extend(['--gencode-gtf', gtf])

        if average_cpm_per_group is not None:
            command.extend(
                ['--average-cpm-per-group',
                 str(average_cpm_per_group)])

        if min_cpm_per_group is not None:
            command.extend(['--min-cpm-per-group', str(min_cpm_per_group)])

        if asm_proportion_gene is not None:
            command.extend(
                ['--asm-proportion-of-gene',
                 str(asm_proportion_gene)])

        if average_reads_per_group is not None:
            command.extend(
                ['--average-reads-per-group',
                 str(average_reads_per_group)])

        log_path = os.path.join(log_dir, 'rmats_long.log')
        run_command(command, log=log_path)
        return rmats_long_out

    def parse_read_compatibility(self, path):
        counts = dict()
        expected_headers = [
            'asm_id', 'gene_id', 'read_id', 'sample_id', 'isoform_id'
        ]
        gene_id = None
        with open(path, 'rt') as handle:
            for line_i, line in enumerate(handle):
                columns = read_tsv_line(line)
                if line_i == 0:
                    headers = columns
                    self.assertEqual(headers, expected_headers)
                    continue

                row = dict(zip(headers, columns))
                row_gene = row['gene_id']
                if (gene_id is None) or (row_gene not in ['', gene_id]):
                    gene_id = row_gene

                asm = row['asm_id']
                read = row['read_id']
                sample = row['sample_id']
                isoform = row['isoform_id']
                by_sample = counts.get(isoform)
                if by_sample is None:
                    by_sample = dict()
                    counts[isoform] = by_sample

                old_count = by_sample.get(sample, 0)
                by_sample[sample] = old_count + 1

        return counts

    def parse_asm_counts(self, path):
        counts = dict()
        with open(path, 'rt') as handle:
            for line_i, line in enumerate(handle):
                columns = read_tsv_line(line)
                if line_i == 0:
                    headers = columns
                    expected_headers = [
                        'asm_id', 'gene_id', 'sample_id', 'isoform_id',
                        'count', 'prop'
                    ]
                    self.assertEqual(headers, expected_headers)
                    continue

                row = dict(zip(headers, columns))
                sample = row['sample_id']
                isoform = row['isoform_id']
                count_str = row['count']
                by_sample = counts.get(isoform)
                if by_sample is None:
                    by_sample = dict()
                    counts[isoform] = by_sample

                by_sample[sample] = float(count_str)

        return counts

    def parse_summary_txt(self, path):
        summary = dict()
        with open(path, 'rt') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue

                parts = line.split(':')
                self.assertEqual(len(parts), 2)
                key, value = parts
                key = key.strip()
                value = value.strip()
                summary[key] = value

        return summary

    def parse_differential_transcripts(self, filt_path, is_asm=False):
        by_transcript = dict()
        if is_asm:
            feature_header = 'isoform_id'
            expected_initial_headers = [
                'asm_id', 'gene_id', 'isoform_id', 'lr', 'df', 'pvalue',
                'adj_pvalue'
            ]
        else:
            feature_header = 'feature_id'
            expected_initial_headers = [
                'gene_id', 'feature_id', 'lr', 'df', 'pvalue', 'adj_pvalue'
            ]

        with open(filt_path, 'rt') as handle:
            for line_i, line in enumerate(handle):
                columns = read_tsv_line(line)
                if line_i == 0:
                    len_initial = len(expected_initial_headers)
                    self.assertEqual(columns[:len_initial],
                                     expected_initial_headers)
                    next_non_sample_header = 'group_1_average_proportion'
                    samples = list()
                    prop_suffix = '_proportion'
                    for header in columns[len_initial:]:
                        if header == next_non_sample_header:
                            break

                        if not header.endswith(prop_suffix):
                            self.fail('Unexpected header: {}'.format(header))

                        sample = header[:-len(prop_suffix)]
                        samples.append(sample)

                    next_header_i = len_initial + len(samples)
                    expected_remaining_headers = [
                        'group_1_average_proportion',
                        'group_2_average_proportion',
                        'delta_isoform_proportion'
                    ]
                    for sample in samples:
                        expected_remaining_headers.append(
                            '{}_count'.format(sample))
                    for sample in samples:
                        expected_remaining_headers.append(
                            '{}_cpm'.format(sample))

                    self.assertEqual(columns[next_header_i:],
                                     expected_remaining_headers)
                    headers = columns
                    continue

                row = dict(zip(headers, columns))
                transcript = row[feature_header]
                by_transcript[transcript] = row

        return by_transcript

    def parse_id_to_chr(self, path):
        id_to_chr = dict()
        expected_headers = ['name', 'id']
        mapping_path = os.path.join(path, 'chr_name_id_mapping.tsv')
        with open(mapping_path, 'rt') as handle:
            for line_i, line in enumerate(handle):
                columns = read_tsv_line(line)
                if line_i == 0:
                    headers = columns
                    self.assertEqual(headers, expected_headers)
                    continue

                row = dict(zip(headers, columns))
                id_str = row['id']
                chr_name = row['name']
                id_to_chr[id_str] = chr_name

        return id_to_chr

    def parse_asm_definitions(self, path):
        asms = dict()
        expected_headers = [
            'asm_id', 'gene_i', 'gene_id', 'event_type', 'strand', 'start',
            'end', 'isoforms', 'isoform_ids', 'start_always_ss',
            'start_never_ss', 'end_always_ss', 'end_never_ss', 'is_strict'
        ]
        id_to_chr = self.parse_id_to_chr(path)
        for id_str, chr_name in id_to_chr.items():
            file_name = 'chr_id_{}.tsv'.format(id_str)
            asm_path = os.path.join(path, file_name)
            with open(asm_path, 'rt') as handle:
                prev_gene_i = None
                prev_gene_id = None
                for line_i, line in enumerate(handle):
                    columns = read_tsv_line(line)
                    if line_i == 0:
                        headers = columns
                        self.assertEqual(headers, expected_headers)
                        continue

                    row = dict(zip(headers, columns))
                    # gene_id is only set for the 1st line for each gene_i
                    if (prev_gene_i is None) or (row['gene_i'] != prev_gene_i):
                        prev_gene_i = row['gene_i']
                        prev_gene_id = row['gene_id']
                    else:
                        row['gene_id'] = prev_gene_id

                    row['chr'] = chr_name
                    isoform_ids_str = row['isoform_ids']
                    row['isoform_ids'] = isoform_ids_str.split(';')
                    isoforms_str = row['isoforms']
                    isoform_strs = isoforms_str.split(';')
                    isoforms_exons = list()
                    row['isoforms'] = isoforms_exons
                    for string in isoform_strs:
                        exons = list()
                        coord_strs = string.split(',')
                        num_exons, extra = divmod(len(coord_strs), 2)
                        self.assertEqual(extra, 0, msg=coord_strs)
                        for exon_i in range(num_exons):
                            start_i = exon_i * 2
                            start_str = coord_strs[start_i]
                            end_str = coord_strs[start_i + 1]
                            start = int(start_str)
                            end = int(end_str)
                            exons.append((start, end))

                        isoforms_exons.append(tuple(exons))

                    asm_key = tuple(isoforms_exons)
                    self.assertNotIn(asm_key, asms)
                    asms[asm_key] = row

        return asms

    def parse_isoform_diff_tsv(self, path):
        isoform_diffs = list()
        expected_headers = [
            'transcript1', 'transcript2', 'event', 'coordinates'
        ]
        with open(path, 'rt') as handle:
            for line_i, line in enumerate(handle):
                columns = read_tsv_line(line)
                if line_i == 0:
                    headers = columns
                    self.assertEqual(headers, expected_headers)
                    continue

                row = dict(zip(headers, columns))
                isoform_diffs.append(row)

        return isoform_diffs


class Cigar:
    def __init__(self):
        self.operations = list()

    def __str__(self):
        parts = list()
        for op in self.operations:
            parts.append(str(op))

        return ''.join(parts)

    def copy(self):
        copy = Cigar()
        for op in self.operations:
            copy.operations.append(op.copy())

        return copy


class CigarOp:
    def __init__(self, num, char):
        self.num = num
        self.char = char

    def __str__(self):
        return '{}{}'.format(self.num, self.char)

    def copy(self):
        copy = CigarOp(self.num, self.char)
        return copy


class Alignment:
    def __init__(self):
        self.chr_name = None
        self.start = None  # 1-based
        self.cigar = Cigar()
        self.mapq = 30  # default to 30 which is 0.999 accuracy
        # even for - strand the sequence written to the SAM is the + strand
        self.sequence = None
        self.quality = None
        self.strand = '+'
        self.read_id = None

    def flag(self):
        is_minus_strand = self.strand == '-'
        return sam_flag_int(is_reversed=is_minus_strand)

    def match_to_isoform(self, isoform):
        self.strand = isoform.strand
        self.start = isoform.exons[0].start + 1
        ops = list()
        self.cigar.operations = ops
        prev_end = None
        for exon in isoform.exons:
            if prev_end is not None:
                gap = (exon.start - prev_end) - 1
                ops.append(CigarOp(gap, 'N'))

            prev_end = exon.end
            length = exon.length()
            ops.append(CigarOp(length, 'M'))

    def copy(self):
        copy = Alignment()
        copy.chr_name = self.chr_name
        copy.start = self.start
        copy.cigar = self.cigar.copy()
        copy.mapq = self.mapq
        copy.sequence = self.sequence
        copy.quality = self.quality
        copy.strand = self.strand
        copy.read_id = self.read_id
        return copy


class Region:
    def __init__(self, start=None, end=None):
        self.start = start  # 0-based
        self.end = end  # 0-based

    def length(self):
        return (self.end - self.start) + 1

    def copy(self):
        return Region(self.start, self.end)


class Isoform:
    def __init__(self):
        self.id = None
        self.name = None
        self.strand = '+'
        self.exons = list()

    def get_start_and_end_from_exons(self):
        isoform_start = None
        isoform_end = None
        for exon in self.exons:
            if (isoform_start is None) or (exon.start < isoform_start):
                isoform_start = exon.start

            if (isoform_end is None) or (exon.end > isoform_end):
                isoform_end = exon.end

        return (isoform_start, isoform_end)

    def copy(self):
        copy = Isoform()
        copy.id = self.id
        copy.name = self.name
        copy.strand = self.strand
        copy.exons = [exon.copy() for exon in self.exons]
        return copy


class Gene:
    def __init__(self):
        self.id = None
        self.name = None
        self.isoforms = list()

    def get_start_end_and_strand_from_isoforms(self):
        gene_start = None
        gene_end = None
        gene_strand = None
        for isoform in self.isoforms:
            if gene_strand is None:
                gene_strand = isoform.strand

            for exon in isoform.exons:
                if (gene_start is None) or (exon.start < gene_start):
                    gene_start = exon.start

                if (gene_end is None) or (exon.end > gene_end):
                    gene_end = exon.end

        return (gene_start, gene_end, gene_strand)


class Chromosome:
    def __init__(self):
        self.name = None
        self.genes = list()
