import os
import os.path
import unittest

import tests.base_test


class BasicEventsBaseTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_dir = os.path.dirname(__file__)
        self._input_dir = None
        self._out_dir = None
        self._log_dir = None
        self._gtf = None
        self._sams_by_sample = None

    def _set_input_out_log_with_prefix(self, prefix):
        self._input_dir = os.path.join(self._test_dir,
                                       '{}_input'.format(prefix))
        self._out_dir = os.path.join(self._test_dir, '{}_out'.format(prefix))
        self._log_dir = os.path.join(self._test_dir, '{}_log'.format(prefix))

    def _initialize_dirs(self):
        dirs = [self._input_dir, self._out_dir, self._log_dir]
        tests.base_test.remove_then_create_directories(dirs)

    def _test(self):
        self._initialize_dirs()
        self._create_test_inputs()
        self._gtf_dir = self.run_organize_gtf(self._gtf, self._out_dir,
                                              self._log_dir)
        self._simplified_by_sample = self.run_simplify_alignments(
            self._sams_by_sample, self._input_dir, self._log_dir)
        self._align_dir = self.run_organize_alignments(
            self._gtf_dir, self._simplified_by_sample, self._input_dir,
            self._out_dir, self._log_dir)
        self._event_dir = self.run_detect_events(self._gtf_dir,
                                                 self._out_dir,
                                                 self._log_dir,
                                                 output_basic_events=True)
        self._asm_counts = self.run_count_reads(self._gtf_dir, self._align_dir,
                                                self._event_dir, self._out_dir,
                                                self._log_dir)
        self._rmats_long_out = self.run_rmats_long(self._group_1,
                                                   self._group_2,
                                                   self._align_dir,
                                                   self._event_dir,
                                                   self._asm_counts,
                                                   self._out_dir,
                                                   self._log_dir,
                                                   gtf=self._gtf)
        self._check_output()

    def _create_test_inputs(self):
        self._create_gtf()
        self._create_alignments()
        self._group_1 = tests.base_test.create_group_1_file(
            ['sample_1_1', 'sample_1_2'], self._input_dir)
        self._group_2 = tests.base_test.create_group_2_file(
            ['sample_2_1', 'sample_2_2'], self._input_dir)

    def _create_gtf(self):
        self.fail('implement in subclass')

    def _create_alignments(self):
        self.fail('implement in subclass')

    def _check_output(self):
        self.fail('implement in subclass')


class SkippedExonTest(BasicEventsBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'se'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1 includes all exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        gene.isoforms.append(isoform)
        # transcript 2 skips the 3rd exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        gene.isoforms.append(isoform)
        # transcript 3 skips the 2nd and 4th exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        gene.isoforms.append(isoform)
        # transcript 4 skips the 2nd, 3rd, and 4th exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        gene.isoforms.append(isoform)
        # gene_2, minus_strand
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1 includes all exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 2500))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3000, 3100))
        isoform.exons.append(tests.base_test.Region(3300, 3400))
        gene.isoforms.append(isoform)
        # transcript 2 skips the 3rd exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 2500))
        isoform.exons.append(tests.base_test.Region(3000, 3100))
        isoform.exons.append(tests.base_test.Region(3300, 3400))
        gene.isoforms.append(isoform)
        # transcript 3 skips the 2nd and 4th exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3300, 3400))
        gene.isoforms.append(isoform)
        # transcript 4 skips the 2nd, 3rd, and 4th exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(3300, 3400))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        self._chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, self._chrs)

    def _create_alignments(self):
        sample_1_1_sam = os.path.join(self._input_dir, 'g_1_1.sam')
        sample_1_2_sam = os.path.join(self._input_dir, 'g_1_2.sam')
        sample_2_1_sam = os.path.join(self._input_dir, 'g_2_1.sam')
        sample_2_2_sam = os.path.join(self._input_dir, 'g_2_2.sam')
        sams_by_sample = dict()
        sams_by_sample['sample_1_1'] = [sample_1_1_sam]
        sams_by_sample['sample_1_2'] = [sample_1_2_sam]
        sams_by_sample['sample_2_1'] = [sample_2_1_sam]
        sams_by_sample['sample_2_2'] = [sample_2_2_sam]
        self._sams_by_sample = sams_by_sample

        chr_1 = self._chrs[0]
        # 1,2,3,4,5
        all_isoform = chr_1.genes[0].isoforms[0]
        all_align = tests.base_test.Alignment()
        all_align.chr_name = chr_1.name
        all_align.match_to_isoform(all_isoform)

        # 1,2,4,5
        no_3_isoform = chr_1.genes[0].isoforms[1]
        no_3_align = tests.base_test.Alignment()
        no_3_align.chr_name = chr_1.name
        no_3_align.match_to_isoform(no_3_isoform)

        # 1,3,5
        no_2_4_isoform = chr_1.genes[0].isoforms[2]
        no_2_4_align = tests.base_test.Alignment()
        no_2_4_align.chr_name = chr_1.name
        no_2_4_align.match_to_isoform(no_2_4_isoform)

        # 1,5
        no_2_3_4_isoform = chr_1.genes[0].isoforms[3]
        no_2_3_4_align = tests.base_test.Alignment()
        no_2_3_4_align.chr_name = chr_1.name
        no_2_3_4_align.match_to_isoform(no_2_3_4_isoform)

        # gene_2
        # 1,2,3,4,5
        g2_all_isoform = chr_1.genes[1].isoforms[0]
        g2_all_align = tests.base_test.Alignment()
        g2_all_align.chr_name = chr_1.name
        g2_all_align.match_to_isoform(g2_all_isoform)

        # 1,2,4,5
        g2_no_3_isoform = chr_1.genes[1].isoforms[1]
        g2_no_3_align = tests.base_test.Alignment()
        g2_no_3_align.chr_name = chr_1.name
        g2_no_3_align.match_to_isoform(g2_no_3_isoform)

        # 1,3,5
        g2_no_2_4_isoform = chr_1.genes[1].isoforms[2]
        g2_no_2_4_align = tests.base_test.Alignment()
        g2_no_2_4_align.chr_name = chr_1.name
        g2_no_2_4_align.match_to_isoform(g2_no_2_4_isoform)

        # 1,5
        g2_no_2_3_4_isoform = chr_1.genes[1].isoforms[3]
        g2_no_2_3_4_align = tests.base_test.Alignment()
        g2_no_2_3_4_align.chr_name = chr_1.name
        g2_no_2_3_4_align.match_to_isoform(g2_no_2_3_4_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(all_align, 10, aligns)
        tests.base_test.append_copies(no_3_align, 20, aligns)
        tests.base_test.append_copies(no_2_4_align, 30, aligns)
        tests.base_test.append_copies(no_2_3_4_align, 40, aligns)
        tests.base_test.append_copies(g2_all_align, 10, aligns)
        tests.base_test.append_copies(g2_no_3_align, 20, aligns)
        tests.base_test.append_copies(g2_no_2_4_align, 30, aligns)
        tests.base_test.append_copies(g2_no_2_3_4_align, 40, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 10, aligns)
        tests.base_test.append_copies(no_3_align, 20, aligns)
        tests.base_test.append_copies(no_2_4_align, 30, aligns)
        tests.base_test.append_copies(no_2_3_4_align, 40, aligns)
        tests.base_test.append_copies(g2_all_align, 10, aligns)
        tests.base_test.append_copies(g2_no_3_align, 20, aligns)
        tests.base_test.append_copies(g2_no_2_4_align, 30, aligns)
        tests.base_test.append_copies(g2_no_2_3_4_align, 40, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 40, aligns)
        tests.base_test.append_copies(no_3_align, 10, aligns)
        tests.base_test.append_copies(no_2_4_align, 20, aligns)
        tests.base_test.append_copies(no_2_3_4_align, 30, aligns)
        tests.base_test.append_copies(g2_all_align, 40, aligns)
        tests.base_test.append_copies(g2_no_3_align, 10, aligns)
        tests.base_test.append_copies(g2_no_2_4_align, 20, aligns)
        tests.base_test.append_copies(g2_no_2_3_4_align, 30, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 40, aligns)
        tests.base_test.append_copies(no_3_align, 10, aligns)
        tests.base_test.append_copies(no_2_4_align, 20, aligns)
        tests.base_test.append_copies(no_2_3_4_align, 30, aligns)
        tests.base_test.append_copies(g2_all_align, 40, aligns)
        tests.base_test.append_copies(g2_no_3_align, 10, aligns)
        tests.base_test.append_copies(g2_no_2_4_align, 20, aligns)
        tests.base_test.append_copies(g2_no_2_3_4_align, 30, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    def _check_output(self):
        unfilt_path = os.path.join(self._rmats_long_out,
                                   'differential_isoforms.tsv')
        chr_id = '0'
        event_i = 0
        asm_123_id = '{}_{}'.format(chr_id, event_i)
        inc_123_isoform_id = '{}_0'.format(asm_123_id)
        skip_123_isoform_id = '{}_1'.format(asm_123_id)
        event_i += 1
        asm_135_id = '{}_{}'.format(chr_id, event_i)
        inc_135_isoform_id = '{}_0'.format(asm_135_id)
        skip_135_isoform_id = '{}_1'.format(asm_135_id)
        event_i += 1
        asm_234_id = '{}_{}'.format(chr_id, event_i)
        inc_234_isoform_id = '{}_0'.format(asm_234_id)
        skip_234_isoform_id = '{}_1'.format(asm_234_id)
        event_i += 1
        asm_345_id = '{}_{}'.format(chr_id, event_i)
        inc_345_isoform_id = '{}_0'.format(asm_345_id)
        skip_345_isoform_id = '{}_1'.format(asm_345_id)
        event_i += 1
        # There are 2 MXE events which are detected, but
        # only have reads matching 1 isoform. (1,2,4;1,3,4)(2,3,5;2,4,5)
        event_i += 2
        # gene_2
        g2_asm_123_id = '{}_{}'.format(chr_id, event_i)
        g2_inc_123_isoform_id = '{}_0'.format(g2_asm_123_id)
        g2_skip_123_isoform_id = '{}_1'.format(g2_asm_123_id)
        event_i += 1
        g2_asm_135_id = '{}_{}'.format(chr_id, event_i)
        g2_inc_135_isoform_id = '{}_0'.format(g2_asm_135_id)
        g2_skip_135_isoform_id = '{}_1'.format(g2_asm_135_id)
        event_i += 1
        g2_asm_234_id = '{}_{}'.format(chr_id, event_i)
        g2_inc_234_isoform_id = '{}_0'.format(g2_asm_234_id)
        g2_skip_234_isoform_id = '{}_1'.format(g2_asm_234_id)
        event_i += 1
        g2_asm_345_id = '{}_{}'.format(chr_id, event_i)
        g2_inc_345_isoform_id = '{}_0'.format(g2_asm_345_id)
        g2_skip_345_isoform_id = '{}_1'.format(g2_asm_345_id)
        event_i += 1

        unfilt_transcripts = self.parse_differential_transcripts(unfilt_path,
                                                                 is_asm=True)
        self.assertIn(inc_123_isoform_id, unfilt_transcripts)
        self.assertIn(skip_123_isoform_id, unfilt_transcripts)
        self.assertIn(inc_135_isoform_id, unfilt_transcripts)
        self.assertIn(skip_135_isoform_id, unfilt_transcripts)
        self.assertIn(inc_234_isoform_id, unfilt_transcripts)
        self.assertIn(skip_234_isoform_id, unfilt_transcripts)
        self.assertIn(inc_345_isoform_id, unfilt_transcripts)
        self.assertIn(skip_345_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_123_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_123_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_135_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_135_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_234_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_234_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_345_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_345_isoform_id, unfilt_transcripts)

        row = unfilt_transcripts.pop(inc_123_isoform_id)
        self.assertEqual(row['asm_id'], asm_123_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(skip_123_isoform_id)
        self.assertEqual(row['asm_id'], asm_123_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.416,
                               delta=0.001)

        row = unfilt_transcripts.pop(inc_135_isoform_id)
        self.assertEqual(row['asm_id'], asm_135_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.028,
                               delta=0.001)
        row = unfilt_transcripts.pop(skip_135_isoform_id)
        self.assertEqual(row['asm_id'], asm_135_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.028,
                               delta=0.001)

        row = unfilt_transcripts.pop(inc_234_isoform_id)
        self.assertEqual(row['asm_id'], asm_234_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.466,
                               delta=0.001)
        row = unfilt_transcripts.pop(skip_234_isoform_id)
        self.assertEqual(row['asm_id'], asm_234_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.466,
                               delta=0.001)

        row = unfilt_transcripts.pop(inc_345_isoform_id)
        self.assertEqual(row['asm_id'], asm_345_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(skip_345_isoform_id)
        self.assertEqual(row['asm_id'], asm_345_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.416,
                               delta=0.001)

        # gene_2
        row = unfilt_transcripts.pop(g2_inc_123_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_123_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_skip_123_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_123_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.416,
                               delta=0.001)

        row = unfilt_transcripts.pop(g2_inc_135_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_135_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.028,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_skip_135_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_135_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.028,
                               delta=0.001)

        row = unfilt_transcripts.pop(g2_inc_234_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_234_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.466,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_skip_234_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_234_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.466,
                               delta=0.001)

        row = unfilt_transcripts.pop(g2_inc_345_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_345_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_skip_345_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_345_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.416,
                               delta=0.001)

        self.assertEqual(unfilt_transcripts, dict())


class Alt5SpliceSiteTest(BasicEventsBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'a5ss'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # short
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(600, 700))
        gene.isoforms.append(isoform)
        # long
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 300))
        isoform.exons.append(tests.base_test.Region(600, 700))
        gene.isoforms.append(isoform)
        # longer
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 400))
        isoform.exons.append(tests.base_test.Region(600, 700))
        gene.isoforms.append(isoform)
        # gene_2, minus_strand
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # short
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        isoform.exons.append(tests.base_test.Region(1600, 1700))
        gene.isoforms.append(isoform)
        # long
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        isoform.exons.append(tests.base_test.Region(1500, 1700))
        gene.isoforms.append(isoform)
        # longer
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        isoform.exons.append(tests.base_test.Region(1400, 1700))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        self._chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, self._chrs)

    def _create_alignments(self):
        sample_1_1_sam = os.path.join(self._input_dir, 'g_1_1.sam')
        sample_1_2_sam = os.path.join(self._input_dir, 'g_1_2.sam')
        sample_2_1_sam = os.path.join(self._input_dir, 'g_2_1.sam')
        sample_2_2_sam = os.path.join(self._input_dir, 'g_2_2.sam')
        sams_by_sample = dict()
        sams_by_sample['sample_1_1'] = [sample_1_1_sam]
        sams_by_sample['sample_1_2'] = [sample_1_2_sam]
        sams_by_sample['sample_2_1'] = [sample_2_1_sam]
        sams_by_sample['sample_2_2'] = [sample_2_2_sam]
        self._sams_by_sample = sams_by_sample

        chr_1 = self._chrs[0]
        # 1s,2
        short_isoform = chr_1.genes[0].isoforms[0]
        short_align = tests.base_test.Alignment()
        short_align.chr_name = chr_1.name
        short_align.match_to_isoform(short_isoform)
        # 1l,2
        long_isoform = chr_1.genes[0].isoforms[1]
        long_align = tests.base_test.Alignment()
        long_align.chr_name = chr_1.name
        long_align.match_to_isoform(long_isoform)
        # 1ll,2
        longer_isoform = chr_1.genes[0].isoforms[2]
        longer_align = tests.base_test.Alignment()
        longer_align.chr_name = chr_1.name
        longer_align.match_to_isoform(longer_isoform)

        # gene_2
        # 1,s2
        g2_short_isoform = chr_1.genes[1].isoforms[0]
        g2_short_align = tests.base_test.Alignment()
        g2_short_align.chr_name = chr_1.name
        g2_short_align.match_to_isoform(g2_short_isoform)
        # 1,l2
        g2_long_isoform = chr_1.genes[1].isoforms[1]
        g2_long_align = tests.base_test.Alignment()
        g2_long_align.chr_name = chr_1.name
        g2_long_align.match_to_isoform(g2_long_isoform)
        # 1,ll2
        g2_longer_isoform = chr_1.genes[1].isoforms[2]
        g2_longer_align = tests.base_test.Alignment()
        g2_longer_align.chr_name = chr_1.name
        g2_longer_align.match_to_isoform(g2_longer_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(short_align, 10, aligns)
        tests.base_test.append_copies(long_align, 20, aligns)
        tests.base_test.append_copies(longer_align, 30, aligns)
        tests.base_test.append_copies(g2_short_align, 10, aligns)
        tests.base_test.append_copies(g2_long_align, 20, aligns)
        tests.base_test.append_copies(g2_longer_align, 30, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(short_align, 10, aligns)
        tests.base_test.append_copies(long_align, 20, aligns)
        tests.base_test.append_copies(longer_align, 30, aligns)
        tests.base_test.append_copies(g2_short_align, 10, aligns)
        tests.base_test.append_copies(g2_long_align, 20, aligns)
        tests.base_test.append_copies(g2_longer_align, 30, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(short_align, 20, aligns)
        tests.base_test.append_copies(long_align, 30, aligns)
        tests.base_test.append_copies(longer_align, 10, aligns)
        tests.base_test.append_copies(g2_short_align, 20, aligns)
        tests.base_test.append_copies(g2_long_align, 30, aligns)
        tests.base_test.append_copies(g2_longer_align, 10, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(short_align, 20, aligns)
        tests.base_test.append_copies(long_align, 30, aligns)
        tests.base_test.append_copies(longer_align, 10, aligns)
        tests.base_test.append_copies(g2_short_align, 20, aligns)
        tests.base_test.append_copies(g2_long_align, 30, aligns)
        tests.base_test.append_copies(g2_longer_align, 10, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += 1

    def _check_output(self):
        unfilt_path = os.path.join(self._rmats_long_out,
                                   'differential_isoforms.tsv')
        chr_id = '0'
        event_i = 0
        asm_short_long_id = '{}_{}'.format(chr_id, event_i)
        skip_short_long_isoform_id = '{}_0'.format(asm_short_long_id)
        inc_short_long_isoform_id = '{}_1'.format(asm_short_long_id)
        event_i += 1
        asm_short_longer_id = '{}_{}'.format(chr_id, event_i)
        skip_short_longer_isoform_id = '{}_0'.format(asm_short_longer_id)
        inc_short_longer_isoform_id = '{}_1'.format(asm_short_longer_id)
        event_i += 1
        asm_long_longer_id = '{}_{}'.format(chr_id, event_i)
        skip_long_longer_isoform_id = '{}_0'.format(asm_long_longer_id)
        inc_long_longer_isoform_id = '{}_1'.format(asm_long_longer_id)
        event_i += 1
        # gene_2
        g2_asm_short_long_id = '{}_{}'.format(chr_id, event_i)
        g2_inc_short_long_isoform_id = '{}_0'.format(g2_asm_short_long_id)
        g2_skip_short_long_isoform_id = '{}_1'.format(g2_asm_short_long_id)
        event_i += 1
        g2_asm_short_longer_id = '{}_{}'.format(chr_id, event_i)
        g2_inc_short_longer_isoform_id = '{}_0'.format(g2_asm_short_longer_id)
        g2_skip_short_longer_isoform_id = '{}_1'.format(g2_asm_short_longer_id)
        event_i += 1
        g2_asm_long_longer_id = '{}_{}'.format(chr_id, event_i)
        g2_inc_long_longer_isoform_id = '{}_0'.format(g2_asm_long_longer_id)
        g2_skip_long_longer_isoform_id = '{}_1'.format(g2_asm_long_longer_id)

        unfilt_transcripts = self.parse_differential_transcripts(unfilt_path,
                                                                 is_asm=True)
        self.assertIn(skip_short_long_isoform_id, unfilt_transcripts)
        self.assertIn(inc_short_long_isoform_id, unfilt_transcripts)
        self.assertIn(skip_short_longer_isoform_id, unfilt_transcripts)
        self.assertIn(inc_short_longer_isoform_id, unfilt_transcripts)
        self.assertIn(skip_long_longer_isoform_id, unfilt_transcripts)
        self.assertIn(inc_long_longer_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_short_long_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_short_long_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_short_longer_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_short_longer_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_long_longer_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_long_longer_isoform_id, unfilt_transcripts)

        row = unfilt_transcripts.pop(skip_short_long_isoform_id)
        self.assertEqual(row['asm_id'], asm_short_long_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.066,
                               delta=0.001)
        row = unfilt_transcripts.pop(inc_short_long_isoform_id)
        self.assertEqual(row['asm_id'], asm_short_long_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.066,
                               delta=0.001)
        row = unfilt_transcripts.pop(skip_short_longer_isoform_id)
        self.assertEqual(row['asm_id'], asm_short_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(inc_short_longer_isoform_id)
        self.assertEqual(row['asm_id'], asm_short_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(skip_long_longer_isoform_id)
        self.assertEqual(row['asm_id'], asm_long_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.35,
                               delta=0.001)
        row = unfilt_transcripts.pop(inc_long_longer_isoform_id)
        self.assertEqual(row['asm_id'], asm_long_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.35,
                               delta=0.001)

        # gene_2
        row = unfilt_transcripts.pop(g2_skip_short_long_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_short_long_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.066,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_inc_short_long_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_short_long_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.066,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_skip_short_longer_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_short_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_inc_short_longer_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_short_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_skip_long_longer_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_long_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.35,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_inc_long_longer_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_long_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.35,
                               delta=0.001)

        self.assertEqual(unfilt_transcripts, dict())


class Alt3SpliceSiteTest(BasicEventsBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'a3ss'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # short
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(600, 700))
        gene.isoforms.append(isoform)
        # long
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(500, 700))
        gene.isoforms.append(isoform)
        # longer
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 700))
        gene.isoforms.append(isoform)
        # gene_2, minus_strand
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # short
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        isoform.exons.append(tests.base_test.Region(1600, 1700))
        gene.isoforms.append(isoform)
        # long
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(1100, 1300))
        isoform.exons.append(tests.base_test.Region(1600, 1700))
        gene.isoforms.append(isoform)
        # longer
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(1100, 1400))
        isoform.exons.append(tests.base_test.Region(1600, 1700))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        self._chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, self._chrs)

    def _create_alignments(self):
        sample_1_1_sam = os.path.join(self._input_dir, 'g_1_1.sam')
        sample_1_2_sam = os.path.join(self._input_dir, 'g_1_2.sam')
        sample_2_1_sam = os.path.join(self._input_dir, 'g_2_1.sam')
        sample_2_2_sam = os.path.join(self._input_dir, 'g_2_2.sam')
        sams_by_sample = dict()
        sams_by_sample['sample_1_1'] = [sample_1_1_sam]
        sams_by_sample['sample_1_2'] = [sample_1_2_sam]
        sams_by_sample['sample_2_1'] = [sample_2_1_sam]
        sams_by_sample['sample_2_2'] = [sample_2_2_sam]
        self._sams_by_sample = sams_by_sample

        chr_1 = self._chrs[0]
        # 1,2s
        short_isoform = chr_1.genes[0].isoforms[0]
        short_align = tests.base_test.Alignment()
        short_align.chr_name = chr_1.name
        short_align.match_to_isoform(short_isoform)
        # 1,2l
        long_isoform = chr_1.genes[0].isoforms[1]
        long_align = tests.base_test.Alignment()
        long_align.chr_name = chr_1.name
        long_align.match_to_isoform(long_isoform)
        # 1,2ll
        longer_isoform = chr_1.genes[0].isoforms[2]
        longer_align = tests.base_test.Alignment()
        longer_align.chr_name = chr_1.name
        longer_align.match_to_isoform(longer_isoform)

        # gene_2
        # 1s,2
        g2_short_isoform = chr_1.genes[1].isoforms[0]
        g2_short_align = tests.base_test.Alignment()
        g2_short_align.chr_name = chr_1.name
        g2_short_align.match_to_isoform(g2_short_isoform)
        # 1l,2
        g2_long_isoform = chr_1.genes[1].isoforms[1]
        g2_long_align = tests.base_test.Alignment()
        g2_long_align.chr_name = chr_1.name
        g2_long_align.match_to_isoform(g2_long_isoform)
        # 1ll,2
        g2_longer_isoform = chr_1.genes[1].isoforms[2]
        g2_longer_align = tests.base_test.Alignment()
        g2_longer_align.chr_name = chr_1.name
        g2_longer_align.match_to_isoform(g2_longer_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(short_align, 10, aligns)
        tests.base_test.append_copies(long_align, 20, aligns)
        tests.base_test.append_copies(longer_align, 30, aligns)
        tests.base_test.append_copies(g2_short_align, 10, aligns)
        tests.base_test.append_copies(g2_long_align, 20, aligns)
        tests.base_test.append_copies(g2_longer_align, 30, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(short_align, 10, aligns)
        tests.base_test.append_copies(long_align, 20, aligns)
        tests.base_test.append_copies(longer_align, 30, aligns)
        tests.base_test.append_copies(g2_short_align, 10, aligns)
        tests.base_test.append_copies(g2_long_align, 20, aligns)
        tests.base_test.append_copies(g2_longer_align, 30, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(short_align, 20, aligns)
        tests.base_test.append_copies(long_align, 30, aligns)
        tests.base_test.append_copies(longer_align, 10, aligns)
        tests.base_test.append_copies(g2_short_align, 20, aligns)
        tests.base_test.append_copies(g2_long_align, 30, aligns)
        tests.base_test.append_copies(g2_longer_align, 10, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(short_align, 20, aligns)
        tests.base_test.append_copies(long_align, 30, aligns)
        tests.base_test.append_copies(longer_align, 10, aligns)
        tests.base_test.append_copies(g2_short_align, 20, aligns)
        tests.base_test.append_copies(g2_long_align, 30, aligns)
        tests.base_test.append_copies(g2_longer_align, 10, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += 1

    def _check_output(self):
        unfilt_path = os.path.join(self._rmats_long_out,
                                   'differential_isoforms.tsv')
        chr_id = '0'
        event_i = 0
        asm_long_longer_id = '{}_{}'.format(chr_id, event_i)
        inc_long_longer_isoform_id = '{}_0'.format(asm_long_longer_id)
        skip_long_longer_isoform_id = '{}_1'.format(asm_long_longer_id)
        event_i += 1
        asm_short_longer_id = '{}_{}'.format(chr_id, event_i)
        inc_short_longer_isoform_id = '{}_0'.format(asm_short_longer_id)
        skip_short_longer_isoform_id = '{}_1'.format(asm_short_longer_id)
        event_i += 1
        asm_short_long_id = '{}_{}'.format(chr_id, event_i)
        inc_short_long_isoform_id = '{}_0'.format(asm_short_long_id)
        skip_short_long_isoform_id = '{}_1'.format(asm_short_long_id)
        event_i += 1
        # gene_2
        g2_asm_long_longer_id = '{}_{}'.format(chr_id, event_i)
        g2_skip_long_longer_isoform_id = '{}_0'.format(g2_asm_long_longer_id)
        g2_inc_long_longer_isoform_id = '{}_1'.format(g2_asm_long_longer_id)
        event_i += 1
        g2_asm_short_longer_id = '{}_{}'.format(chr_id, event_i)
        g2_skip_short_longer_isoform_id = '{}_0'.format(g2_asm_short_longer_id)
        g2_inc_short_longer_isoform_id = '{}_1'.format(g2_asm_short_longer_id)
        event_i += 1
        g2_asm_short_long_id = '{}_{}'.format(chr_id, event_i)
        g2_skip_short_long_isoform_id = '{}_0'.format(g2_asm_short_long_id)
        g2_inc_short_long_isoform_id = '{}_1'.format(g2_asm_short_long_id)

        unfilt_transcripts = self.parse_differential_transcripts(unfilt_path,
                                                                 is_asm=True)
        self.assertIn(inc_long_longer_isoform_id, unfilt_transcripts)
        self.assertIn(skip_long_longer_isoform_id, unfilt_transcripts)
        self.assertIn(inc_short_longer_isoform_id, unfilt_transcripts)
        self.assertIn(skip_short_longer_isoform_id, unfilt_transcripts)
        self.assertIn(inc_short_long_isoform_id, unfilt_transcripts)
        self.assertIn(skip_short_long_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_long_longer_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_long_longer_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_short_longer_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_short_longer_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_short_long_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_short_long_isoform_id, unfilt_transcripts)

        row = unfilt_transcripts.pop(inc_long_longer_isoform_id)
        self.assertEqual(row['asm_id'], asm_long_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.35,
                               delta=0.001)
        row = unfilt_transcripts.pop(skip_long_longer_isoform_id)
        self.assertEqual(row['asm_id'], asm_long_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.35,
                               delta=0.001)
        row = unfilt_transcripts.pop(inc_short_longer_isoform_id)
        self.assertEqual(row['asm_id'], asm_short_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(skip_short_longer_isoform_id)
        self.assertEqual(row['asm_id'], asm_short_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(inc_short_long_isoform_id)
        self.assertEqual(row['asm_id'], asm_short_long_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.066,
                               delta=0.001)
        row = unfilt_transcripts.pop(skip_short_long_isoform_id)
        self.assertEqual(row['asm_id'], asm_short_long_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.066,
                               delta=0.001)

        # gene_2
        row = unfilt_transcripts.pop(g2_inc_long_longer_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_long_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.35,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_skip_long_longer_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_long_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.35,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_inc_short_longer_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_short_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_skip_short_longer_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_short_longer_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_inc_short_long_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_short_long_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.066,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_skip_short_long_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_short_long_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.066,
                               delta=0.001)

        self.assertEqual(unfilt_transcripts, dict())


class MutuallyExclusiveExonsTest(BasicEventsBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'mxe'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # 1,2,5
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        gene.isoforms.append(isoform)
        # 1,3,5
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        gene.isoforms.append(isoform)
        # 1,4,5
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        gene.isoforms.append(isoform)
        # gene_2, minus_strand
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # 1,2,5
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(3000, 3100))
        isoform.exons.append(tests.base_test.Region(3300, 3400))
        gene.isoforms.append(isoform)
        # 1,3,5
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3300, 3400))
        gene.isoforms.append(isoform)
        # 1,4,5
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 2500))
        isoform.exons.append(tests.base_test.Region(3300, 3400))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        self._chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, self._chrs)

    def _create_alignments(self):
        sample_1_1_sam = os.path.join(self._input_dir, 'g_1_1.sam')
        sample_1_2_sam = os.path.join(self._input_dir, 'g_1_2.sam')
        sample_2_1_sam = os.path.join(self._input_dir, 'g_2_1.sam')
        sample_2_2_sam = os.path.join(self._input_dir, 'g_2_2.sam')
        sams_by_sample = dict()
        sams_by_sample['sample_1_1'] = [sample_1_1_sam]
        sams_by_sample['sample_1_2'] = [sample_1_2_sam]
        sams_by_sample['sample_2_1'] = [sample_2_1_sam]
        sams_by_sample['sample_2_2'] = [sample_2_2_sam]
        self._sams_by_sample = sams_by_sample

        chr_1 = self._chrs[0]
        # 1,2,5
        use_second_isoform = chr_1.genes[0].isoforms[0]
        use_second_align = tests.base_test.Alignment()
        use_second_align.chr_name = chr_1.name
        use_second_align.match_to_isoform(use_second_isoform)
        # 1,3,5
        use_third_isoform = chr_1.genes[0].isoforms[1]
        use_third_align = tests.base_test.Alignment()
        use_third_align.chr_name = chr_1.name
        use_third_align.match_to_isoform(use_third_isoform)
        # 1,4,5
        use_fourth_isoform = chr_1.genes[0].isoforms[2]
        use_fourth_align = tests.base_test.Alignment()
        use_fourth_align.chr_name = chr_1.name
        use_fourth_align.match_to_isoform(use_fourth_isoform)

        # gene_2
        # 1,2,5
        g2_use_second_isoform = chr_1.genes[1].isoforms[0]
        g2_use_second_align = tests.base_test.Alignment()
        g2_use_second_align.chr_name = chr_1.name
        g2_use_second_align.match_to_isoform(g2_use_second_isoform)
        # 1,3,5
        g2_use_third_isoform = chr_1.genes[1].isoforms[1]
        g2_use_third_align = tests.base_test.Alignment()
        g2_use_third_align.chr_name = chr_1.name
        g2_use_third_align.match_to_isoform(g2_use_third_isoform)
        # 1,4,5
        g2_use_fourth_isoform = chr_1.genes[1].isoforms[2]
        g2_use_fourth_align = tests.base_test.Alignment()
        g2_use_fourth_align.chr_name = chr_1.name
        g2_use_fourth_align.match_to_isoform(g2_use_fourth_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(use_second_align, 10, aligns)
        tests.base_test.append_copies(use_third_align, 20, aligns)
        tests.base_test.append_copies(use_fourth_align, 30, aligns)
        tests.base_test.append_copies(g2_use_second_align, 10, aligns)
        tests.base_test.append_copies(g2_use_third_align, 20, aligns)
        tests.base_test.append_copies(g2_use_fourth_align, 30, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(use_second_align, 10, aligns)
        tests.base_test.append_copies(use_third_align, 20, aligns)
        tests.base_test.append_copies(use_fourth_align, 30, aligns)
        tests.base_test.append_copies(g2_use_second_align, 10, aligns)
        tests.base_test.append_copies(g2_use_third_align, 20, aligns)
        tests.base_test.append_copies(g2_use_fourth_align, 30, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(use_second_align, 20, aligns)
        tests.base_test.append_copies(use_third_align, 30, aligns)
        tests.base_test.append_copies(use_fourth_align, 10, aligns)
        tests.base_test.append_copies(g2_use_second_align, 20, aligns)
        tests.base_test.append_copies(g2_use_third_align, 30, aligns)
        tests.base_test.append_copies(g2_use_fourth_align, 10, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(use_second_align, 20, aligns)
        tests.base_test.append_copies(use_third_align, 30, aligns)
        tests.base_test.append_copies(use_fourth_align, 10, aligns)
        tests.base_test.append_copies(g2_use_second_align, 20, aligns)
        tests.base_test.append_copies(g2_use_third_align, 30, aligns)
        tests.base_test.append_copies(g2_use_fourth_align, 10, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += 1

    def _check_output(self):
        unfilt_path = os.path.join(self._rmats_long_out,
                                   'differential_isoforms.tsv')
        chr_id = '0'
        event_i = 0
        asm_2_3_id = '{}_{}'.format(chr_id, event_i)
        inc_2_3_isoform_id = '{}_0'.format(asm_2_3_id)
        skip_2_3_isoform_id = '{}_1'.format(asm_2_3_id)
        event_i += 1
        asm_2_4_id = '{}_{}'.format(chr_id, event_i)
        inc_2_4_isoform_id = '{}_0'.format(asm_2_4_id)
        skip_2_4_isoform_id = '{}_1'.format(asm_2_4_id)
        event_i += 1
        asm_3_4_id = '{}_{}'.format(chr_id, event_i)
        inc_3_4_isoform_id = '{}_0'.format(asm_3_4_id)
        skip_3_4_isoform_id = '{}_1'.format(asm_3_4_id)
        event_i += 1
        # gene_2
        g2_asm_2_3_id = '{}_{}'.format(chr_id, event_i)
        g2_skip_2_3_isoform_id = '{}_0'.format(g2_asm_2_3_id)
        g2_inc_2_3_isoform_id = '{}_1'.format(g2_asm_2_3_id)
        event_i += 1
        g2_asm_2_4_id = '{}_{}'.format(chr_id, event_i)
        g2_skip_2_4_isoform_id = '{}_0'.format(g2_asm_2_4_id)
        g2_inc_2_4_isoform_id = '{}_1'.format(g2_asm_2_4_id)
        event_i += 1
        g2_asm_3_4_id = '{}_{}'.format(chr_id, event_i)
        g2_skip_3_4_isoform_id = '{}_0'.format(g2_asm_3_4_id)
        g2_inc_3_4_isoform_id = '{}_1'.format(g2_asm_3_4_id)

        unfilt_transcripts = self.parse_differential_transcripts(unfilt_path,
                                                                 is_asm=True)
        self.assertIn(inc_2_3_isoform_id, unfilt_transcripts)
        self.assertIn(skip_2_3_isoform_id, unfilt_transcripts)
        self.assertIn(inc_2_4_isoform_id, unfilt_transcripts)
        self.assertIn(skip_2_4_isoform_id, unfilt_transcripts)
        self.assertIn(inc_3_4_isoform_id, unfilt_transcripts)
        self.assertIn(skip_3_4_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_2_3_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_2_3_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_2_4_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_2_4_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_3_4_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_3_4_isoform_id, unfilt_transcripts)

        row = unfilt_transcripts.pop(inc_2_3_isoform_id)
        self.assertEqual(row['asm_id'], asm_2_3_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.066,
                               delta=0.001)
        row = unfilt_transcripts.pop(skip_2_3_isoform_id)
        self.assertEqual(row['asm_id'], asm_2_3_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.066,
                               delta=0.001)
        row = unfilt_transcripts.pop(inc_2_4_isoform_id)
        self.assertEqual(row['asm_id'], asm_2_4_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(skip_2_4_isoform_id)
        self.assertEqual(row['asm_id'], asm_2_4_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(inc_3_4_isoform_id)
        self.assertEqual(row['asm_id'], asm_3_4_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.35,
                               delta=0.001)
        row = unfilt_transcripts.pop(skip_3_4_isoform_id)
        self.assertEqual(row['asm_id'], asm_3_4_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.35,
                               delta=0.001)

        # gene_2
        row = unfilt_transcripts.pop(g2_inc_2_3_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_2_3_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.066,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_skip_2_3_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_2_3_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.066,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_inc_2_4_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_2_4_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_skip_2_4_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_2_4_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_inc_3_4_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_3_4_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.35,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_skip_3_4_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_3_4_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.35,
                               delta=0.001)

        self.assertEqual(unfilt_transcripts, dict())


class RetainedIntronTest(BasicEventsBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'ri'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # 1,2,3,4
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        gene.isoforms.append(isoform)
        # 1-3,4
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 800))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        gene.isoforms.append(isoform)
        # 1,2,3-4
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(700, 1100))
        gene.isoforms.append(isoform)
        # gene_2, minus_strand
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # 1,2,3,4
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 2500))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3000, 3100))
        gene.isoforms.append(isoform)
        # 1-3,4
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 3100))
        gene.isoforms.append(isoform)
        # 1,2,3-4
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2500))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3000, 3100))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        self._chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, self._chrs)

    def _create_alignments(self):
        sample_1_1_sam = os.path.join(self._input_dir, 'g_1_1.sam')
        sample_1_2_sam = os.path.join(self._input_dir, 'g_1_2.sam')
        sample_2_1_sam = os.path.join(self._input_dir, 'g_2_1.sam')
        sample_2_2_sam = os.path.join(self._input_dir, 'g_2_2.sam')
        sams_by_sample = dict()
        sams_by_sample['sample_1_1'] = [sample_1_1_sam]
        sams_by_sample['sample_1_2'] = [sample_1_2_sam]
        sams_by_sample['sample_2_1'] = [sample_2_1_sam]
        sams_by_sample['sample_2_2'] = [sample_2_2_sam]
        self._sams_by_sample = sams_by_sample

        chr_1 = self._chrs[0]
        # 1,2,3,4
        all_exons_isoform = chr_1.genes[0].isoforms[0]
        all_exons_align = tests.base_test.Alignment()
        all_exons_align.chr_name = chr_1.name
        all_exons_align.match_to_isoform(all_exons_isoform)
        # 1-3,4
        retain_1_3_isoform = chr_1.genes[0].isoforms[1]
        retain_1_3_align = tests.base_test.Alignment()
        retain_1_3_align.chr_name = chr_1.name
        retain_1_3_align.match_to_isoform(retain_1_3_isoform)
        # 1,2,3-4
        retain_3_4_isoform = chr_1.genes[0].isoforms[2]
        retain_3_4_align = tests.base_test.Alignment()
        retain_3_4_align.chr_name = chr_1.name
        retain_3_4_align.match_to_isoform(retain_3_4_isoform)

        # gene_2
        # 1,2,3,4
        g2_all_exons_isoform = chr_1.genes[1].isoforms[0]
        g2_all_exons_align = tests.base_test.Alignment()
        g2_all_exons_align.chr_name = chr_1.name
        g2_all_exons_align.match_to_isoform(g2_all_exons_isoform)
        # 1-3,4
        g2_retain_1_3_isoform = chr_1.genes[1].isoforms[1]
        g2_retain_1_3_align = tests.base_test.Alignment()
        g2_retain_1_3_align.chr_name = chr_1.name
        g2_retain_1_3_align.match_to_isoform(g2_retain_1_3_isoform)
        # 1,2,3-4
        g2_retain_3_4_isoform = chr_1.genes[1].isoforms[2]
        g2_retain_3_4_align = tests.base_test.Alignment()
        g2_retain_3_4_align.chr_name = chr_1.name
        g2_retain_3_4_align.match_to_isoform(g2_retain_3_4_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(all_exons_align, 10, aligns)
        tests.base_test.append_copies(retain_1_3_align, 20, aligns)
        tests.base_test.append_copies(retain_3_4_align, 30, aligns)
        tests.base_test.append_copies(g2_all_exons_align, 10, aligns)
        tests.base_test.append_copies(g2_retain_1_3_align, 20, aligns)
        tests.base_test.append_copies(g2_retain_3_4_align, 30, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(all_exons_align, 10, aligns)
        tests.base_test.append_copies(retain_1_3_align, 20, aligns)
        tests.base_test.append_copies(retain_3_4_align, 30, aligns)
        tests.base_test.append_copies(g2_all_exons_align, 10, aligns)
        tests.base_test.append_copies(g2_retain_1_3_align, 20, aligns)
        tests.base_test.append_copies(g2_retain_3_4_align, 30, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(all_exons_align, 20, aligns)
        tests.base_test.append_copies(retain_1_3_align, 30, aligns)
        tests.base_test.append_copies(retain_3_4_align, 10, aligns)
        tests.base_test.append_copies(g2_all_exons_align, 20, aligns)
        tests.base_test.append_copies(g2_retain_1_3_align, 30, aligns)
        tests.base_test.append_copies(g2_retain_3_4_align, 10, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(all_exons_align, 20, aligns)
        tests.base_test.append_copies(retain_1_3_align, 30, aligns)
        tests.base_test.append_copies(retain_3_4_align, 10, aligns)
        tests.base_test.append_copies(g2_all_exons_align, 20, aligns)
        tests.base_test.append_copies(g2_retain_1_3_align, 30, aligns)
        tests.base_test.append_copies(g2_retain_3_4_align, 10, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += 1

    def _check_output(self):
        unfilt_path = os.path.join(self._rmats_long_out,
                                   'differential_isoforms.tsv')
        chr_id = '0'
        event_i = 0
        asm_3_4_id = '{}_{}'.format(chr_id, event_i)
        skip_3_4_isoform_id = '{}_0'.format(asm_3_4_id)
        inc_3_4_isoform_id = '{}_1'.format(asm_3_4_id)
        event_i += 1
        # gene_2
        g2_asm_3_4_id = '{}_{}'.format(chr_id, event_i)
        g2_skip_3_4_isoform_id = '{}_0'.format(g2_asm_3_4_id)
        g2_inc_3_4_isoform_id = '{}_1'.format(g2_asm_3_4_id)

        unfilt_transcripts = self.parse_differential_transcripts(unfilt_path,
                                                                 is_asm=True)
        self.assertIn(inc_3_4_isoform_id, unfilt_transcripts)
        self.assertIn(skip_3_4_isoform_id, unfilt_transcripts)
        self.assertIn(g2_inc_3_4_isoform_id, unfilt_transcripts)
        self.assertIn(g2_skip_3_4_isoform_id, unfilt_transcripts)

        row = unfilt_transcripts.pop(skip_3_4_isoform_id)
        self.assertEqual(row['asm_id'], asm_3_4_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(inc_3_4_isoform_id)
        self.assertEqual(row['asm_id'], asm_3_4_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.416,
                               delta=0.001)

        # gene_2
        row = unfilt_transcripts.pop(g2_skip_3_4_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_3_4_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.416,
                               delta=0.001)
        row = unfilt_transcripts.pop(g2_inc_3_4_isoform_id)
        self.assertEqual(row['asm_id'], g2_asm_3_4_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.416,
                               delta=0.001)

        self.assertEqual(unfilt_transcripts, dict())


class Both5And3Test(BasicEventsBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'both_5_and_3'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # Flank exons are fixed.
        # The middle exon has 3 each of A5SS and A3SS.
        # 1,1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 700))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)
        # 1,2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 800))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)
        # 1,3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 900))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)
        # 2,1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(500, 700))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)
        # 2,2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_5'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(500, 800))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)
        # 2,3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_6'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(500, 900))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)
        # 3,1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_7'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(600, 700))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)
        # 3,2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_8'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(600, 800))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)
        # 3,3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_9'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(600, 900))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        self._chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, self._chrs)

    def _create_alignments(self):
        sample_1_1_sam = os.path.join(self._input_dir, 'g_1_1.sam')
        sample_1_2_sam = os.path.join(self._input_dir, 'g_1_2.sam')
        sample_2_1_sam = os.path.join(self._input_dir, 'g_2_1.sam')
        sample_2_2_sam = os.path.join(self._input_dir, 'g_2_2.sam')
        sams_by_sample = dict()
        sams_by_sample['sample_1_1'] = [sample_1_1_sam]
        sams_by_sample['sample_1_2'] = [sample_1_2_sam]
        sams_by_sample['sample_2_1'] = [sample_2_1_sam]
        sams_by_sample['sample_2_2'] = [sample_2_2_sam]
        self._sams_by_sample = sams_by_sample

        chr_1 = self._chrs[0]
        # 1,1
        one_1_isoform = chr_1.genes[0].isoforms[0]
        one_1_align = tests.base_test.Alignment()
        one_1_align.chr_name = chr_1.name
        one_1_align.match_to_isoform(one_1_isoform)
        # 1,2
        one_2_isoform = chr_1.genes[0].isoforms[1]
        one_2_align = tests.base_test.Alignment()
        one_2_align.chr_name = chr_1.name
        one_2_align.match_to_isoform(one_2_isoform)
        # 1,3
        one_3_isoform = chr_1.genes[0].isoforms[2]
        one_3_align = tests.base_test.Alignment()
        one_3_align.chr_name = chr_1.name
        one_3_align.match_to_isoform(one_3_isoform)
        # 2,1
        two_1_isoform = chr_1.genes[0].isoforms[3]
        two_1_align = tests.base_test.Alignment()
        two_1_align.chr_name = chr_1.name
        two_1_align.match_to_isoform(two_1_isoform)
        # 2,2
        two_2_isoform = chr_1.genes[0].isoforms[4]
        two_2_align = tests.base_test.Alignment()
        two_2_align.chr_name = chr_1.name
        two_2_align.match_to_isoform(two_2_isoform)
        # 2,3
        two_3_isoform = chr_1.genes[0].isoforms[5]
        two_3_align = tests.base_test.Alignment()
        two_3_align.chr_name = chr_1.name
        two_3_align.match_to_isoform(two_3_isoform)
        # 3,1
        three_1_isoform = chr_1.genes[0].isoforms[6]
        three_1_align = tests.base_test.Alignment()
        three_1_align.chr_name = chr_1.name
        three_1_align.match_to_isoform(three_1_isoform)
        # 3,2
        three_2_isoform = chr_1.genes[0].isoforms[7]
        three_2_align = tests.base_test.Alignment()
        three_2_align.chr_name = chr_1.name
        three_2_align.match_to_isoform(three_2_isoform)
        # 3,3
        three_3_isoform = chr_1.genes[0].isoforms[8]
        three_3_align = tests.base_test.Alignment()
        three_3_align.chr_name = chr_1.name
        three_3_align.match_to_isoform(three_3_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(one_1_align, 5, aligns)
        tests.base_test.append_copies(one_2_align, 10, aligns)
        tests.base_test.append_copies(one_3_align, 15, aligns)
        tests.base_test.append_copies(two_1_align, 20, aligns)
        tests.base_test.append_copies(two_2_align, 25, aligns)
        tests.base_test.append_copies(two_3_align, 30, aligns)
        tests.base_test.append_copies(three_1_align, 35, aligns)
        tests.base_test.append_copies(three_2_align, 40, aligns)
        tests.base_test.append_copies(three_3_align, 45, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(one_1_align, 5, aligns)
        tests.base_test.append_copies(one_2_align, 10, aligns)
        tests.base_test.append_copies(one_3_align, 15, aligns)
        tests.base_test.append_copies(two_1_align, 20, aligns)
        tests.base_test.append_copies(two_2_align, 25, aligns)
        tests.base_test.append_copies(two_3_align, 30, aligns)
        tests.base_test.append_copies(three_1_align, 35, aligns)
        tests.base_test.append_copies(three_2_align, 40, aligns)
        tests.base_test.append_copies(three_3_align, 45, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(one_1_align, 45, aligns)
        tests.base_test.append_copies(one_2_align, 40, aligns)
        tests.base_test.append_copies(one_3_align, 35, aligns)
        tests.base_test.append_copies(two_1_align, 30, aligns)
        tests.base_test.append_copies(two_2_align, 25, aligns)
        tests.base_test.append_copies(two_3_align, 20, aligns)
        tests.base_test.append_copies(three_1_align, 15, aligns)
        tests.base_test.append_copies(three_2_align, 10, aligns)
        tests.base_test.append_copies(three_3_align, 5, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += 1

        aligns = list()
        tests.base_test.append_copies(one_1_align, 45, aligns)
        tests.base_test.append_copies(one_2_align, 40, aligns)
        tests.base_test.append_copies(one_3_align, 35, aligns)
        tests.base_test.append_copies(two_1_align, 30, aligns)
        tests.base_test.append_copies(two_2_align, 25, aligns)
        tests.base_test.append_copies(two_3_align, 20, aligns)
        tests.base_test.append_copies(three_1_align, 15, aligns)
        tests.base_test.append_copies(three_2_align, 10, aligns)
        tests.base_test.append_copies(three_3_align, 5, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += 1

    def _check_output(self):
        unfilt_path = os.path.join(self._rmats_long_out,
                                   'differential_isoforms.tsv')
        unfilt_transcripts = self.parse_differential_transcripts(unfilt_path,
                                                                 is_asm=True)
        self.assertEqual(len(unfilt_transcripts), 36)

        asms = self.parse_asm_definitions(self._event_dir)
        defs_by_type = dict()
        for asm in asms.values():
            asm_type = asm['event_type']
            old_count = defs_by_type.get(asm_type, 0)
            defs_by_type[asm_type] = old_count + 1

        expected_defs_by_type = {
            "alternative 5'-splice site usage": 9,
            "alternative 3'-splice site usage": 9,
        }
        self.assertEqual(defs_by_type, expected_defs_by_type)

        by_gene_dir = os.path.join(self._rmats_long_out, 'results_by_gene')
        gene_id = 'gene_1'
        gene_dir = os.path.join(by_gene_dir, gene_id)
        file_names = os.listdir(gene_dir)
        events_by_type = dict()
        for file_name in file_names:
            if '_isoform_differences_' not in file_name:
                continue

            diff_path = os.path.join(gene_dir, file_name)
            isoform_diffs = self.parse_isoform_diff_tsv(diff_path)
            self.assertEqual(len(isoform_diffs), 1)
            diff_row = isoform_diffs[0]
            event_type = diff_row['event']
            old_count = events_by_type.get(event_type, 0)
            events_by_type[event_type] = old_count + 1

        self.assertEqual(set(events_by_type.keys()), {'A5SS', 'A3SS'})


if __name__ == '__main__':
    unittest.main(verbosity=2)
