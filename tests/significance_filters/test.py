import os
import os.path
import unittest

import tests.base_test


class SignificanceFiltersBaseTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_dir = os.path.dirname(__file__)
        self._input_dir = None
        self._out_dir = None
        self._log_dir = None
        self._gtf = None
        self._sams_by_sample = None
        self._avg_cpm_group = None
        self._min_cpm_group = None
        self._asm_prop_gene = None
        self._avg_reads_group = None

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
        self._event_dir = self.run_detect_events(self._gtf_dir, self._out_dir,
                                                 self._log_dir)
        self._asm_counts = self.run_count_reads(self._gtf_dir, self._align_dir,
                                                self._event_dir, self._out_dir,
                                                self._log_dir)
        self._rmats_long_out = self.run_rmats_long(
            self._group_1,
            self._group_2,
            self._align_dir,
            self._event_dir,
            self._asm_counts,
            self._out_dir,
            self._log_dir,
            gtf=self._gtf,
            average_cpm_per_group=self._avg_cpm_group,
            min_cpm_per_group=self._min_cpm_group,
            asm_proportion_gene=self._asm_prop_gene,
            average_reads_per_group=self._avg_reads_group)
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


class SignificanceFiltersTest(SignificanceFiltersBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'sig_filters'
        self._set_input_out_log_with_prefix(self._prefix)
        self._asm_prop_gene = 0.1

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
        isoform.exons.append(tests.base_test.Region(1600, 1700))
        isoform.exons.append(tests.base_test.Region(1900, 2000))
        gene.isoforms.append(isoform)
        # transcript 2 skips the 5th exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        isoform.exons.append(tests.base_test.Region(1600, 1700))
        isoform.exons.append(tests.base_test.Region(1900, 2000))
        gene.isoforms.append(isoform)
        # transcript 3 skips exons 4,5,6
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1900, 2000))
        gene.isoforms.append(isoform)
        # transcript 4 skips exon 2 and also 4,5,6
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1900, 2000))
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
        # 1,2,3,4,5,6,7
        all_isoform = chr_1.genes[0].isoforms[0]
        all_align = tests.base_test.Alignment()
        all_align.chr_name = chr_1.name
        all_align.match_to_isoform(all_isoform)

        # 1,2,3,4,6,7
        no_5_isoform = chr_1.genes[0].isoforms[1]
        no_5_align = tests.base_test.Alignment()
        no_5_align.chr_name = chr_1.name
        no_5_align.match_to_isoform(no_5_isoform)

        # 1,2,3,7
        no_456_isoform = chr_1.genes[0].isoforms[2]
        no_456_align = tests.base_test.Alignment()
        no_456_align.chr_name = chr_1.name
        no_456_align.match_to_isoform(no_456_isoform)

        # 1,3,7
        only_137_isoform = chr_1.genes[0].isoforms[3]
        only_137_align = tests.base_test.Alignment()
        only_137_align.chr_name = chr_1.name
        only_137_align.match_to_isoform(only_137_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(all_align, 25, aligns)
        tests.base_test.append_copies(no_5_align, 25, aligns)
        tests.base_test.append_copies(no_456_align, 700, aligns)
        tests.base_test.append_copies(only_137_align, 250, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 25, aligns)
        tests.base_test.append_copies(no_5_align, 25, aligns)
        tests.base_test.append_copies(no_456_align, 700, aligns)
        tests.base_test.append_copies(only_137_align, 250, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 10, aligns)
        tests.base_test.append_copies(no_5_align, 40, aligns)
        tests.base_test.append_copies(no_456_align, 600, aligns)
        tests.base_test.append_copies(only_137_align, 350, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 10, aligns)
        tests.base_test.append_copies(no_5_align, 40, aligns)
        tests.base_test.append_copies(no_456_align, 600, aligns)
        tests.base_test.append_copies(only_137_align, 350, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    def _check_output(self):
        unfilt_path = os.path.join(self._rmats_long_out,
                                   'differential_isoforms.tsv')
        filt_path = os.path.join(self._rmats_long_out,
                                 'differential_isoforms_filtered.tsv')
        chr_id = '0'
        event_i = 0
        four_5_6_asm_id = '{}_{}'.format(chr_id, event_i)
        inc_456_isoform_id = '{}_0'.format(four_5_6_asm_id)
        inc_46_isoform_id = '{}_1'.format(four_5_6_asm_id)
        skip_456_isoform_id = '{}_2'.format(four_5_6_asm_id)

        event_i += 1
        exon_2_asm_id = '{}_{}'.format(chr_id, event_i)
        inc_2_isoform_id = '{}_0'.format(exon_2_asm_id)
        skip_2_isoform_id = '{}_1'.format(exon_2_asm_id)

        event_i += 1
        exon_5_asm_id = '{}_{}'.format(chr_id, event_i)
        inc_5_isoform_id = '{}_0'.format(exon_5_asm_id)
        skip_5_isoform_id = '{}_1'.format(exon_5_asm_id)

        unfilt_transcripts = self.parse_differential_transcripts(unfilt_path,
                                                                 is_asm=True)
        filt_transcripts = self.parse_differential_transcripts(filt_path,
                                                               is_asm=True)
        self.assertIn(inc_2_isoform_id, unfilt_transcripts)
        self.assertIn(skip_2_isoform_id, unfilt_transcripts)
        self.assertIn(inc_456_isoform_id, unfilt_transcripts)
        self.assertIn(inc_46_isoform_id, unfilt_transcripts)
        self.assertIn(skip_456_isoform_id, unfilt_transcripts)
        self.assertIn(inc_5_isoform_id, unfilt_transcripts)
        self.assertIn(skip_5_isoform_id, unfilt_transcripts)

        self.assertIn(inc_2_isoform_id, filt_transcripts)
        row = unfilt_transcripts[inc_2_isoform_id]
        self.assertEqual(row['asm_id'], exon_2_asm_id)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               0.1,
                               delta=0.001)

        self.assertNotIn(inc_46_isoform_id, filt_transcripts)
        row = unfilt_transcripts[inc_46_isoform_id]
        self.assertEqual(row['asm_id'], four_5_6_asm_id)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.015,
                               delta=0.001)

        self.assertNotIn(skip_5_isoform_id, filt_transcripts)
        row = unfilt_transcripts[skip_5_isoform_id]
        self.assertEqual(row['asm_id'], exon_5_asm_id)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.3,
                               delta=0.001)


class AverageCpmTest(SignificanceFiltersBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'avg_cpm'
        self._set_input_out_log_with_prefix(self._prefix)
        self._avg_cpm_group = 100_000  # 10%

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1 has 4 exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        gene.isoforms.append(isoform)
        # transcript 2 skips exon 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        gene.isoforms.append(isoform)
        # transcript 3 has skips exon 3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        gene.isoforms.append(isoform)

        # gene 2
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1 has 4 exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 2500))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3000, 3100))
        gene.isoforms.append(isoform)
        # transcript 2 skips exon 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3000, 3100))
        gene.isoforms.append(isoform)
        # transcript 3 has skips exon 3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 2500))
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
        all_isoform = chr_1.genes[0].isoforms[0]
        all_align = tests.base_test.Alignment()
        all_align.chr_name = chr_1.name
        all_align.match_to_isoform(all_isoform)

        skip_2_isoform = chr_1.genes[0].isoforms[1]
        skip_2_align = tests.base_test.Alignment()
        skip_2_align.chr_name = chr_1.name
        skip_2_align.match_to_isoform(skip_2_isoform)

        skip_3_isoform = chr_1.genes[0].isoforms[2]
        skip_3_align = tests.base_test.Alignment()
        skip_3_align.chr_name = chr_1.name
        skip_3_align.match_to_isoform(skip_3_isoform)

        # gene 2
        g_2_all_isoform = chr_1.genes[1].isoforms[0]
        g_2_all_align = tests.base_test.Alignment()
        g_2_all_align.chr_name = chr_1.name
        g_2_all_align.match_to_isoform(g_2_all_isoform)

        g_2_skip_2_isoform = chr_1.genes[1].isoforms[1]
        g_2_skip_2_align = tests.base_test.Alignment()
        g_2_skip_2_align.chr_name = chr_1.name
        g_2_skip_2_align.match_to_isoform(g_2_skip_2_isoform)

        g_2_skip_3_isoform = chr_1.genes[1].isoforms[2]
        g_2_skip_3_align = tests.base_test.Alignment()
        g_2_skip_3_align.chr_name = chr_1.name
        g_2_skip_3_align.match_to_isoform(g_2_skip_3_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(all_align, 10, aligns)
        tests.base_test.append_copies(skip_2_align, 10, aligns)
        tests.base_test.append_copies(skip_3_align, 10, aligns)
        tests.base_test.append_copies(g_2_all_align, 100, aligns)
        tests.base_test.append_copies(g_2_skip_2_align, 100, aligns)
        tests.base_test.append_copies(g_2_skip_3_align, 100, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 11, aligns)
        tests.base_test.append_copies(skip_2_align, 10, aligns)
        tests.base_test.append_copies(skip_3_align, 9, aligns)
        tests.base_test.append_copies(g_2_all_align, 110, aligns)
        tests.base_test.append_copies(g_2_skip_2_align, 100, aligns)
        tests.base_test.append_copies(g_2_skip_3_align, 90, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 6, aligns)
        tests.base_test.append_copies(skip_2_align, 7, aligns)
        tests.base_test.append_copies(skip_3_align, 17, aligns)
        tests.base_test.append_copies(g_2_all_align, 60, aligns)
        tests.base_test.append_copies(g_2_skip_2_align, 70, aligns)
        tests.base_test.append_copies(g_2_skip_3_align, 170, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 5, aligns)
        tests.base_test.append_copies(skip_2_align, 6, aligns)
        tests.base_test.append_copies(skip_3_align, 19, aligns)
        tests.base_test.append_copies(g_2_all_align, 50, aligns)
        tests.base_test.append_copies(g_2_skip_2_align, 60, aligns)
        tests.base_test.append_copies(g_2_skip_3_align, 190, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    def _check_output(self):
        unfilt_path = os.path.join(self._rmats_long_out,
                                   'differential_isoforms.tsv')
        filt_path = os.path.join(self._rmats_long_out,
                                 'differential_isoforms_filtered.tsv')
        chr_id = '0'
        event_i = 0
        g_1_asm_id = '{}_{}'.format(chr_id, event_i)
        g_1_all_isoform_id = '{}_0'.format(g_1_asm_id)
        g_1_skip_3_isoform_id = '{}_1'.format(g_1_asm_id)
        g_1_skip_2_isoform_id = '{}_2'.format(g_1_asm_id)

        event_i += 1
        # g_1_skip_2_asm_id = '{}_{}'.format(chr_id, event_i)
        # g_1_skip_2_inc_isoform_id = '{}_0'.format(g_1_skip_2_asm_id)
        # g_1_skip_2_skip_isoform_id = '{}_1'.format(g_1_skip_2_asm_id)

        event_i += 1
        # g_1_skip_3_asm_id = '{}_{}'.format(chr_id, event_i)
        # g_1_skip_3_inc_isoform_id = '{}_0'.format(g_1_skip_3_asm_id)
        # g_1_skip_3_skip_isoform_id = '{}_1'.format(g_1_skip_3_asm_id)

        event_i += 1
        g_2_asm_id = '{}_{}'.format(chr_id, event_i)
        g_2_all_isoform_id = '{}_0'.format(g_2_asm_id)
        g_2_skip_3_isoform_id = '{}_1'.format(g_2_asm_id)
        g_2_skip_2_isoform_id = '{}_2'.format(g_2_asm_id)

        event_i += 1
        # g_2_skip_2_asm_id = '{}_{}'.format(chr_id, event_i)
        # g_2_skip_2_inc_isoform_id = '{}_0'.format(g_2_skip_2_asm_id)
        # g_2_skip_2_skip_isoform_id = '{}_1'.format(g_2_skip_2_asm_id)

        event_i += 1
        # g_2_skip_3_asm_id = '{}_{}'.format(chr_id, event_i)
        # g_2_skip_3_inc_isoform_id = '{}_0'.format(g_2_skip_3_asm_id)
        # g_2_skip_3_skip_isoform_id = '{}_1'.format(g_2_skip_3_asm_id)

        unfilt_transcripts = self.parse_differential_transcripts(unfilt_path,
                                                                 is_asm=True)
        filt_transcripts = self.parse_differential_transcripts(filt_path,
                                                               is_asm=True)
        self.assertIn(g_1_all_isoform_id, unfilt_transcripts)
        self.assertIn(g_1_skip_2_isoform_id, unfilt_transcripts)
        self.assertIn(g_1_skip_3_isoform_id, unfilt_transcripts)
        self.assertIn(g_2_all_isoform_id, unfilt_transcripts)
        self.assertIn(g_2_skip_2_isoform_id, unfilt_transcripts)
        self.assertIn(g_2_skip_3_isoform_id, unfilt_transcripts)

        self.assertNotIn(g_1_skip_3_isoform_id, filt_transcripts)
        row = unfilt_transcripts[g_1_skip_3_isoform_id]
        self.assertEqual(row['asm_id'], g_1_asm_id)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.283,
                               delta=0.001)

        self.assertIn(g_2_skip_3_isoform_id, filt_transcripts)
        row = unfilt_transcripts[g_2_skip_3_isoform_id]
        self.assertEqual(row['asm_id'], g_2_asm_id)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.283,
                               delta=0.001)


class MinCpmTest(SignificanceFiltersBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'min_cpm'
        self._set_input_out_log_with_prefix(self._prefix)
        self._min_cpm_group = 100_000  # 10%

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1 has 4 exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        gene.isoforms.append(isoform)
        # transcript 2 skips exon 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        gene.isoforms.append(isoform)
        # transcript 3 has skips exon 3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        gene.isoforms.append(isoform)

        # gene 2
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1 has 4 exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 2500))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3000, 3100))
        gene.isoforms.append(isoform)
        # transcript 2 skips exon 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3000, 3100))
        gene.isoforms.append(isoform)
        # transcript 3 has skips exon 3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 2500))
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
        all_isoform = chr_1.genes[0].isoforms[0]
        all_align = tests.base_test.Alignment()
        all_align.chr_name = chr_1.name
        all_align.match_to_isoform(all_isoform)

        skip_2_isoform = chr_1.genes[0].isoforms[1]
        skip_2_align = tests.base_test.Alignment()
        skip_2_align.chr_name = chr_1.name
        skip_2_align.match_to_isoform(skip_2_isoform)

        skip_3_isoform = chr_1.genes[0].isoforms[2]
        skip_3_align = tests.base_test.Alignment()
        skip_3_align.chr_name = chr_1.name
        skip_3_align.match_to_isoform(skip_3_isoform)

        # gene 2
        g_2_all_isoform = chr_1.genes[1].isoforms[0]
        g_2_all_align = tests.base_test.Alignment()
        g_2_all_align.chr_name = chr_1.name
        g_2_all_align.match_to_isoform(g_2_all_isoform)

        g_2_skip_2_isoform = chr_1.genes[1].isoforms[1]
        g_2_skip_2_align = tests.base_test.Alignment()
        g_2_skip_2_align.chr_name = chr_1.name
        g_2_skip_2_align.match_to_isoform(g_2_skip_2_isoform)

        g_2_skip_3_isoform = chr_1.genes[1].isoforms[2]
        g_2_skip_3_align = tests.base_test.Alignment()
        g_2_skip_3_align.chr_name = chr_1.name
        g_2_skip_3_align.match_to_isoform(g_2_skip_3_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(all_align, 10, aligns)
        tests.base_test.append_copies(skip_2_align, 10, aligns)
        tests.base_test.append_copies(skip_3_align, 10, aligns)
        tests.base_test.append_copies(g_2_all_align, 200, aligns)
        tests.base_test.append_copies(g_2_skip_2_align, 200, aligns)
        tests.base_test.append_copies(g_2_skip_3_align, 200, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 90, aligns)
        tests.base_test.append_copies(skip_2_align, 100, aligns)
        tests.base_test.append_copies(skip_3_align, 110, aligns)
        tests.base_test.append_copies(g_2_all_align, 90, aligns)
        tests.base_test.append_copies(g_2_skip_2_align, 100, aligns)
        tests.base_test.append_copies(g_2_skip_3_align, 110, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 50, aligns)
        tests.base_test.append_copies(skip_2_align, 50, aligns)
        tests.base_test.append_copies(skip_3_align, 200, aligns)
        tests.base_test.append_copies(g_2_all_align, 50, aligns)
        tests.base_test.append_copies(g_2_skip_2_align, 50, aligns)
        tests.base_test.append_copies(g_2_skip_3_align, 200, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 45, aligns)
        tests.base_test.append_copies(skip_2_align, 55, aligns)
        tests.base_test.append_copies(skip_3_align, 210, aligns)
        tests.base_test.append_copies(g_2_all_align, 45, aligns)
        tests.base_test.append_copies(g_2_skip_2_align, 55, aligns)
        tests.base_test.append_copies(g_2_skip_3_align, 210, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    def _check_output(self):
        unfilt_path = os.path.join(self._rmats_long_out,
                                   'differential_isoforms.tsv')
        filt_path = os.path.join(self._rmats_long_out,
                                 'differential_isoforms_filtered.tsv')
        chr_id = '0'
        event_i = 0
        g_1_asm_id = '{}_{}'.format(chr_id, event_i)
        g_1_all_isoform_id = '{}_0'.format(g_1_asm_id)
        g_1_skip_3_isoform_id = '{}_1'.format(g_1_asm_id)
        g_1_skip_2_isoform_id = '{}_2'.format(g_1_asm_id)

        event_i += 1
        # g_1_skip_2_asm_id = '{}_{}'.format(chr_id, event_i)
        # g_1_skip_2_inc_isoform_id = '{}_0'.format(g_1_skip_2_asm_id)
        # g_1_skip_2_skip_isoform_id = '{}_1'.format(g_1_skip_2_asm_id)

        event_i += 1
        # g_1_skip_3_asm_id = '{}_{}'.format(chr_id, event_i)
        # g_1_skip_3_inc_isoform_id = '{}_0'.format(g_1_skip_3_asm_id)
        # g_1_skip_3_skip_isoform_id = '{}_1'.format(g_1_skip_3_asm_id)

        event_i += 1
        g_2_asm_id = '{}_{}'.format(chr_id, event_i)
        g_2_all_isoform_id = '{}_0'.format(g_2_asm_id)
        g_2_skip_3_isoform_id = '{}_1'.format(g_2_asm_id)
        g_2_skip_2_isoform_id = '{}_2'.format(g_2_asm_id)

        event_i += 1
        # g_2_skip_2_asm_id = '{}_{}'.format(chr_id, event_i)
        # g_2_skip_2_inc_isoform_id = '{}_0'.format(g_2_skip_2_asm_id)
        # g_2_skip_2_skip_isoform_id = '{}_1'.format(g_2_skip_2_asm_id)

        event_i += 1
        # g_2_skip_3_asm_id = '{}_{}'.format(chr_id, event_i)
        # g_2_skip_3_inc_isoform_id = '{}_0'.format(g_2_skip_3_asm_id)
        # g_2_skip_3_skip_isoform_id = '{}_1'.format(g_2_skip_3_asm_id)

        unfilt_transcripts = self.parse_differential_transcripts(unfilt_path,
                                                                 is_asm=True)
        filt_transcripts = self.parse_differential_transcripts(filt_path,
                                                               is_asm=True)
        self.assertIn(g_1_all_isoform_id, unfilt_transcripts)
        self.assertIn(g_1_skip_2_isoform_id, unfilt_transcripts)
        self.assertIn(g_1_skip_3_isoform_id, unfilt_transcripts)
        self.assertIn(g_2_all_isoform_id, unfilt_transcripts)
        self.assertIn(g_2_skip_2_isoform_id, unfilt_transcripts)
        self.assertIn(g_2_skip_3_isoform_id, unfilt_transcripts)

        self.assertNotIn(g_1_skip_3_isoform_id, filt_transcripts)
        row = unfilt_transcripts[g_1_skip_3_isoform_id]
        self.assertEqual(row['asm_id'], g_1_asm_id)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.322,
                               delta=0.001)

        self.assertIn(g_2_skip_3_isoform_id, filt_transcripts)
        row = unfilt_transcripts[g_2_skip_3_isoform_id]
        self.assertEqual(row['asm_id'], g_2_asm_id)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.322,
                               delta=0.001)


class AverageReadsTest(SignificanceFiltersBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'avg_reads'
        self._set_input_out_log_with_prefix(self._prefix)
        self._avg_reads_group = 10

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1 includes the exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(700, 800))
        gene.isoforms.append(isoform)
        # transcript 2 skips the exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(700, 800))
        gene.isoforms.append(isoform)

        # gene 2
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1 includes the exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 2500))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        gene.isoforms.append(isoform)
        # transcript 2 skips the exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
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
        inc_isoform = chr_1.genes[0].isoforms[0]
        inc_align = tests.base_test.Alignment()
        inc_align.chr_name = chr_1.name
        inc_align.match_to_isoform(inc_isoform)

        skip_isoform = chr_1.genes[0].isoforms[1]
        skip_align = tests.base_test.Alignment()
        skip_align.chr_name = chr_1.name
        skip_align.match_to_isoform(skip_isoform)

        # gene 2
        g_2_inc_isoform = chr_1.genes[1].isoforms[0]
        g_2_inc_align = tests.base_test.Alignment()
        g_2_inc_align.chr_name = chr_1.name
        g_2_inc_align.match_to_isoform(g_2_inc_isoform)

        g_2_skip_isoform = chr_1.genes[1].isoforms[1]
        g_2_skip_align = tests.base_test.Alignment()
        g_2_skip_align.chr_name = chr_1.name
        g_2_skip_align.match_to_isoform(g_2_skip_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(inc_align, 6, aligns)
        tests.base_test.append_copies(skip_align, 2, aligns)
        tests.base_test.append_copies(g_2_inc_align, 12, aligns)
        tests.base_test.append_copies(g_2_skip_align, 4, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 7, aligns)
        tests.base_test.append_copies(skip_align, 3, aligns)
        tests.base_test.append_copies(g_2_inc_align, 14, aligns)
        tests.base_test.append_copies(g_2_skip_align, 6, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 2, aligns)
        tests.base_test.append_copies(skip_align, 6, aligns)
        tests.base_test.append_copies(g_2_inc_align, 4, aligns)
        tests.base_test.append_copies(g_2_skip_align, 12, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 3, aligns)
        tests.base_test.append_copies(skip_align, 7, aligns)
        tests.base_test.append_copies(g_2_inc_align, 6, aligns)
        tests.base_test.append_copies(g_2_skip_align, 14, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    def _check_output(self):
        unfilt_path = os.path.join(self._rmats_long_out,
                                   'differential_isoforms.tsv')
        filt_path = os.path.join(self._rmats_long_out,
                                 'differential_isoforms_filtered.tsv')
        chr_id = '0'
        event_i = 0
        g_1_asm_id = '{}_{}'.format(chr_id, event_i)
        g_1_inc_isoform_id = '{}_0'.format(g_1_asm_id)
        g_1_skip_isoform_id = '{}_1'.format(g_1_asm_id)

        event_i += 1
        g_2_asm_id = '{}_{}'.format(chr_id, event_i)
        g_2_inc_isoform_id = '{}_0'.format(g_2_asm_id)
        g_2_skip_isoform_id = '{}_1'.format(g_2_asm_id)

        unfilt_transcripts = self.parse_differential_transcripts(unfilt_path,
                                                                 is_asm=True)
        filt_transcripts = self.parse_differential_transcripts(filt_path,
                                                               is_asm=True)
        self.assertIn(g_1_inc_isoform_id, unfilt_transcripts)
        self.assertIn(g_1_skip_isoform_id, unfilt_transcripts)
        self.assertIn(g_2_inc_isoform_id, unfilt_transcripts)
        self.assertIn(g_2_skip_isoform_id, unfilt_transcripts)

        self.assertNotIn(g_1_skip_isoform_id, filt_transcripts)
        row = unfilt_transcripts[g_1_skip_isoform_id]
        self.assertEqual(row['asm_id'], g_1_asm_id)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.45,
                               delta=0.001)

        self.assertIn(g_2_skip_isoform_id, filt_transcripts)
        row = unfilt_transcripts[g_2_skip_isoform_id]
        self.assertEqual(row['asm_id'], g_2_asm_id)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']),
                               -0.45,
                               delta=0.001)


if __name__ == '__main__':
    unittest.main(verbosity=2)
