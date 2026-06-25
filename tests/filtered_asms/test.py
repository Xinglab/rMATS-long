import os
import os.path
import unittest

import tests.base_test


class FilteredAsmsTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_dir = os.path.dirname(__file__)
        self._input_dir = os.path.join(self._test_dir, 'input')
        self._out_dir = os.path.join(self._test_dir, 'out')
        self._log_dir = os.path.join(self._test_dir, 'log')

    def test(self):
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
                                                 align_dir=self._align_dir,
                                                 min_reads=100)
        self._asm_counts = self.run_count_reads(self._gtf_dir, self._align_dir,
                                                self._event_dir, self._out_dir,
                                                self._log_dir)
        self._rmats_long_out = self.run_rmats_long(
            self._group_1, self._group_2, self._align_dir, self._event_dir,
            self._asm_counts, self._out_dir, self._log_dir)
        self._check_output()

    def _initialize_dirs(self):
        dirs = [self._input_dir, self._out_dir, self._log_dir]
        tests.base_test.remove_then_create_directories(dirs)

    def _create_test_inputs(self):
        self._create_gtf()
        self._create_alignments()
        self._group_1 = tests.base_test.create_group_1_file(
            ['sample_1_1', 'sample_1_2'], self._input_dir)
        self._group_2 = tests.base_test.create_group_2_file(
            ['sample_2_1', 'sample_2_2'], self._input_dir)

    # The full splice graph allows these isoforms:
    # 1,2,3,4,5,6,7
    # 1,2,3,4,5,7
    # 1,2,3,4,6,7
    # 1,2,3,5,6,7
    # 1,2,3,5,7
    # 1,2,4,5,6,7
    # 1,2,4,5,7
    # 1,2,4,6,7
    # 1,3,4,5,6,7
    # 1,3,4,5,7
    # 1,3,4,6,7
    # 1,3,5,6,7
    # 1,3,5,7
    # There are 20 reads for each isoform. The junctions have this many reads:
    # 1,2: 160
    # 1,3: 100
    # 2,3: 100
    # 2,4: 60
    # 3,4: 120
    # 3,5: 80
    # 4,5: 120
    # 4,6: 60
    # 5,6: 100
    # 5,7: 100
    # 6,7: 160
    # Isoforms with 2,4 or 3,5 or 4,6 will not be detectable in the filtered splice graph.
    # The remaining isoforms are
    # 1,2,3,4,5,6,7
    # 1,2,3,4,5,7
    # 1,3,4,5,6,7
    # 1,3,4,5,7
    # The events are:
    # 1,2,3;1,3
    # 5,6,7;5,7
    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1: 1,2,3,4,5,6,7
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

        # transcript 2: 1,3,5,7
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        isoform.exons.append(tests.base_test.Region(1900, 2000))
        gene.isoforms.append(isoform)

        # transcript 3: 1,2,4,6,7
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        isoform.exons.append(tests.base_test.Region(1600, 1700))
        isoform.exons.append(tests.base_test.Region(1900, 2000))
        gene.isoforms.append(isoform)

        # gene_2 has three isoforms.
        # After SJs are filtered, a source to sink ASM has two isoforms.
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1: 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(3100, 3200))
        gene.isoforms.append(isoform)

        # transcript 2: 2,3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(3400, 3500))
        isoform.exons.append(tests.base_test.Region(3700, 3800))
        gene.isoforms.append(isoform)

        # transcript 3: 4,5,6
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(4000, 4100))
        isoform.exons.append(tests.base_test.Region(4300, 4400))
        isoform.exons.append(tests.base_test.Region(4600, 4700))
        gene.isoforms.append(isoform)

        # gene_3 has three isoforms:
        # 1,2,             4,5
        # 1,2(early start),4,5
        # 1,3,             4,5
        # The 4-5 SJ will be filtered
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_3'

        # transcript 1: 1,2,4,5
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(6100, 6200))
        isoform.exons.append(tests.base_test.Region(6500, 6600))
        isoform.exons.append(tests.base_test.Region(7100, 7200))
        isoform.exons.append(tests.base_test.Region(7400, 7500))
        gene.isoforms.append(isoform)

        # transcript 2: 1,2(early start),4,5
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(6100, 6200))
        isoform.exons.append(tests.base_test.Region(6400, 6600))
        isoform.exons.append(tests.base_test.Region(7100, 7200))
        isoform.exons.append(tests.base_test.Region(7400, 7500))
        gene.isoforms.append(isoform)

        # transcript 3: 1,3,4,5
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(6100, 6200))
        isoform.exons.append(tests.base_test.Region(6800, 6900))
        isoform.exons.append(tests.base_test.Region(7100, 7200))
        isoform.exons.append(tests.base_test.Region(7400, 7500))
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
        iso_1234567 = chr_1.genes[0].isoforms[0]
        align_1234567 = tests.base_test.Alignment()
        align_1234567.chr_name = chr_1.name
        align_1234567.match_to_isoform(iso_1234567)

        iso_123457 = iso_1234567.copy()
        iso_123457.exons.pop(5)
        align_123457 = tests.base_test.Alignment()
        align_123457.chr_name = chr_1.name
        align_123457.match_to_isoform(iso_123457)

        iso_123467 = iso_1234567.copy()
        iso_123467.exons.pop(4)
        align_123467 = tests.base_test.Alignment()
        align_123467.chr_name = chr_1.name
        align_123467.match_to_isoform(iso_123467)

        iso_123567 = iso_1234567.copy()
        iso_123567.exons.pop(3)
        align_123567 = tests.base_test.Alignment()
        align_123567.chr_name = chr_1.name
        align_123567.match_to_isoform(iso_123567)

        iso_124567 = iso_1234567.copy()
        iso_124567.exons.pop(2)
        align_124567 = tests.base_test.Alignment()
        align_124567.chr_name = chr_1.name
        align_124567.match_to_isoform(iso_124567)

        iso_134567 = iso_1234567.copy()
        iso_134567.exons.pop(1)
        align_134567 = tests.base_test.Alignment()
        align_134567.chr_name = chr_1.name
        align_134567.match_to_isoform(iso_134567)

        iso_12357 = iso_123457.copy()
        iso_12357.exons.pop(3)
        align_12357 = tests.base_test.Alignment()
        align_12357.chr_name = chr_1.name
        align_12357.match_to_isoform(iso_12357)

        iso_12457 = iso_123457.copy()
        iso_12457.exons.pop(2)
        align_12457 = tests.base_test.Alignment()
        align_12457.chr_name = chr_1.name
        align_12457.match_to_isoform(iso_12457)

        iso_12467 = iso_123467.copy()
        iso_12467.exons.pop(2)
        align_12467 = tests.base_test.Alignment()
        align_12467.chr_name = chr_1.name
        align_12467.match_to_isoform(iso_12467)

        iso_13457 = iso_123457.copy()
        iso_13457.exons.pop(1)
        align_13457 = tests.base_test.Alignment()
        align_13457.chr_name = chr_1.name
        align_13457.match_to_isoform(iso_13457)

        iso_13467 = iso_123467.copy()
        iso_13467.exons.pop(1)
        align_13467 = tests.base_test.Alignment()
        align_13467.chr_name = chr_1.name
        align_13467.match_to_isoform(iso_13467)

        iso_13567 = iso_123567.copy()
        iso_13567.exons.pop(1)
        align_13567 = tests.base_test.Alignment()
        align_13567.chr_name = chr_1.name
        align_13567.match_to_isoform(iso_13567)

        iso_1357 = iso_12357.copy()
        iso_1357.exons.pop(1)
        align_1357 = tests.base_test.Alignment()
        align_1357.chr_name = chr_1.name
        align_1357.match_to_isoform(iso_1357)

        # gene_2
        g2_iso_1 = chr_1.genes[1].isoforms[0]
        g2_align_1 = tests.base_test.Alignment()
        g2_align_1.chr_name = chr_1.name
        g2_align_1.match_to_isoform(g2_iso_1)

        g2_iso_23 = chr_1.genes[1].isoforms[1]
        g2_align_23 = tests.base_test.Alignment()
        g2_align_23.chr_name = chr_1.name
        g2_align_23.match_to_isoform(g2_iso_23)

        g2_iso_456 = chr_1.genes[1].isoforms[2]
        g2_align_456 = tests.base_test.Alignment()
        g2_align_456.chr_name = chr_1.name
        g2_align_456.match_to_isoform(g2_iso_456)

        # gene_3
        g3_iso_1 = chr_1.genes[2].isoforms[0]
        g3_iso_1_no_5 = g3_iso_1.copy()
        g3_iso_1_no_5.exons.pop()
        g3_align_1_no_5 = tests.base_test.Alignment()
        g3_align_1_no_5.chr_name = chr_1.name
        g3_align_1_no_5.match_to_isoform(g3_iso_1_no_5)

        g3_iso_2 = chr_1.genes[2].isoforms[1]
        g3_iso_2_no_5 = g3_iso_2.copy()
        g3_iso_2_no_5.exons.pop()
        g3_align_2_no_5 = tests.base_test.Alignment()
        g3_align_2_no_5.chr_name = chr_1.name
        g3_align_2_no_5.match_to_isoform(g3_iso_2_no_5)

        g3_iso_3 = chr_1.genes[2].isoforms[2]
        g3_iso_3_no_5 = g3_iso_3.copy()
        g3_iso_3_no_5.exons.pop()
        g3_align_3_no_5 = tests.base_test.Alignment()
        g3_align_3_no_5.chr_name = chr_1.name
        g3_align_3_no_5.match_to_isoform(g3_iso_3_no_5)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(align_1234567, 10, aligns)
        tests.base_test.append_copies(align_123467, 10, aligns)
        tests.base_test.append_copies(align_12357, 10, aligns)
        tests.base_test.append_copies(align_12457, 10, aligns)
        tests.base_test.append_copies(align_134567, 10, aligns)
        tests.base_test.append_copies(align_13467, 10, aligns)
        tests.base_test.append_copies(align_1357, 10, aligns)
        tests.base_test.append_copies(g2_align_1, 25, aligns)
        tests.base_test.append_copies(g2_align_23, 25, aligns)
        tests.base_test.append_copies(g2_align_456, 1, aligns)
        tests.base_test.append_copies(g3_align_1_no_5, 25, aligns)
        tests.base_test.append_copies(g3_align_2_no_5, 25, aligns)
        tests.base_test.append_copies(g3_align_3_no_5, 25, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(align_123457, 10, aligns)
        tests.base_test.append_copies(align_123567, 10, aligns)
        tests.base_test.append_copies(align_124567, 10, aligns)
        tests.base_test.append_copies(align_12467, 10, aligns)
        tests.base_test.append_copies(align_13457, 10, aligns)
        tests.base_test.append_copies(align_13567, 10, aligns)
        tests.base_test.append_copies(g2_align_1, 25, aligns)
        tests.base_test.append_copies(g2_align_23, 25, aligns)
        tests.base_test.append_copies(g2_align_456, 1, aligns)
        tests.base_test.append_copies(g3_align_1_no_5, 25, aligns)
        tests.base_test.append_copies(g3_align_2_no_5, 25, aligns)
        tests.base_test.append_copies(g3_align_3_no_5, 25, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    def _check_output(self):
        chr_id = '0'
        skip_2_event_i = '0'
        skip_2_id = '{}_{}'.format(chr_id, skip_2_event_i)
        skip_6_event_i = '1'
        skip_6_id = '{}_{}'.format(chr_id, skip_6_event_i)
        inc_2_iso_id = '{}_0'.format(skip_2_id)
        skip_2_iso_id = '{}_1'.format(skip_2_id)
        inc_6_iso_id = '{}_0'.format(skip_6_id)
        skip_6_iso_id = '{}_1'.format(skip_6_id)
        g2_event_i = '2'
        g2_id = '{}_{}'.format(chr_id, g2_event_i)
        g2_iso_1_id = '{}_0'.format(g2_id)
        g2_iso_2_id = '{}_1'.format(g2_id)
        g3_alt_2_event_i = '3'
        g3_alt_2_id = '{}_{}'.format(chr_id, g3_alt_2_event_i)
        g3_alt_2_long_id = '{}_0'.format(g3_alt_2_id)
        g3_alt_2_short_id = '{}_1'.format(g3_alt_2_id)
        g3_full_event_i = '4'
        g3_full_id = '{}_{}'.format(chr_id, g3_full_event_i)
        g3_full_long_id = '{}_0'.format(g3_full_id)
        g3_full_short_id = '{}_1'.format(g3_full_id)
        g3_full_3_id = '{}_2'.format(g3_full_id)
        count_path = os.path.join(self._rmats_long_out, 'count.tsv')
        counts = self.parse_asm_counts(count_path)
        self.assertEqual(sorted(counts.keys()), [
            inc_2_iso_id, skip_2_iso_id, inc_6_iso_id, skip_6_iso_id,
            g2_iso_1_id, g2_iso_2_id, g3_alt_2_long_id, g3_alt_2_short_id,
            g3_full_long_id, g3_full_short_id, g3_full_3_id
        ])
        self.assertEqual(counts[inc_2_iso_id]['sample_1_1'], 30)
        self.assertEqual(counts[inc_2_iso_id]['sample_1_2'], 30)
        self.assertEqual(counts[inc_2_iso_id]['sample_2_1'], 20)
        self.assertEqual(counts[inc_2_iso_id]['sample_2_2'], 20)
        self.assertEqual(counts[skip_2_iso_id]['sample_1_1'], 30)
        self.assertEqual(counts[skip_2_iso_id]['sample_1_2'], 30)
        self.assertEqual(counts[skip_2_iso_id]['sample_2_1'], 20)
        self.assertEqual(counts[skip_2_iso_id]['sample_2_2'], 20)
        self.assertEqual(counts[inc_6_iso_id]['sample_1_1'], 20)
        self.assertEqual(counts[inc_6_iso_id]['sample_1_2'], 20)
        self.assertEqual(counts[inc_6_iso_id]['sample_2_1'], 30)
        self.assertEqual(counts[inc_6_iso_id]['sample_2_2'], 30)
        self.assertEqual(counts[skip_6_iso_id]['sample_1_1'], 30)
        self.assertEqual(counts[skip_6_iso_id]['sample_1_2'], 30)
        self.assertEqual(counts[skip_6_iso_id]['sample_2_1'], 20)
        self.assertEqual(counts[skip_6_iso_id]['sample_2_2'], 20)
        self.assertEqual(counts[g2_iso_1_id]['sample_1_1'], 25)
        self.assertEqual(counts[g2_iso_1_id]['sample_1_2'], 25)
        self.assertEqual(counts[g2_iso_1_id]['sample_2_1'], 25)
        self.assertEqual(counts[g2_iso_1_id]['sample_2_2'], 25)
        self.assertEqual(counts[g2_iso_2_id]['sample_1_1'], 25)
        self.assertEqual(counts[g2_iso_2_id]['sample_1_2'], 25)
        self.assertEqual(counts[g2_iso_2_id]['sample_2_1'], 25)
        self.assertEqual(counts[g2_iso_2_id]['sample_2_2'], 25)
        self.assertEqual(counts[g3_alt_2_long_id]['sample_1_1'], 25)
        self.assertEqual(counts[g3_alt_2_long_id]['sample_1_2'], 25)
        self.assertEqual(counts[g3_alt_2_long_id]['sample_2_1'], 25)
        self.assertEqual(counts[g3_alt_2_long_id]['sample_2_2'], 25)
        self.assertEqual(counts[g3_alt_2_short_id]['sample_1_1'], 25)
        self.assertEqual(counts[g3_alt_2_short_id]['sample_1_2'], 25)
        self.assertEqual(counts[g3_alt_2_short_id]['sample_2_1'], 25)
        self.assertEqual(counts[g3_alt_2_short_id]['sample_2_2'], 25)
        self.assertEqual(counts[g3_full_long_id]['sample_1_1'], 25)
        self.assertEqual(counts[g3_full_long_id]['sample_1_2'], 25)
        self.assertEqual(counts[g3_full_long_id]['sample_2_1'], 25)
        self.assertEqual(counts[g3_full_long_id]['sample_2_2'], 25)
        self.assertEqual(counts[g3_full_short_id]['sample_1_1'], 25)
        self.assertEqual(counts[g3_full_short_id]['sample_1_2'], 25)
        self.assertEqual(counts[g3_full_short_id]['sample_2_1'], 25)
        self.assertEqual(counts[g3_full_short_id]['sample_2_2'], 25)
        self.assertEqual(counts[g3_full_3_id]['sample_1_1'], 25)
        self.assertEqual(counts[g3_full_3_id]['sample_1_2'], 25)
        self.assertEqual(counts[g3_full_3_id]['sample_2_1'], 25)
        self.assertEqual(counts[g3_full_3_id]['sample_2_2'], 25)

        summary_txt_path = os.path.join(self._rmats_long_out, 'summary.txt')
        self.assert_exists(summary_txt_path)
        summary = self.parse_summary_txt(summary_txt_path)
        self.assertEqual(summary['total ASMs with significant isoforms'], '1')
        self.assertEqual(summary['exon skipping'], '1')
        filt_path = os.path.join(self._rmats_long_out,
                                 'differential_isoforms_filtered.tsv')
        self.assert_exists(filt_path)
        diff_transcripts = self.parse_differential_transcripts(filt_path,
                                                               is_asm=True)
        self.assertIn(inc_6_iso_id, diff_transcripts)
        row = diff_transcripts[inc_6_iso_id]
        self.assertEqual(row['asm_id'], skip_6_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), -0.2)
        self.assertLess(float(row['adj_pvalue']), 0.05)

        asms = self.parse_asm_definitions(self._event_dir)
        exons_1 = ((3101, 3201), )
        exons_23 = ((3401, 3501), (3701, 3801))
        key = tuple(sorted([exons_1, exons_23]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        exons_1 = ((6101, 6201), (6501, 6601))
        exons_2 = ((6101, 6201), (6401, 6601))
        key = tuple(sorted([exons_1, exons_2]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], "alternative 3'-splice site usage")
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_3')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        exons_1 = ((6101, 6201), (6501, 6601), (7101, 7201))
        exons_2 = ((6101, 6201), (6401, 6601), (7101, 7201))
        exons_3 = ((6101, 6201), (6801, 6901), (7101, 7201))
        key = tuple(sorted([exons_1, exons_2, exons_3]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_3')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')


if __name__ == '__main__':
    unittest.main(verbosity=2)
