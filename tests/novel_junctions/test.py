import os
import os.path
import unittest

import tests.base_test


class NovelJunctionsTest(tests.base_test.BaseTest):
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
                                                 novel_junctions=True)
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

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1: 1,2,3,4
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1: 1,2,3,4
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 2500))
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
        iso_all = chr_1.genes[0].isoforms[0]
        align_all = tests.base_test.Alignment()
        align_all.chr_name = chr_1.name
        align_all.match_to_isoform(iso_all)

        iso_skip_2 = iso_all.copy()
        iso_skip_2.exons.pop(1)
        align_skip_2 = tests.base_test.Alignment()
        align_skip_2.chr_name = chr_1.name
        align_skip_2.match_to_isoform(iso_skip_2)

        iso_skip_23 = iso_skip_2.copy()
        iso_skip_23.exons.pop(1)
        align_skip_23 = tests.base_test.Alignment()
        align_skip_23.chr_name = chr_1.name
        align_skip_23.match_to_isoform(iso_skip_23)

        iso_skip_3 = iso_all.copy()
        iso_skip_3.exons.pop(2)
        align_skip_3 = tests.base_test.Alignment()
        align_skip_3.chr_name = chr_1.name
        align_skip_3.match_to_isoform(iso_skip_3)

        iso_ri_3 = iso_all.copy()
        iso_ri_3.exons[2].end = iso_ri_3.exons[3].end
        iso_ri_3.exons.pop(3)
        align_ri_3 = tests.base_test.Alignment()
        align_ri_3.chr_name = chr_1.name
        align_ri_3.match_to_isoform(iso_ri_3)

        # gene_2
        g2_iso_all = chr_1.genes[1].isoforms[0]
        g2_align_all = tests.base_test.Alignment()
        g2_align_all.chr_name = chr_1.name
        g2_align_all.match_to_isoform(g2_iso_all)

        g2_iso_skip_2 = g2_iso_all.copy()
        g2_iso_skip_2.exons.pop(1)
        g2_align_skip_2 = tests.base_test.Alignment()
        g2_align_skip_2.chr_name = chr_1.name
        g2_align_skip_2.match_to_isoform(g2_iso_skip_2)

        g2_iso_skip_23 = g2_iso_skip_2.copy()
        g2_iso_skip_23.exons.pop(1)
        g2_align_skip_23 = tests.base_test.Alignment()
        g2_align_skip_23.chr_name = chr_1.name
        g2_align_skip_23.match_to_isoform(g2_iso_skip_23)

        g2_iso_skip_3 = g2_iso_all.copy()
        g2_iso_skip_3.exons.pop(2)
        g2_align_skip_3 = tests.base_test.Alignment()
        g2_align_skip_3.chr_name = chr_1.name
        g2_align_skip_3.match_to_isoform(g2_iso_skip_3)

        g2_iso_ri_3 = g2_iso_all.copy()
        g2_iso_ri_3.exons[2].end = g2_iso_ri_3.exons[3].end
        g2_iso_ri_3.exons.pop(3)
        g2_align_ri_3 = tests.base_test.Alignment()
        g2_align_ri_3.chr_name = chr_1.name
        g2_align_ri_3.match_to_isoform(g2_iso_ri_3)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(align_all, 8, aligns)
        tests.base_test.append_copies(align_skip_2, 11, aligns)
        tests.base_test.append_copies(align_skip_23, 12, aligns)
        tests.base_test.append_copies(align_skip_3, 1, aligns)
        tests.base_test.append_copies(align_ri_3, 2, aligns)
        tests.base_test.append_copies(g2_align_all, 8, aligns)
        tests.base_test.append_copies(g2_align_skip_2, 11, aligns)
        tests.base_test.append_copies(g2_align_skip_23, 12, aligns)
        tests.base_test.append_copies(g2_align_skip_3, 1, aligns)
        tests.base_test.append_copies(g2_align_ri_3, 2, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(align_all, 3, aligns)
        tests.base_test.append_copies(align_skip_2, 15, aligns)
        tests.base_test.append_copies(align_skip_23, 30, aligns)
        tests.base_test.append_copies(align_skip_3, 1, aligns)
        tests.base_test.append_copies(align_ri_3, 2, aligns)
        tests.base_test.append_copies(g2_align_all, 3, aligns)
        tests.base_test.append_copies(g2_align_skip_2, 15, aligns)
        tests.base_test.append_copies(g2_align_skip_23, 30, aligns)
        tests.base_test.append_copies(g2_align_skip_3, 1, aligns)
        tests.base_test.append_copies(g2_align_ri_3, 2, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    def _check_output(self):
        chr_id = '0'
        skip_23_event_i = '0'
        skip_23_id = '{}_{}'.format(chr_id, skip_23_event_i)
        skip_23_inc_iso_id = '{}_0'.format(skip_23_id)
        skip_23_skip_2_iso_id = '{}_1'.format(skip_23_id)
        skip_23_skip_23_iso_id = '{}_2'.format(skip_23_id)

        skip_2_event_i = '1'
        skip_2_id = '{}_{}'.format(chr_id, skip_2_event_i)
        skip_2_inc_iso_id = '{}_0'.format(skip_2_id)
        skip_2_skip_iso_id = '{}_1'.format(skip_2_id)

        g2_skip_23_event_i = '2'
        g2_skip_23_id = '{}_{}'.format(chr_id, g2_skip_23_event_i)
        g2_skip_23_inc_iso_id = '{}_0'.format(g2_skip_23_id)
        g2_skip_23_skip_2_iso_id = '{}_1'.format(g2_skip_23_id)
        g2_skip_23_skip_23_iso_id = '{}_2'.format(g2_skip_23_id)

        g2_skip_2_event_i = '3'
        g2_skip_2_id = '{}_{}'.format(chr_id, g2_skip_2_event_i)
        g2_skip_2_inc_iso_id = '{}_0'.format(g2_skip_2_id)
        g2_skip_2_skip_iso_id = '{}_1'.format(g2_skip_2_id)

        count_path = os.path.join(self._rmats_long_out, 'count.tsv')
        counts = self.parse_asm_counts(count_path)
        self.assertEqual(sorted(counts.keys()), [
            skip_23_inc_iso_id,
            skip_23_skip_2_iso_id,
            skip_23_skip_23_iso_id,
            skip_2_inc_iso_id,
            skip_2_skip_iso_id,
            g2_skip_23_inc_iso_id,
            g2_skip_23_skip_2_iso_id,
            g2_skip_23_skip_23_iso_id,
            g2_skip_2_inc_iso_id,
            g2_skip_2_skip_iso_id,
        ])
        self.assertEqual(counts[skip_2_inc_iso_id]['sample_1_1'], 10)
        self.assertEqual(counts[skip_2_inc_iso_id]['sample_1_2'], 10)
        self.assertEqual(counts[skip_2_inc_iso_id]['sample_2_1'], 5)
        self.assertEqual(counts[skip_2_inc_iso_id]['sample_2_2'], 5)
        self.assertEqual(counts[skip_2_skip_iso_id]['sample_1_1'], 11)
        self.assertEqual(counts[skip_2_skip_iso_id]['sample_1_2'], 11)
        self.assertEqual(counts[skip_2_skip_iso_id]['sample_2_1'], 15)
        self.assertEqual(counts[skip_2_skip_iso_id]['sample_2_2'], 15)
        self.assertEqual(counts[skip_23_inc_iso_id]['sample_1_1'], 8)
        self.assertEqual(counts[skip_23_inc_iso_id]['sample_1_2'], 8)
        self.assertEqual(counts[skip_23_inc_iso_id]['sample_2_1'], 3)
        self.assertEqual(counts[skip_23_inc_iso_id]['sample_2_2'], 3)
        self.assertEqual(counts[skip_23_skip_2_iso_id]['sample_1_1'], 11)
        self.assertEqual(counts[skip_23_skip_2_iso_id]['sample_1_2'], 11)
        self.assertEqual(counts[skip_23_skip_2_iso_id]['sample_2_1'], 15)
        self.assertEqual(counts[skip_23_skip_2_iso_id]['sample_2_2'], 15)
        self.assertEqual(counts[skip_23_skip_23_iso_id]['sample_1_1'], 12)
        self.assertEqual(counts[skip_23_skip_23_iso_id]['sample_1_2'], 12)
        self.assertEqual(counts[skip_23_skip_23_iso_id]['sample_2_1'], 30)
        self.assertEqual(counts[skip_23_skip_23_iso_id]['sample_2_2'], 30)

        summary_txt_path = os.path.join(self._rmats_long_out, 'summary.txt')
        self.assert_exists(summary_txt_path)
        summary = self.parse_summary_txt(summary_txt_path)
        self.assertEqual(summary['total ASMs with significant isoforms'], '4')
        self.assertEqual(summary['exon skipping'], '2')
        self.assertEqual(summary['complex'], '2')
        filt_path = os.path.join(self._rmats_long_out,
                                 'differential_isoforms_filtered.tsv')
        self.assert_exists(filt_path)
        diff_transcripts = self.parse_differential_transcripts(filt_path,
                                                               is_asm=True)
        self.assertIn(skip_2_skip_iso_id, diff_transcripts)
        row = diff_transcripts[skip_2_skip_iso_id]
        self.assertEqual(row['asm_id'], skip_2_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), -0.2262)
        self.assertLess(float(row['adj_pvalue']), 0.05)

        self.assertIn(skip_23_skip_23_iso_id, diff_transcripts)
        row = diff_transcripts[skip_23_skip_23_iso_id]
        self.assertEqual(row['asm_id'], skip_23_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), -0.2379)
        self.assertLess(float(row['adj_pvalue']), 0.05)

        self.assertIn(g2_skip_2_skip_iso_id, diff_transcripts)
        row = diff_transcripts[g2_skip_2_skip_iso_id]
        self.assertEqual(row['asm_id'], g2_skip_2_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), -0.2262)
        self.assertLess(float(row['adj_pvalue']), 0.05)

        self.assertIn(g2_skip_23_skip_23_iso_id, diff_transcripts)
        row = diff_transcripts[g2_skip_23_skip_23_iso_id]
        self.assertEqual(row['asm_id'], g2_skip_23_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), -0.2379)
        self.assertLess(float(row['adj_pvalue']), 0.05)


if __name__ == '__main__':
    unittest.main(verbosity=2)
