import os
import os.path
import unittest

import tests.base_test


class SkippedExonAsmTest(tests.base_test.BaseTest):
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
        self._event_dir = self.run_detect_events(self._gtf_dir, self._out_dir,
                                                 self._log_dir)
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

        # transcript 1 includes the exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        gene.isoforms.append(isoform)
        # transcript 2 skips the exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        gene.isoforms.append(isoform)

        gene_2 = tests.base_test.Gene()
        chrom.genes.append(gene_2)
        gene_2.id = 'gene_2'

        # transcript 1 includes the exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene_2.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(10100, 10200))
        isoform.exons.append(tests.base_test.Region(10400, 10500))
        isoform.exons.append(tests.base_test.Region(10700, 10800))
        isoform.exons.append(tests.base_test.Region(11000, 11100))
        isoform.exons.append(tests.base_test.Region(11300, 11400))
        gene_2.isoforms.append(isoform)
        # transcript 2 skips the exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene_2.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(10100, 10200))
        isoform.exons.append(tests.base_test.Region(10400, 10500))
        isoform.exons.append(tests.base_test.Region(11000, 11100))
        isoform.exons.append(tests.base_test.Region(11300, 11400))
        gene_2.isoforms.append(isoform)

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
        inc_isoform = chr_1.genes[0].isoforms[0]
        # 1,2,4,5
        skip_isoform = chr_1.genes[0].isoforms[1]
        inc_align = tests.base_test.Alignment()
        inc_align.chr_name = chr_1.name
        inc_align.match_to_isoform(inc_isoform)
        skip_align = tests.base_test.Alignment()
        skip_align.chr_name = chr_1.name
        skip_align.match_to_isoform(skip_isoform)

        shortened_inc = inc_isoform.copy()
        # 1,2,3,4
        shortened_inc.exons.pop()
        short_1_inc_align = tests.base_test.Alignment()
        short_1_inc_align.chr_name = chr_1.name
        short_1_inc_align.match_to_isoform(shortened_inc)
        # 1,2,3
        shortened_inc.exons.pop()
        short_2_inc_align = tests.base_test.Alignment()
        short_2_inc_align.chr_name = chr_1.name
        short_2_inc_align.match_to_isoform(shortened_inc)

        shortened_skip = skip_isoform.copy()
        # 1,2,4
        shortened_skip.exons.pop()
        short_1_skip_align = tests.base_test.Alignment()
        short_1_skip_align.chr_name = chr_1.name
        short_1_skip_align.match_to_isoform(shortened_skip)

        # gene_2
        # 1,2,3,4,5
        gene_2_inc_isoform = chr_1.genes[1].isoforms[0]
        gene_2_inc_align = tests.base_test.Alignment()
        gene_2_inc_align.chr_name = chr_1.name
        gene_2_inc_align.match_to_isoform(gene_2_inc_isoform)
        # 1,2,4,5
        gene_2_skip_isoform = chr_1.genes[1].isoforms[1]
        gene_2_skip_align = tests.base_test.Alignment()
        gene_2_skip_align.chr_name = chr_1.name
        gene_2_skip_align.match_to_isoform(gene_2_skip_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(inc_align, 13, aligns)
        tests.base_test.append_copies(short_1_inc_align, 9, aligns)
        tests.base_test.append_copies(short_2_inc_align, 8, aligns)
        tests.base_test.append_copies(skip_align, 7, aligns)
        tests.base_test.append_copies(short_1_skip_align, 3, aligns)
        tests.base_test.append_copies(gene_2_inc_align, 10, aligns)
        tests.base_test.append_copies(gene_2_skip_align, 0, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 14, aligns)
        tests.base_test.append_copies(short_1_inc_align, 11, aligns)
        tests.base_test.append_copies(short_2_inc_align, 7, aligns)
        tests.base_test.append_copies(skip_align, 6, aligns)
        tests.base_test.append_copies(short_1_skip_align, 3, aligns)
        tests.base_test.append_copies(gene_2_inc_align, 15, aligns)
        tests.base_test.append_copies(gene_2_skip_align, 0, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 9, aligns)
        tests.base_test.append_copies(short_1_inc_align, 7, aligns)
        tests.base_test.append_copies(short_2_inc_align, 4, aligns)
        tests.base_test.append_copies(skip_align, 11, aligns)
        tests.base_test.append_copies(short_1_skip_align, 9, aligns)
        tests.base_test.append_copies(gene_2_inc_align, 0, aligns)
        tests.base_test.append_copies(gene_2_skip_align, 10, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 8, aligns)
        tests.base_test.append_copies(short_1_inc_align, 5, aligns)
        tests.base_test.append_copies(short_2_inc_align, 5, aligns)
        tests.base_test.append_copies(skip_align, 11, aligns)
        tests.base_test.append_copies(short_1_skip_align, 10, aligns)
        tests.base_test.append_copies(gene_2_inc_align, 0, aligns)
        tests.base_test.append_copies(gene_2_skip_align, 15, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    def _check_output(self):
        summary_txt_path = os.path.join(self._rmats_long_out, 'summary.txt')
        summary_plot_path = os.path.join(self._rmats_long_out,
                                         'summary_plot.png')
        filt_path = os.path.join(self._rmats_long_out,
                                 'differential_isoforms_filtered.tsv')
        chr_id = '0'
        event_i = '0'
        gene_2_event_i = '1'
        gene_id = 'gene_1'
        by_gene_dir = os.path.join(self._rmats_long_out, 'results_by_gene')
        gene_dir = os.path.join(by_gene_dir, gene_id)
        asm_id = '{}_{}'.format(chr_id, event_i)
        gene_2_asm_id = '{}_{}'.format(chr_id, gene_2_event_i)
        structure_path = os.path.join(gene_dir,
                                      '{}_structure.png'.format(asm_id))
        abundance_path = os.path.join(gene_dir,
                                      '{}_abundance.png'.format(asm_id))
        inc_isoform_id = '{}_0'.format(asm_id)
        skip_isoform_id = '{}_1'.format(asm_id)
        gene_2_inc_isoform_id = '{}_0'.format(gene_2_asm_id)
        gene_2_skip_isoform_id = '{}_1'.format(gene_2_asm_id)

        count_path = os.path.join(self._rmats_long_out, 'count.tsv')
        counts = self.parse_asm_counts(count_path)
        self.assertEqual(sorted(counts.keys()), [
            inc_isoform_id, skip_isoform_id, gene_2_inc_isoform_id,
            gene_2_skip_isoform_id
        ])
        self.assertEqual(counts[inc_isoform_id]['sample_1_1'], 30)
        self.assertEqual(counts[skip_isoform_id]['sample_1_1'], 10)
        self.assertEqual(counts[gene_2_inc_isoform_id]['sample_1_1'], 10)
        self.assertEqual(counts[gene_2_skip_isoform_id].get('sample_1_1', 0),
                         0)
        self.assertEqual(counts[inc_isoform_id]['sample_1_2'], 32)
        self.assertEqual(counts[skip_isoform_id]['sample_1_2'], 9)
        self.assertEqual(counts[gene_2_inc_isoform_id]['sample_1_2'], 15)
        self.assertEqual(counts[gene_2_skip_isoform_id].get('sample_1_2', 0),
                         0)
        self.assertEqual(counts[inc_isoform_id]['sample_2_1'], 20)
        self.assertEqual(counts[skip_isoform_id]['sample_2_1'], 20)
        self.assertEqual(counts[gene_2_inc_isoform_id].get('sample_2_1', 0), 0)
        self.assertEqual(counts[gene_2_skip_isoform_id]['sample_2_1'], 10)
        self.assertEqual(counts[inc_isoform_id]['sample_2_2'], 18)
        self.assertEqual(counts[skip_isoform_id]['sample_2_2'], 21)
        self.assertEqual(counts[gene_2_inc_isoform_id].get('sample_2_2', 0), 0)
        self.assertEqual(counts[gene_2_skip_isoform_id]['sample_2_2'], 15)
        self.assert_exists(summary_txt_path)
        self.assert_exists(summary_plot_path)
        summary = self.parse_summary_txt(summary_txt_path)
        self.assertEqual(summary['total ASMs with significant isoforms'], '2')
        self.assertEqual(summary['exon skipping'], '2')
        self.assert_exists(filt_path)
        diff_transcripts = self.parse_differential_transcripts(filt_path,
                                                               is_asm=True)
        self.assertIn(inc_isoform_id, diff_transcripts)
        row = diff_transcripts[inc_isoform_id]
        self.assertEqual(row['asm_id'], asm_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), 0.2845)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertGreater(float(row['lr']), 0)
        self.assertGreater(float(row['df']), 0)
        self.assert_exists(structure_path)
        self.assert_exists(abundance_path)

        self.assertIn(gene_2_inc_isoform_id, diff_transcripts)
        row = diff_transcripts[gene_2_inc_isoform_id]
        self.assertEqual(row['asm_id'], gene_2_asm_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), 1)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertGreater(float(row['lr']), 0)
        self.assertGreater(float(row['df']), 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
