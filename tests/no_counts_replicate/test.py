import os
import os.path
import unittest

import tests.base_test


class NoCountsReplicateTest(tests.base_test.BaseTest):
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
            ['sample_1_1', 'sample_1_2', 'sample_1_3'], self._input_dir)
        self._group_2 = tests.base_test.create_group_2_file(
            ['sample_2_1', 'sample_2_2', 'sample_2_3'], self._input_dir)

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

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        self._chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, self._chrs)

    def _create_alignments(self):
        sample_1_1_sam = os.path.join(self._input_dir, 'g_1_1.sam')
        sample_1_2_sam = os.path.join(self._input_dir, 'g_1_2.sam')
        sample_1_3_sam = os.path.join(self._input_dir, 'g_1_3.sam')
        sample_2_1_sam = os.path.join(self._input_dir, 'g_2_1.sam')
        sample_2_2_sam = os.path.join(self._input_dir, 'g_2_2.sam')
        sample_2_3_sam = os.path.join(self._input_dir, 'g_2_3.sam')
        sams_by_sample = dict()
        sams_by_sample['sample_1_1'] = [sample_1_1_sam]
        sams_by_sample['sample_1_2'] = [sample_1_2_sam]
        sams_by_sample['sample_1_3'] = [sample_1_3_sam]
        sams_by_sample['sample_2_1'] = [sample_2_1_sam]
        sams_by_sample['sample_2_2'] = [sample_2_2_sam]
        sams_by_sample['sample_2_3'] = [sample_2_3_sam]
        self._sams_by_sample = sams_by_sample

        chr_1 = self._chrs[0]
        # 1,2,3
        inc_isoform = chr_1.genes[0].isoforms[0]
        inc_align = tests.base_test.Alignment()
        inc_align.chr_name = chr_1.name
        inc_align.match_to_isoform(inc_isoform)
        # 1,3
        skip_isoform = chr_1.genes[0].isoforms[1]
        skip_align = tests.base_test.Alignment()
        skip_align.chr_name = chr_1.name
        skip_align.match_to_isoform(skip_isoform)
        # 4
        other_align = tests.base_test.Alignment()
        other_align.chr_name = chr_1.name
        other_align.strand = skip_align.strand
        other_align.start = 1001
        other_align.cigar.operations = [tests.base_test.CigarOp(100, 'M')]

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(inc_align, 20, aligns)
        tests.base_test.append_copies(skip_align, 40, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 19, aligns)
        tests.base_test.append_copies(skip_align, 39, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 21, aligns)
        tests.base_test.append_copies(skip_align, 41, aligns)
        tests.base_test.write_sam(sample_1_3_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 20, aligns)
        tests.base_test.append_copies(skip_align, 20, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 19, aligns)
        tests.base_test.append_copies(skip_align, 19, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(other_align, 1, aligns)
        tests.base_test.write_sam(sample_2_3_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    def _check_output(self):
        filt_path = os.path.join(self._rmats_long_out,
                                 'differential_isoforms_filtered.tsv')
        chr_id = '0'
        event_i = '0'
        asm_id = '{}_{}'.format(chr_id, event_i)
        inc_isoform_id = '{}_0'.format(asm_id)
        skip_isoform_id = '{}_1'.format(asm_id)

        count_path = os.path.join(self._rmats_long_out, 'count.tsv')
        counts = self.parse_asm_counts(count_path)
        self.assertEqual(sorted(counts.keys()),
                         [inc_isoform_id, skip_isoform_id])
        self.assertEqual(counts[inc_isoform_id]['sample_1_1'], 20)
        self.assertEqual(counts[skip_isoform_id]['sample_1_1'], 40)
        self.assertEqual(counts[inc_isoform_id]['sample_1_2'], 19)
        self.assertEqual(counts[skip_isoform_id]['sample_1_2'], 39)
        self.assertEqual(counts[inc_isoform_id]['sample_1_3'], 21)
        self.assertEqual(counts[skip_isoform_id]['sample_1_3'], 41)
        self.assertEqual(counts[inc_isoform_id]['sample_2_1'], 20)
        self.assertEqual(counts[skip_isoform_id]['sample_2_1'], 20)
        self.assertEqual(counts[inc_isoform_id]['sample_2_2'], 19)
        self.assertEqual(counts[skip_isoform_id]['sample_2_2'], 19)
        self.assertNotIn('sample_2_3', counts[inc_isoform_id])
        self.assertNotIn('sample_2_3', counts[skip_isoform_id])
        self.assert_exists(filt_path)
        diff_transcripts = self.parse_differential_transcripts(filt_path,
                                                               is_asm=True)
        self.assertIn(inc_isoform_id, diff_transcripts)
        row = diff_transcripts[inc_isoform_id]
        self.assertEqual(row['asm_id'], asm_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), -0.1668)
        self.assertLess(float(row['adj_pvalue']), 0.05)


if __name__ == '__main__':
    unittest.main(verbosity=2)
