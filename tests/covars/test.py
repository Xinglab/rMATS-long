import os
import os.path
import unittest

import tests.base_test


class CovariatesTest(tests.base_test.BaseTest):
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
                                                   gtf=self._gtf,
                                                   covars=self._covars)

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
        self._create_covars_tsv()

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
        # transcript 2 skips exons 2 and 3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
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
        # 1,2,3,4,5
        inc_isoform = chr_1.genes[0].isoforms[0]
        # 1,2,4,5
        skip_isoform = chr_1.genes[0].isoforms[1]
        # 1,4,5
        skip_2_isoform = chr_1.genes[0].isoforms[2]
        inc_align = tests.base_test.Alignment()
        inc_align.chr_name = chr_1.name
        inc_align.match_to_isoform(inc_isoform)
        skip_align = tests.base_test.Alignment()
        skip_align.chr_name = chr_1.name
        skip_align.match_to_isoform(skip_isoform)
        skip_2_align = tests.base_test.Alignment()
        skip_2_align.chr_name = chr_1.name
        skip_2_align.match_to_isoform(skip_2_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(inc_align, 30, aligns)
        tests.base_test.append_copies(skip_align, 10, aligns)
        tests.base_test.append_copies(skip_2_align, 1, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 10, aligns)
        tests.base_test.append_copies(skip_align, 30, aligns)
        tests.base_test.append_copies(skip_2_align, 1, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 9, aligns)
        tests.base_test.append_copies(skip_align, 31, aligns)
        tests.base_test.append_copies(skip_2_align, 1, aligns)
        tests.base_test.write_sam(sample_1_3_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 10, aligns)
        tests.base_test.append_copies(skip_align, 30, aligns)
        tests.base_test.append_copies(skip_2_align, 1, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 30, aligns)
        tests.base_test.append_copies(skip_align, 10, aligns)
        tests.base_test.append_copies(skip_2_align, 1, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 31, aligns)
        tests.base_test.append_copies(skip_align, 9, aligns)
        tests.base_test.append_copies(skip_2_align, 0, aligns)
        tests.base_test.write_sam(sample_2_3_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    def _create_covars_tsv(self):
        self._covars = os.path.join(self._input_dir, 'covars.tsv')
        rows = [
            ['sample_id', 'Batch'],
            ['sample_1_1', 'First'],
            ['sample_1_2', 'Second'],
            ['sample_1_3', 'Second'],
            ['sample_2_1', 'Second'],
            ['sample_2_2', 'First'],
            ['sample_2_3', 'First'],
        ]
        with open(self._covars, 'wt') as handle:
            for row in rows:
                handle.write('{}\n'.format('\t'.join(row)))

    def _check_output(self):
        diff_path = os.path.join(self._rmats_long_out,
                                 'differential_isoforms.tsv')
        sig_path = os.path.join(self._rmats_long_out,
                                'differential_isoforms_filtered.tsv')
        chr_id = '0'
        event_i = '0'
        asm_id = '{}_{}'.format(chr_id, event_i)
        inc_isoform_id = '{}_0'.format(asm_id)
        skip_isoform_id = '{}_1'.format(asm_id)
        skip_2_isoform_id = '{}_2'.format(asm_id)
        diff_transcripts = self.parse_differential_transcripts(diff_path,
                                                               is_asm=True)
        sig_transcripts = self.parse_differential_transcripts(sig_path,
                                                              is_asm=True)
        self.assertIn(inc_isoform_id, diff_transcripts)
        self.assertIn(skip_isoform_id, diff_transcripts)
        self.assertIn(skip_2_isoform_id, diff_transcripts)
        self.assertNotIn(inc_isoform_id, sig_transcripts)
        self.assertNotIn(skip_isoform_id, sig_transcripts)
        self.assertNotIn(skip_2_isoform_id, sig_transcripts)
        row = diff_transcripts[inc_isoform_id]
        self.assertEqual(row['asm_id'], asm_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), -0.1852)
        self.assertGreater(float(row['adj_pvalue']), 0.05)


if __name__ == '__main__':
    unittest.main(verbosity=2)
