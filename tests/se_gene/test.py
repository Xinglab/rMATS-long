import os
import os.path
import unittest

import tests.base_test


class SkippedExonGeneTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_dir = os.path.dirname(__file__)
        self._input_dir = os.path.join(self._test_dir, 'input')
        self._out_dir = os.path.join(self._test_dir, 'out')
        self._log_dir = os.path.join(self._test_dir, 'log')

    def test(self):
        self._initialize_dirs()
        self._create_test_inputs()
        self._run_rmats_long()
        self._check_output()

    def _initialize_dirs(self):
        dirs = [self._input_dir, self._out_dir, self._log_dir]
        tests.base_test.remove_then_create_directories(dirs)

    def _create_test_inputs(self):
        self._chr = 'chr_1'
        self._gene_id = 'gene_1'
        self._transcript_id_1 = 'gene_1_t_1'
        self._transcript_id_2 = 'gene_1_t_2'
        self._samples = [
            'sample_1_1', 'sample_1_2', 'sample_2_1', 'sample_2_2'
        ]
        self._create_abundance()
        self._create_gtf()
        self._create_group_1()
        self._create_group_2()

    def _create_abundance(self):
        self._abundance = os.path.join(self._input_dir, 'abun.esp')
        with open(self._abundance, 'wt') as handle:
            headers = ['transcript_ID', 'transcript_name', 'gene_ID'
                       ] + self._samples
            tests.base_test.write_tsv_line(handle, headers)
            t1_columns = [
                self._transcript_id_1, 'NA', self._gene_id, '30', '32', '20',
                '18'
            ]
            tests.base_test.write_tsv_line(handle, t1_columns)
            t2_columns = [
                self._transcript_id_2, 'NA', self._gene_id, '10', '9', '20',
                '21'
            ]
            tests.base_test.write_tsv_line(handle, t2_columns)

    def _create_gtf(self):
        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        strand = '+'
        with open(self._gtf, 'wt') as handle:
            tests.base_test.write_gtf_gene_line(handle, self._chr,
                                                self._gene_id, 100, 1000,
                                                strand)
            # transcript 1 includes the exon
            tests.base_test.write_gtf_transcript_line(handle, self._chr,
                                                      self._gene_id,
                                                      self._transcript_id_1,
                                                      100, 1000, strand)
            tests.base_test.write_gtf_exon_line(handle, self._chr,
                                                self._gene_id,
                                                self._transcript_id_1, 100,
                                                200, strand)
            tests.base_test.write_gtf_exon_line(handle, self._chr,
                                                self._gene_id,
                                                self._transcript_id_1, 400,
                                                500, strand)
            tests.base_test.write_gtf_exon_line(handle, self._chr,
                                                self._gene_id,
                                                self._transcript_id_1, 900,
                                                1000, strand)
            # transcript 2 skips the exon
            tests.base_test.write_gtf_transcript_line(handle, self._chr,
                                                      self._gene_id,
                                                      self._transcript_id_2,
                                                      100, 1000, strand)
            tests.base_test.write_gtf_exon_line(handle, self._chr,
                                                self._gene_id,
                                                self._transcript_id_2, 100,
                                                200, strand)
            tests.base_test.write_gtf_exon_line(handle, self._chr,
                                                self._gene_id,
                                                self._transcript_id_2, 900,
                                                1000, strand)

    def _create_group_1(self):
        self._group_1 = os.path.join(self._input_dir, 'group_1.txt')
        with open(self._group_1, 'wt') as handle:
            handle.write('{},{}\n'.format(self._samples[0], self._samples[1]))

    def _create_group_2(self):
        self._group_2 = os.path.join(self._input_dir, 'group_2.txt')
        with open(self._group_2, 'wt') as handle:
            handle.write('{},{}\n'.format(self._samples[2], self._samples[3]))

    def _run_rmats_long(self):
        command = [
            self._python_exe, self._rmats_long_py, '--abundance',
            self._abundance, '--updated-gtf', self._gtf, '--group-1',
            self._group_1, '--group-2', self._group_2, '--out-dir',
            self._out_dir
        ]
        log_path = os.path.join(self._log_dir, 'rmats_long.log')
        tests.base_test.run_command(command, log=log_path)

    def _check_output(self):
        summary_txt_path = os.path.join(self._out_dir, 'summary.txt')
        summary_plot_path = os.path.join(self._out_dir, 'summary_plot.png')
        filt_path = os.path.join(self._out_dir,
                                 'differential_transcripts_filtered.tsv')
        by_gene_dir = os.path.join(self._out_dir, 'results_by_gene')
        gene_dir = os.path.join(by_gene_dir, self._gene_id)
        structure_path = os.path.join(gene_dir,
                                      '{}_structure.png'.format(self._gene_id))
        abundance_path = os.path.join(gene_dir,
                                      '{}_abundance.png'.format(self._gene_id))
        self.assert_exists(summary_txt_path)
        self.assert_exists(summary_plot_path)
        summary = self.parse_summary_txt(summary_txt_path)
        self.assertEqual(summary['total genes with significant isoforms'], '1')
        self.assertEqual(summary['exon skipping'], '1')
        self.assert_exists(filt_path)
        diff_transcripts = self.parse_differential_transcripts(filt_path)
        self.assertIn(self._transcript_id_1, diff_transcripts)
        row = diff_transcripts[self._transcript_id_1]
        self.assertEqual(row['gene_id'], self._gene_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), 0.2845)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assert_exists(structure_path)
        self.assert_exists(abundance_path)


if __name__ == '__main__':
    unittest.main(verbosity=2)
