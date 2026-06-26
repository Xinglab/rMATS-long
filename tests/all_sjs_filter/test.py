import os
import os.path
import unittest

import tests.base_test


class AllSjsFilterTest(tests.base_test.BaseTest):
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
        initial_event_dir = self.run_detect_events(self._gtf_dir,
                                                   self._out_dir,
                                                   self._log_dir)
        initial_asm_counts = self.run_count_reads(self._gtf_dir,
                                                  self._align_dir,
                                                  initial_event_dir,
                                                  self._out_dir, self._log_dir)
        updated_event_details = self.run_filter_isoforms(initial_event_dir,
                                                         initial_asm_counts,
                                                         self._out_dir,
                                                         self._log_dir,
                                                         all_sjs_count=2)
        self._event_dir = updated_event_details['event_dir']
        self._asm_counts = updated_event_details['counts_dir']
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

        # transcript 1: 1,2,3,4,5
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1000, 1100))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        gene.isoforms.append(isoform)
        # transcript 2: 1,3,5
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        gene.isoforms.append(isoform)
        # transcript 3: 1,4,5
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
        # 1,3,5
        no_2_4_isoform = chr_1.genes[0].isoforms[1]
        no_2_4_align = tests.base_test.Alignment()
        no_2_4_align.chr_name = chr_1.name
        no_2_4_align.match_to_isoform(no_2_4_isoform)
        # 1,4,5
        no_2_3_isoform = chr_1.genes[0].isoforms[2]
        no_2_3_align = tests.base_test.Alignment()
        no_2_3_align.chr_name = chr_1.name
        no_2_3_align.match_to_isoform(no_2_3_isoform)
        # 1,2,3,5
        no_4_isoform = all_isoform.copy()
        no_4_isoform.exons.pop(3)
        no_4_align = tests.base_test.Alignment()
        no_4_align.chr_name = chr_1.name
        no_4_align.match_to_isoform(no_4_isoform)
        # 1,2
        sj_1_2_isoform = all_isoform.copy()
        sj_1_2_isoform.exons = sj_1_2_isoform.exons[:2]
        sj_1_2_align = tests.base_test.Alignment()
        sj_1_2_align.chr_name = chr_1.name
        sj_1_2_align.match_to_isoform(sj_1_2_isoform)
        # 1,3
        sj_1_3_isoform = no_2_4_isoform.copy()
        sj_1_3_isoform.exons.pop(2)
        sj_1_3_align = tests.base_test.Alignment()
        sj_1_3_align.chr_name = chr_1.name
        sj_1_3_align.match_to_isoform(sj_1_3_isoform)
        # 1,4
        sj_1_4_isoform = no_2_3_isoform.copy()
        sj_1_4_isoform.exons.pop(2)
        sj_1_4_align = tests.base_test.Alignment()
        sj_1_4_align.chr_name = chr_1.name
        sj_1_4_align.match_to_isoform(sj_1_4_isoform)
        # 2,3
        sj_2_3_isoform = all_isoform.copy()
        sj_2_3_isoform.exons = sj_2_3_isoform.exons[1:3]
        sj_2_3_align = tests.base_test.Alignment()
        sj_2_3_align.chr_name = chr_1.name
        sj_2_3_align.match_to_isoform(sj_2_3_isoform)
        # 3,4
        sj_3_4_isoform = all_isoform.copy()
        sj_3_4_isoform.exons = sj_3_4_isoform.exons[2:4]
        sj_3_4_align = tests.base_test.Alignment()
        sj_3_4_align.chr_name = chr_1.name
        sj_3_4_align.match_to_isoform(sj_3_4_isoform)
        # 3,5
        sj_3_5_isoform = no_2_4_isoform.copy()
        sj_3_5_isoform.exons.pop(0)
        sj_3_5_align = tests.base_test.Alignment()
        sj_3_5_align.chr_name = chr_1.name
        sj_3_5_align.match_to_isoform(sj_3_5_isoform)
        # 4,5
        sj_4_5_isoform = no_2_3_isoform.copy()
        sj_4_5_isoform.exons.pop(0)
        sj_4_5_align = tests.base_test.Alignment()
        sj_4_5_align.chr_name = chr_1.name
        sj_4_5_align.match_to_isoform(sj_4_5_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(all_align, 5, aligns)
        tests.base_test.append_copies(no_2_4_align, 10, aligns)
        tests.base_test.append_copies(no_2_3_align, 1, aligns)
        tests.base_test.append_copies(no_4_align, 0, aligns)
        tests.base_test.append_copies(sj_1_2_align, 2, aligns)
        tests.base_test.append_copies(sj_1_3_align, 5, aligns)
        tests.base_test.append_copies(sj_1_4_align, 0, aligns)
        tests.base_test.append_copies(sj_2_3_align, 2, aligns)
        tests.base_test.append_copies(sj_3_4_align, 2, aligns)
        tests.base_test.append_copies(sj_3_5_align, 5, aligns)
        tests.base_test.append_copies(sj_4_5_align, 2, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 6, aligns)
        tests.base_test.append_copies(no_2_4_align, 11, aligns)
        tests.base_test.append_copies(no_2_3_align, 0, aligns)
        tests.base_test.append_copies(no_4_align, 0, aligns)
        tests.base_test.append_copies(sj_1_2_align, 1, aligns)
        tests.base_test.append_copies(sj_1_3_align, 4, aligns)
        tests.base_test.append_copies(sj_1_4_align, 0, aligns)
        tests.base_test.append_copies(sj_2_3_align, 1, aligns)
        tests.base_test.append_copies(sj_3_4_align, 1, aligns)
        tests.base_test.append_copies(sj_3_5_align, 4, aligns)
        tests.base_test.append_copies(sj_4_5_align, 1, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 5, aligns)
        tests.base_test.append_copies(no_2_4_align, 3, aligns)
        tests.base_test.append_copies(no_2_3_align, 0, aligns)
        tests.base_test.append_copies(no_4_align, 2, aligns)
        tests.base_test.append_copies(sj_1_2_align, 5, aligns)
        tests.base_test.append_copies(sj_1_3_align, 2, aligns)
        tests.base_test.append_copies(sj_1_4_align, 0, aligns)
        tests.base_test.append_copies(sj_2_3_align, 5, aligns)
        tests.base_test.append_copies(sj_3_4_align, 5, aligns)
        tests.base_test.append_copies(sj_3_5_align, 2, aligns)
        tests.base_test.append_copies(sj_4_5_align, 5, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(all_align, 6, aligns)
        tests.base_test.append_copies(no_2_4_align, 2, aligns)
        tests.base_test.append_copies(no_2_3_align, 0, aligns)
        tests.base_test.append_copies(no_4_align, 0, aligns)
        tests.base_test.append_copies(sj_1_2_align, 6, aligns)
        tests.base_test.append_copies(sj_1_3_align, 1, aligns)
        tests.base_test.append_copies(sj_1_4_align, 0, aligns)
        tests.base_test.append_copies(sj_2_3_align, 6, aligns)
        tests.base_test.append_copies(sj_3_4_align, 6, aligns)
        tests.base_test.append_copies(sj_3_5_align, 1, aligns)
        tests.base_test.append_copies(sj_4_5_align, 6, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    # require 2 reads for each isoform
    # ASMs
    # skip 2,3 (filtered out):
    # 1,2,3,4
    # 1,3,4 (no full reads)
    # 1,4 (1 full read)
    #
    # entire gene:
    # 1,2,3,4,5
    # 1,2,3,5 (2 full reads)
    # 1,3,4,5 (no full reads)
    # 1,3,5
    # 1,4,5 (1 full read)
    #
    # SE 2:
    # 1,2,3
    # 1,3
    #
    # SE 4:
    # 3,4,5
    # 3,5
    def _check_output(self):
        diff_path = os.path.join(self._rmats_long_out,
                                 'differential_isoforms.tsv')
        chr_id = '0'
        event_i = '1'
        asm_id = '{}_{}'.format(chr_id, event_i)
        all_isoform_id = '{}_0'.format(asm_id)
        no_4_isoform_id = '{}_1'.format(asm_id)
        # no_2_isoform_id = '{}_2'.format(asm_id)
        no_2_4_isoform_id = '{}_3'.format(asm_id)
        # no_2_3_isoform_id = '{}_4'.format(asm_id)
        event_i = '2'
        asm_id = '{}_{}'.format(chr_id, event_i)
        inc_2_isoform_id = '{}_0'.format(asm_id)
        skip_2_isoform_id = '{}_1'.format(asm_id)
        event_i = '3'
        asm_id = '{}_{}'.format(chr_id, event_i)
        inc_4_isoform_id = '{}_0'.format(asm_id)
        skip_4_isoform_id = '{}_1'.format(asm_id)
        diff_transcripts = self.parse_differential_transcripts(diff_path,
                                                               is_asm=True)
        self.assertIn(all_isoform_id, diff_transcripts)
        self.assertIn(no_4_isoform_id, diff_transcripts)
        self.assertIn(no_2_4_isoform_id, diff_transcripts)
        self.assertIn(inc_2_isoform_id, diff_transcripts)
        self.assertIn(skip_2_isoform_id, diff_transcripts)
        self.assertIn(inc_4_isoform_id, diff_transcripts)
        self.assertIn(skip_4_isoform_id, diff_transcripts)
        self.assertEqual(len(diff_transcripts), 7)


if __name__ == '__main__':
    unittest.main(verbosity=2)
