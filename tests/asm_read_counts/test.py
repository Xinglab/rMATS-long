import os
import os.path
import unittest

import tests.base_test


class AsmReadCountsTest(tests.base_test.BaseTest):
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
        # Only the counts are checked. run_rmats_long() is not called
        self._check_output()

    def _initialize_dirs(self):
        dirs = [self._input_dir, self._out_dir, self._log_dir]
        tests.base_test.remove_then_create_directories(dirs)

    def _create_test_inputs(self):
        self._create_gtf()
        self._create_alignments()

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
        isoform.exons.append(tests.base_test.Region(200, 300))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(800, 900))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)
        # transcript 2 retains the 2nd intron
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(200, 300))
        isoform.exons.append(tests.base_test.Region(500, 900))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)
        # transcript 3 has the 1st 3 exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(200, 300))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(800, 900))
        gene.isoforms.append(isoform)
        # transcript 4 has an exon that overlaps exon 1 and an exon that overlaps exon 4
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(175, 225))
        isoform.exons.append(tests.base_test.Region(1175, 1225))
        gene.isoforms.append(isoform)

        # gene_2
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1 has 4 exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(10200, 10300))
        isoform.exons.append(tests.base_test.Region(10500, 10600))
        isoform.exons.append(tests.base_test.Region(10800, 10900))
        isoform.exons.append(tests.base_test.Region(11100, 11200))
        gene.isoforms.append(isoform)
        # transcript 2 retains the 2nd intron
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(10200, 10300))
        isoform.exons.append(tests.base_test.Region(10500, 10900))
        isoform.exons.append(tests.base_test.Region(11100, 11200))
        gene.isoforms.append(isoform)

        # gene_3
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_3'

        # transcript 1 has 3 exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(20200, 20300))
        isoform.exons.append(tests.base_test.Region(20500, 20600))
        isoform.exons.append(tests.base_test.Region(20800, 20900))
        gene.isoforms.append(isoform)
        # transcript 2 retains the 1st intron
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(20200, 20300))
        isoform.exons.append(tests.base_test.Region(20500, 20900))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        self._chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, self._chrs)

    def _create_alignments(self):
        sample_1_0_sam = os.path.join(self._input_dir, 'g_1_0.sam')
        sample_1_1_sam = os.path.join(self._input_dir, 'g_1_1.sam')
        sample_1_2_sam = os.path.join(self._input_dir, 'g_1_2.sam')
        sample_1_3_sam = os.path.join(self._input_dir, 'g_1_3.sam')
        sample_1_4_sam = os.path.join(self._input_dir, 'g_1_4.sam')
        sample_1_5_sam = os.path.join(self._input_dir, 'g_1_5.sam')
        sample_1_6_sam = os.path.join(self._input_dir, 'g_1_6.sam')
        sample_1_7_sam = os.path.join(self._input_dir, 'g_1_7.sam')
        sample_1_8_sam = os.path.join(self._input_dir, 'g_1_8.sam')
        sample_1_9_sam = os.path.join(self._input_dir, 'g_1_9.sam')
        sample_1_10_sam = os.path.join(self._input_dir, 'g_1_10.sam')
        sample_1_11_sam = os.path.join(self._input_dir, 'g_1_11.sam')
        sample_1_12_sam = os.path.join(self._input_dir, 'g_1_12.sam')
        sample_1_13_sam = os.path.join(self._input_dir, 'g_1_13.sam')
        sample_1_14_sam = os.path.join(self._input_dir, 'g_1_14.sam')
        sample_1_15_sam = os.path.join(self._input_dir, 'g_1_15.sam')
        sample_1_16_sam = os.path.join(self._input_dir, 'g_1_16.sam')
        sample_1_17_sam = os.path.join(self._input_dir, 'g_1_17.sam')
        sample_1_18_sam = os.path.join(self._input_dir, 'g_1_18.sam')
        sample_1_19_sam = os.path.join(self._input_dir, 'g_1_19.sam')
        sample_1_20_sam = os.path.join(self._input_dir, 'g_1_20.sam')
        sample_1_21_sam = os.path.join(self._input_dir, 'g_1_21.sam')
        sample_1_22_sam = os.path.join(self._input_dir, 'g_1_22.sam')
        sample_1_23_sam = os.path.join(self._input_dir, 'g_1_23.sam')
        sample_1_24_sam = os.path.join(self._input_dir, 'g_1_24.sam')
        sample_1_25_sam = os.path.join(self._input_dir, 'g_1_25.sam')
        sample_1_26_sam = os.path.join(self._input_dir, 'g_1_26.sam')
        sample_1_27_sam = os.path.join(self._input_dir, 'g_1_27.sam')
        sample_1_28_sam = os.path.join(self._input_dir, 'g_1_28.sam')
        sample_1_29_sam = os.path.join(self._input_dir, 'g_1_29.sam')
        sample_1_30_sam = os.path.join(self._input_dir, 'g_1_30.sam')
        sample_1_31_sam = os.path.join(self._input_dir, 'g_1_31.sam')
        sample_1_32_sam = os.path.join(self._input_dir, 'g_1_32.sam')
        sample_1_33_sam = os.path.join(self._input_dir, 'g_1_33.sam')
        sample_1_34_sam = os.path.join(self._input_dir, 'g_1_34.sam')
        sample_1_35_sam = os.path.join(self._input_dir, 'g_1_35.sam')
        sample_1_36_sam = os.path.join(self._input_dir, 'g_1_36.sam')
        sample_1_37_sam = os.path.join(self._input_dir, 'g_1_37.sam')
        sample_1_38_sam = os.path.join(self._input_dir, 'g_1_38.sam')
        sample_1_39_sam = os.path.join(self._input_dir, 'g_1_39.sam')
        sample_1_40_sam = os.path.join(self._input_dir, 'g_1_40.sam')
        sample_1_41_sam = os.path.join(self._input_dir, 'g_1_41.sam')
        sams_by_sample = dict()
        sams_by_sample['sample_1_0'] = [sample_1_0_sam]
        sams_by_sample['sample_1_1'] = [sample_1_1_sam]
        sams_by_sample['sample_1_2'] = [sample_1_2_sam]
        sams_by_sample['sample_1_3'] = [sample_1_3_sam]
        sams_by_sample['sample_1_4'] = [sample_1_4_sam]
        sams_by_sample['sample_1_5'] = [sample_1_5_sam]
        sams_by_sample['sample_1_6'] = [sample_1_6_sam]
        sams_by_sample['sample_1_7'] = [sample_1_7_sam]
        sams_by_sample['sample_1_8'] = [sample_1_8_sam]
        sams_by_sample['sample_1_9'] = [sample_1_9_sam]
        sams_by_sample['sample_1_10'] = [sample_1_10_sam]
        sams_by_sample['sample_1_11'] = [sample_1_11_sam]
        sams_by_sample['sample_1_12'] = [sample_1_12_sam]
        sams_by_sample['sample_1_13'] = [sample_1_13_sam]
        sams_by_sample['sample_1_14'] = [sample_1_14_sam]
        sams_by_sample['sample_1_15'] = [sample_1_15_sam]
        sams_by_sample['sample_1_16'] = [sample_1_16_sam]
        sams_by_sample['sample_1_17'] = [sample_1_17_sam]
        sams_by_sample['sample_1_18'] = [sample_1_18_sam]
        sams_by_sample['sample_1_19'] = [sample_1_19_sam]
        sams_by_sample['sample_1_20'] = [sample_1_20_sam]
        sams_by_sample['sample_1_21'] = [sample_1_21_sam]
        sams_by_sample['sample_1_22'] = [sample_1_22_sam]
        sams_by_sample['sample_1_23'] = [sample_1_23_sam]
        sams_by_sample['sample_1_24'] = [sample_1_24_sam]
        sams_by_sample['sample_1_25'] = [sample_1_25_sam]
        sams_by_sample['sample_1_26'] = [sample_1_26_sam]
        sams_by_sample['sample_1_27'] = [sample_1_27_sam]
        sams_by_sample['sample_1_28'] = [sample_1_28_sam]
        sams_by_sample['sample_1_29'] = [sample_1_29_sam]
        sams_by_sample['sample_1_30'] = [sample_1_30_sam]
        sams_by_sample['sample_1_31'] = [sample_1_31_sam]
        sams_by_sample['sample_1_32'] = [sample_1_32_sam]
        sams_by_sample['sample_1_33'] = [sample_1_33_sam]
        sams_by_sample['sample_1_34'] = [sample_1_34_sam]
        sams_by_sample['sample_1_35'] = [sample_1_35_sam]
        sams_by_sample['sample_1_36'] = [sample_1_36_sam]
        sams_by_sample['sample_1_37'] = [sample_1_37_sam]
        sams_by_sample['sample_1_38'] = [sample_1_38_sam]
        sams_by_sample['sample_1_39'] = [sample_1_39_sam]
        sams_by_sample['sample_1_40'] = [sample_1_40_sam]
        sams_by_sample['sample_1_41'] = [sample_1_41_sam]
        self._sams_by_sample = sams_by_sample

        isoforms = list()
        chr_1 = self._chrs[0]
        one_to_4 = chr_1.genes[0].isoforms[0].copy()
        isoforms.append(one_to_4)  # 0
        retained_2 = chr_1.genes[0].isoforms[1].copy()
        isoforms.append(retained_2)  # 1
        one_to_3 = chr_1.genes[0].isoforms[2].copy()
        isoforms.append(one_to_3)  # 2
        overlapped_1_and_4 = chr_1.genes[0].isoforms[3].copy()
        isoforms.append(overlapped_1_and_4)  # 3
        two_to_4 = one_to_4.copy()
        two_to_4.exons.pop(0)
        isoforms.append(two_to_4)  # 4
        one_to_5 = one_to_4.copy()
        one_to_5.exons.append(tests.base_test.Region(1400, 1500))
        isoforms.append(one_to_5)  # 5
        zero_to_3 = one_to_4.copy()
        zero_to_3.exons.pop()
        zero_to_3.exons = [tests.base_test.Region(50, 150)] + zero_to_3.exons
        isoforms.append(zero_to_3)  # 6
        two_3_early_4 = two_to_4.copy()
        two_3_early_4.exons[2].start -= 25
        isoforms.append(two_3_early_4)  # 7
        one_then_3 = one_to_3.copy()
        one_then_3.exons.pop(1)
        isoforms.append(one_then_3)  # 8
        retained_2_then_5 = retained_2.copy()
        retained_2_then_5.exons.append(tests.base_test.Region(1400, 1500))
        isoforms.append(retained_2_then_5)  # 9
        zero_1_split_2_3 = zero_to_3.copy()
        zero_1_split_2_3.exons.insert(2, tests.base_test.Region(500, 525))
        zero_1_split_2_3.exons[3].start = 575
        isoforms.append(zero_1_split_2_3)  # 10
        two_to_4_long_end = two_to_4.copy()
        two_to_4_long_end.exons[2].end += 25
        isoforms.append(two_to_4_long_end)  # 11
        two_to_4_short_start = two_to_4.copy()
        two_to_4_short_start.exons[0].start += 25
        isoforms.append(two_to_4_short_start)  # 12
        two_to_4_long_start_short_end = two_to_4.copy()
        two_to_4_long_start_short_end.exons[0].start -= 25
        two_to_4_long_start_short_end.exons[2].end -= 25
        isoforms.append(two_to_4_long_start_short_end)  # 13
        short_2_end = two_to_4.copy()
        short_2_end.exons.pop()
        short_2_end.exons.pop()
        short_2_end.exons[0].end -= 25
        isoforms.append(short_2_end)  # 14
        short_retained_2nd = retained_2.copy()
        short_retained_2nd.exons.pop()
        short_retained_2nd.exons.pop(0)
        short_retained_2nd.exons[0].start += 25
        short_retained_2nd.exons[0].end -= 25
        isoforms.append(short_retained_2nd)  # 15
        retained_2_short_4 = retained_2.copy()
        retained_2_short_4.exons.pop(0)
        retained_2_short_4.exons[1].end -= 25
        isoforms.append(retained_2_short_4)  # 16
        two_3 = one_to_3.copy()
        two_3.exons.pop(0)
        isoforms.append(two_3)  # 17
        three_4 = two_to_4.copy()
        three_4.exons.pop(0)
        isoforms.append(three_4)  # 18
        three = three_4.copy()
        three.exons.pop()
        isoforms.append(three)  # 19
        zero_early_one_2_gene_1 = one_to_4.copy()
        zero_early_one_2_gene_1.exons.pop()
        zero_early_one_2_gene_1.exons.pop()
        zero_early_one_2_gene_1.exons.insert(0, tests.base_test.Region(0, 100))
        zero_early_one_2_gene_1.exons[1].start -= 25
        # appended later
        # isoforms.append(zero_early_one_2_gene_1)  # 32
        three_long_4_5_gene_1 = one_to_4.copy()
        three_long_4_5_gene_1.exons.pop(0)
        three_long_4_5_gene_1.exons.pop(0)
        three_long_4_5_gene_1.exons.append(tests.base_test.Region(1300, 1400))
        three_long_4_5_gene_1.exons[1].end += 25
        # appended later
        # isoforms.append(three_long_4_5_gene_1)  # 33
        # gene_2
        one_to_4 = chr_1.genes[1].isoforms[0].copy()
        isoforms.append(one_to_4)  # 20
        retained_2 = chr_1.genes[1].isoforms[1].copy()
        isoforms.append(retained_2)  # 21
        retained_2_only = retained_2.copy()
        retained_2_only.exons.pop()
        retained_2_only.exons.pop(0)
        isoforms.append(retained_2_only)  # 22
        retained_2_read_over_downstream = retained_2_only.copy()
        retained_2_read_over_downstream.exons[0].start -= 25
        retained_2_read_over_downstream.exons[0].end -= 25
        isoforms.append(retained_2_read_over_downstream)  # 23
        retained_2_read_over_upstream = retained_2_only.copy()
        retained_2_read_over_upstream.exons[0].start += 25
        retained_2_read_over_upstream.exons[0].end += 25
        isoforms.append(retained_2_read_over_upstream)  # 24
        retained_2_contained = retained_2_only.copy()
        retained_2_contained.exons[0].start += 25
        retained_2_contained.exons[0].end -= 25
        isoforms.append(retained_2_contained)  # 25
        retained_2_then_up = retained_2.copy()
        retained_2_then_up.exons.pop(0)
        isoforms.append(retained_2_then_up)  # 26
        down_then_retained_2 = retained_2.copy()
        down_then_retained_2.exons.pop()
        isoforms.append(down_then_retained_2)  # 27
        one_2_end_in_3 = one_to_4.copy()
        one_2_end_in_3.exons.pop()
        one_2_end_in_3.exons[2].end -= 25
        isoforms.append(one_2_end_in_3)  # 28
        start_in_2_3_4 = one_to_4.copy()
        start_in_2_3_4.exons.pop(0)
        start_in_2_3_4.exons[0].start += 25
        isoforms.append(start_in_2_3_4)  # 29
        zero_early_one_2 = one_to_4.copy()
        zero_early_one_2.exons.pop()
        zero_early_one_2.exons.pop()
        zero_early_one_2.exons.insert(0, tests.base_test.Region(10000, 10100))
        zero_early_one_2.exons[1].start -= 25
        isoforms.append(zero_early_one_2)  # 30
        three_long_4_5 = one_to_4.copy()
        three_long_4_5.exons.pop(0)
        three_long_4_5.exons.pop(0)
        three_long_4_5.exons.append(tests.base_test.Region(11300, 11400))
        three_long_4_5.exons[1].end += 25
        isoforms.append(three_long_4_5)  # 31
        isoforms.append(zero_early_one_2_gene_1)  # 32
        isoforms.append(three_long_4_5_gene_1)  # 33
        # gene_3
        one_to_3 = chr_1.genes[2].isoforms[0].copy()
        isoforms.append(one_to_3)  # 34
        one_retain_2_3 = chr_1.genes[2].isoforms[1].copy()
        isoforms.append(one_retain_2_3)  # 35
        one_to_2 = one_to_3.copy()
        one_to_2.exons.pop(0)
        isoforms.append(one_to_2)  # 36
        one_retain_2 = one_retain_2_3.copy()
        one_retain_2.exons.pop(0)
        isoforms.append(one_retain_2)  # 37
        early_one_to_2 = one_to_2.copy()
        early_one_to_2.exons[1].end += 25
        isoforms.append(early_one_to_2)  # 38
        early_one_retain_2 = one_retain_2.copy()
        early_one_retain_2.exons[0].end += 25
        isoforms.append(early_one_retain_2)  # 39
        late_one_to_2 = one_to_2.copy()
        late_one_to_2.exons[0].start -= 25
        isoforms.append(late_one_to_2)  # 40
        late_one_retain_2 = one_retain_2.copy()
        late_one_retain_2.exons[0].start -= 25
        isoforms.append(late_one_retain_2)  # 41

        isoform_alignments = list()
        for isoform in isoforms:
            align = tests.base_test.Alignment()
            align.chr_name = chr_1.name
            align.match_to_isoform(isoform)
            isoform_alignments.append(align)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(isoform_alignments[0], 1, aligns)
        tests.base_test.write_sam(sample_1_0_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[1], 1, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[2], 1, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[3], 1, aligns)
        tests.base_test.write_sam(sample_1_3_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[4], 1, aligns)
        tests.base_test.write_sam(sample_1_4_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[5], 1, aligns)
        tests.base_test.write_sam(sample_1_5_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[6], 1, aligns)
        tests.base_test.write_sam(sample_1_6_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[7], 1, aligns)
        tests.base_test.write_sam(sample_1_7_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[8], 1, aligns)
        tests.base_test.write_sam(sample_1_8_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[9], 1, aligns)
        tests.base_test.write_sam(sample_1_9_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[10], 1, aligns)
        tests.base_test.write_sam(sample_1_10_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[11], 1, aligns)
        tests.base_test.write_sam(sample_1_11_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[12], 1, aligns)
        tests.base_test.write_sam(sample_1_12_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[13], 1, aligns)
        tests.base_test.write_sam(sample_1_13_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[14], 1, aligns)
        tests.base_test.write_sam(sample_1_14_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[15], 1, aligns)
        tests.base_test.write_sam(sample_1_15_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[16], 1, aligns)
        tests.base_test.write_sam(sample_1_16_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[17], 1, aligns)
        tests.base_test.write_sam(sample_1_17_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[18], 1, aligns)
        tests.base_test.write_sam(sample_1_18_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[19], 1, aligns)
        tests.base_test.write_sam(sample_1_19_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[20], 1, aligns)
        tests.base_test.write_sam(sample_1_20_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[21], 1, aligns)
        tests.base_test.write_sam(sample_1_21_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[22], 1, aligns)
        tests.base_test.write_sam(sample_1_22_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[23], 1, aligns)
        tests.base_test.write_sam(sample_1_23_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[24], 1, aligns)
        tests.base_test.write_sam(sample_1_24_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[25], 1, aligns)
        tests.base_test.write_sam(sample_1_25_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[26], 1, aligns)
        tests.base_test.write_sam(sample_1_26_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[27], 1, aligns)
        tests.base_test.write_sam(sample_1_27_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[28], 1, aligns)
        tests.base_test.write_sam(sample_1_28_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[29], 1, aligns)
        tests.base_test.write_sam(sample_1_29_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[30], 1, aligns)
        tests.base_test.write_sam(sample_1_30_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[31], 1, aligns)
        tests.base_test.write_sam(sample_1_31_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[32], 1, aligns)
        tests.base_test.write_sam(sample_1_32_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[33], 1, aligns)
        tests.base_test.write_sam(sample_1_33_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[34], 1, aligns)
        tests.base_test.write_sam(sample_1_34_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[35], 1, aligns)
        tests.base_test.write_sam(sample_1_35_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[36], 1, aligns)
        tests.base_test.write_sam(sample_1_36_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[37], 1, aligns)
        tests.base_test.write_sam(sample_1_37_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[38], 1, aligns)
        tests.base_test.write_sam(sample_1_38_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[39], 1, aligns)
        tests.base_test.write_sam(sample_1_39_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[40], 1, aligns)
        tests.base_test.write_sam(sample_1_40_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(isoform_alignments[41], 1, aligns)
        tests.base_test.write_sam(sample_1_41_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    def _check_output(self):
        asms_by_isoforms = self.parse_asm_definitions(self._event_dir)
        expected_gene_asm_isoforms = (
            ((176, 226), (1176, 1226)),
            ((201, 301), (501, 601), (801, 901)),
            ((201, 301), (501, 601), (801, 901), (1101, 1201)),
            ((201, 301), (501, 901)),
            ((201, 301), (501, 901), (1101, 1201)),
        )
        expected_ri_asm_isoforms = (
            ((501, 601), (801, 901)),
            ((501, 901), ),
        )
        expected_alt_end_asm_isoforms = (
            ((801, 901), ),
            ((801, 901), (1101, 1201)),
        )
        expected_gene_2_ri_asm_isoforms = (
            ((10501, 10601), (10801, 10901)),
            ((10501, 10901), ),
        )
        expected_gene_3_ri_asm_isoforms = (
            ((20501, 20601), (20801, 20901)),
            ((20501, 20901), ),
        )
        gene_asm_id = None
        ri_asm_id = None
        alt_end_asm_id = None
        gene_2_ri_asm_id = None
        gene_3_ri_asm_id = None
        for isoforms, asm_row in asms_by_isoforms.items():
            asm_id = asm_row['asm_id']
            if isoforms == expected_gene_asm_isoforms:
                gene_asm_id = asm_id
                self.assertEqual(asm_row['start_always_ss'], 'False')
                self.assertEqual(asm_row['start_never_ss'], 'False')
                self.assertEqual(asm_row['end_always_ss'], 'False')
                self.assertEqual(asm_row['end_never_ss'], 'False')
            elif isoforms == expected_ri_asm_isoforms:
                ri_asm_id = asm_id
                self.assertEqual(asm_row['start_always_ss'], 'True')
                self.assertEqual(asm_row['start_never_ss'], 'False')
                self.assertEqual(asm_row['end_always_ss'], 'False')
                self.assertEqual(asm_row['end_never_ss'], 'False')
            elif isoforms == expected_alt_end_asm_isoforms:
                alt_end_asm_id = asm_id
                self.assertEqual(asm_row['start_always_ss'], 'True')
                self.assertEqual(asm_row['start_never_ss'], 'False')
                self.assertEqual(asm_row['end_always_ss'], 'False')
                self.assertEqual(asm_row['end_never_ss'], 'False')
            elif isoforms == expected_gene_2_ri_asm_isoforms:
                gene_2_ri_asm_id = asm_id
                self.assertEqual(asm_row['start_always_ss'], 'True')
                self.assertEqual(asm_row['start_never_ss'], 'False')
                self.assertEqual(asm_row['end_always_ss'], 'True')
                self.assertEqual(asm_row['end_never_ss'], 'False')
            elif isoforms == expected_gene_3_ri_asm_isoforms:
                gene_3_ri_asm_id = asm_id
                self.assertEqual(asm_row['start_always_ss'], 'False')
                self.assertEqual(asm_row['start_never_ss'], 'True')
                self.assertEqual(asm_row['end_always_ss'], 'True')
                self.assertEqual(asm_row['end_never_ss'], 'False')
            else:
                self.fail('unexpected ASM isoforms: {}'.format(isoforms))

        self.assertIsNotNone(gene_asm_id)
        self.assertIsNotNone(ri_asm_id)
        self.assertIsNotNone(alt_end_asm_id)
        self.assertIsNotNone(gene_2_ri_asm_id)
        self.assertIsNotNone(gene_3_ri_asm_id)
        overlap_1_and_4_id = '{}_0'.format(gene_asm_id)
        one_to_3_id = '{}_1'.format(gene_asm_id)
        one_to_4_id = '{}_2'.format(gene_asm_id)
        ri_early_end_id = '{}_3'.format(gene_asm_id)
        retained_2_id = '{}_4'.format(gene_asm_id)

        ri_skip_id = '{}_0'.format(ri_asm_id)
        ri_inc_id = '{}_1'.format(ri_asm_id)

        alt_end_skip_id = '{}_0'.format(alt_end_asm_id)
        alt_end_inc_id = '{}_1'.format(alt_end_asm_id)

        gene_2_ri_skip_id = '{}_0'.format(gene_2_ri_asm_id)
        gene_2_ri_inc_id = '{}_1'.format(gene_2_ri_asm_id)

        gene_3_ri_skip_id = '{}_0'.format(gene_3_ri_asm_id)
        gene_3_ri_inc_id = '{}_1'.format(gene_3_ri_asm_id)

        count_path = os.path.join(self._out_dir, 'read_counts.tsv')
        counts = self.parse_read_compatibility(count_path)
        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_0', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_0', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_0', 0), 1)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_0', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_0', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_0', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_0', 0), 0)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_0', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_0', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_1', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_1', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_1', 0), 0)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_1', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_1', 0), 1)
        self.assertEqual(counts[ri_skip_id].get('sample_1_1', 0), 0)
        self.assertEqual(counts[ri_inc_id].get('sample_1_1', 0), 1)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_1', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_1', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_2', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_2', 0), 1)
        self.assertEqual(counts[one_to_4_id].get('sample_1_2', 0), 1)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_2', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_2', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_2', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_2', 0), 0)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_2', 0), 1)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_2', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_3', 0), 1)
        self.assertEqual(counts[one_to_3_id].get('sample_1_3', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_3', 0), 0)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_3', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_3', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_3', 0), 0)
        self.assertEqual(counts[ri_inc_id].get('sample_1_3', 0), 0)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_3', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_3', 0), 0)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_4', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_4', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_4', 0), 1)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_4', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_4', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_4', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_4', 0), 0)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_4', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_4', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_5', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_5', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_5', 0), 0)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_5', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_5', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_5', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_5', 0), 0)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_5', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_5', 0), 0)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_6', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_6', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_6', 0), 0)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_6', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_6', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_6', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_6', 0), 0)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_6', 0), 1)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_6', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_7', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_7', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_7', 0), 0)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_7', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_7', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_7', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_7', 0), 0)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_7', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_7', 0), 0)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_8', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_8', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_8', 0), 0)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_8', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_8', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_8', 0), 0)
        self.assertEqual(counts[ri_inc_id].get('sample_1_8', 0), 0)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_8', 0), 1)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_8', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_9', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_9', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_9', 0), 0)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_9', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_9', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_9', 0), 0)
        self.assertEqual(counts[ri_inc_id].get('sample_1_9', 0), 1)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_9', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_9', 0), 0)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_10', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_10', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_10', 0), 0)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_10', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_10', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_10', 0), 0)
        self.assertEqual(counts[ri_inc_id].get('sample_1_10', 0), 0)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_10', 0), 1)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_10', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_11', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_11', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_11', 0), 1)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_11', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_11', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_11', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_11', 0), 0)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_11', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_11', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_12', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_12', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_12', 0), 1)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_12', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_12', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_12', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_12', 0), 0)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_12', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_12', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_13', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_13', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_13', 0), 0)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_13', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_13', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_13', 0), 0)
        self.assertEqual(counts[ri_inc_id].get('sample_1_13', 0), 0)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_13', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_13', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_14', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_14', 0), 1)
        self.assertEqual(counts[one_to_4_id].get('sample_1_14', 0), 1)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_14', 0), 1)
        self.assertEqual(counts[retained_2_id].get('sample_1_14', 0), 1)
        self.assertEqual(counts[ri_skip_id].get('sample_1_14', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_14', 0), 1)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_14', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_14', 0), 0)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_15', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_15', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_15', 0), 0)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_15', 0), 1)
        self.assertEqual(counts[retained_2_id].get('sample_1_15', 0), 1)
        self.assertEqual(counts[ri_skip_id].get('sample_1_15', 0), 0)
        self.assertEqual(counts[ri_inc_id].get('sample_1_15', 0), 1)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_15', 0), 1)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_15', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_16', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_16', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_16', 0), 0)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_16', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_16', 0), 1)
        self.assertEqual(counts[ri_skip_id].get('sample_1_16', 0), 0)
        self.assertEqual(counts[ri_inc_id].get('sample_1_16', 0), 1)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_16', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_16', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_17', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_17', 0), 1)
        self.assertEqual(counts[one_to_4_id].get('sample_1_17', 0), 1)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_17', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_17', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_17', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_17', 0), 0)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_17', 0), 1)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_17', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_18', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_18', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_18', 0), 1)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_18', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_18', 0), 1)
        self.assertEqual(counts[ri_skip_id].get('sample_1_18', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_18', 0), 1)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_18', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_18', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_19', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_19', 0), 1)
        self.assertEqual(counts[one_to_4_id].get('sample_1_19', 0), 1)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_19', 0), 1)
        self.assertEqual(counts[retained_2_id].get('sample_1_19', 0), 1)
        self.assertEqual(counts[ri_skip_id].get('sample_1_19', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_19', 0), 1)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_19', 0), 1)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_19', 0), 1)

        self.assertEqual(counts[gene_2_ri_skip_id].get('sample_1_20', 0), 1)
        self.assertEqual(counts[gene_2_ri_inc_id].get('sample_1_20', 0), 0)

        self.assertEqual(counts[gene_2_ri_skip_id].get('sample_1_21', 0), 0)
        self.assertEqual(counts[gene_2_ri_inc_id].get('sample_1_21', 0), 1)

        self.assertEqual(counts[gene_2_ri_skip_id].get('sample_1_22', 0), 0)
        self.assertEqual(counts[gene_2_ri_inc_id].get('sample_1_22', 0), 1)

        self.assertEqual(counts[gene_2_ri_skip_id].get('sample_1_23', 0), 0)
        self.assertEqual(counts[gene_2_ri_inc_id].get('sample_1_23', 0), 0)

        self.assertEqual(counts[gene_2_ri_skip_id].get('sample_1_24', 0), 0)
        self.assertEqual(counts[gene_2_ri_inc_id].get('sample_1_24', 0), 0)

        self.assertEqual(counts[gene_2_ri_skip_id].get('sample_1_25', 0), 0)
        self.assertEqual(counts[gene_2_ri_inc_id].get('sample_1_25', 0), 1)

        self.assertEqual(counts[gene_2_ri_skip_id].get('sample_1_26', 0), 0)
        self.assertEqual(counts[gene_2_ri_inc_id].get('sample_1_26', 0), 1)

        self.assertEqual(counts[gene_2_ri_skip_id].get('sample_1_27', 0), 0)
        self.assertEqual(counts[gene_2_ri_inc_id].get('sample_1_27', 0), 1)

        self.assertEqual(counts[gene_2_ri_skip_id].get('sample_1_28', 0), 1)
        self.assertEqual(counts[gene_2_ri_inc_id].get('sample_1_28', 0), 0)

        self.assertEqual(counts[gene_2_ri_skip_id].get('sample_1_29', 0), 1)
        self.assertEqual(counts[gene_2_ri_inc_id].get('sample_1_29', 0), 0)

        self.assertEqual(counts[gene_2_ri_skip_id].get('sample_1_30', 0), 1)
        self.assertEqual(counts[gene_2_ri_inc_id].get('sample_1_30', 0), 1)

        self.assertEqual(counts[gene_2_ri_skip_id].get('sample_1_31', 0), 1)
        self.assertEqual(counts[gene_2_ri_inc_id].get('sample_1_31', 0), 1)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_32', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_32', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_32', 0), 0)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_32', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_32', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_32', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_32', 0), 1)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_32', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_32', 0), 0)

        self.assertEqual(counts[overlap_1_and_4_id].get('sample_1_33', 0), 0)
        self.assertEqual(counts[one_to_3_id].get('sample_1_33', 0), 0)
        self.assertEqual(counts[one_to_4_id].get('sample_1_33', 0), 0)
        self.assertEqual(counts[ri_early_end_id].get('sample_1_33', 0), 0)
        self.assertEqual(counts[retained_2_id].get('sample_1_33', 0), 0)
        self.assertEqual(counts[ri_skip_id].get('sample_1_33', 0), 1)
        self.assertEqual(counts[ri_inc_id].get('sample_1_33', 0), 1)
        self.assertEqual(counts[alt_end_skip_id].get('sample_1_33', 0), 0)
        self.assertEqual(counts[alt_end_inc_id].get('sample_1_33', 0), 0)

        self.assertEqual(counts[gene_3_ri_skip_id].get('sample_1_34', 0), 1)
        self.assertEqual(counts[gene_3_ri_inc_id].get('sample_1_34', 0), 0)

        self.assertEqual(counts[gene_3_ri_skip_id].get('sample_1_35', 0), 0)
        self.assertEqual(counts[gene_3_ri_inc_id].get('sample_1_35', 0), 1)

        self.assertEqual(counts[gene_3_ri_skip_id].get('sample_1_36', 0), 1)
        self.assertEqual(counts[gene_3_ri_inc_id].get('sample_1_36', 0), 0)

        self.assertEqual(counts[gene_3_ri_skip_id].get('sample_1_37', 0), 0)
        self.assertEqual(counts[gene_3_ri_inc_id].get('sample_1_37', 0), 1)

        self.assertEqual(counts[gene_3_ri_skip_id].get('sample_1_38', 0), 1)
        self.assertEqual(counts[gene_3_ri_inc_id].get('sample_1_38', 0), 0)

        self.assertEqual(counts[gene_3_ri_skip_id].get('sample_1_39', 0), 0)
        self.assertEqual(counts[gene_3_ri_inc_id].get('sample_1_39', 0), 1)

        self.assertEqual(counts[gene_3_ri_skip_id].get('sample_1_40', 0), 0)
        self.assertEqual(counts[gene_3_ri_inc_id].get('sample_1_40', 0), 0)

        self.assertEqual(counts[gene_3_ri_skip_id].get('sample_1_41', 0), 0)
        self.assertEqual(counts[gene_3_ri_inc_id].get('sample_1_41', 0), 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
