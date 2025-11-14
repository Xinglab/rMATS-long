import os
import os.path
import unittest

import tests.base_test


class IsoformDifferencesTest(tests.base_test.BaseTest):
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
                                                 output_full_gene_asm=True)
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
                                                   average_reads_per_group=0)
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

        # transcript 1 exon 1 starts at 200, exon 3 ends at 900
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(200, 300))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(800, 900))
        gene.isoforms.append(isoform)
        # transcript 2 exon 1 ends at 200, exon 3 ends at 1000
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(800, 1000))
        gene.isoforms.append(isoform)
        # transcript 3 all different exons
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(225, 325))
        isoform.exons.append(tests.base_test.Region(525, 625))
        isoform.exons.append(tests.base_test.Region(825, 925))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(10200, 10300))
        isoform.exons.append(tests.base_test.Region(10500, 10600))
        isoform.exons.append(tests.base_test.Region(10900, 11000))
        gene.isoforms.append(isoform)
        # transcript 2 exon 2 overlaps transcript 1 exon 1 at 1bp
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(10200, 10300))
        isoform.exons.append(tests.base_test.Region(10600, 10700))
        isoform.exons.append(tests.base_test.Region(10900, 11000))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_3'

        # transcript 1 single bp exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(20200, 20300))
        isoform.exons.append(tests.base_test.Region(20500, 20500))
        isoform.exons.append(tests.base_test.Region(20800, 20900))
        gene.isoforms.append(isoform)
        # transcript 2 different last exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(20200, 20300))
        isoform.exons.append(tests.base_test.Region(20500, 20500))
        isoform.exons.append(tests.base_test.Region(21100, 21200))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_4'

        # transcript 1 has early transcript start and late transcript end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(30100, 30300))
        isoform.exons.append(tests.base_test.Region(30500, 30600))
        isoform.exons.append(tests.base_test.Region(30800, 30900))
        isoform.exons.append(tests.base_test.Region(31100, 31300))
        gene.isoforms.append(isoform)
        # transcript 2 has late start and early end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(30200, 30300))
        isoform.exons.append(tests.base_test.Region(30500, 30600))
        isoform.exons.append(tests.base_test.Region(30800, 30900))
        isoform.exons.append(tests.base_test.Region(31100, 31200))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_5'

        # transcript 1 is a single exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(40200, 40300))
        gene.isoforms.append(isoform)
        # transcript 2 has late transcript end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(40200, 40400))
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
        isoform_1 = chr_1.genes[0].isoforms[0]
        align_1 = tests.base_test.Alignment()
        align_1.chr_name = chr_1.name
        align_1.match_to_isoform(isoform_1)
        isoform_2 = chr_1.genes[0].isoforms[1]
        align_2 = tests.base_test.Alignment()
        align_2.chr_name = chr_1.name
        align_2.match_to_isoform(isoform_2)
        isoform_3 = chr_1.genes[0].isoforms[2]
        align_3 = tests.base_test.Alignment()
        align_3.chr_name = chr_1.name
        align_3.match_to_isoform(isoform_3)

        isoform_1_no_first = isoform_1.copy()
        isoform_1_no_first.exons.pop(0)
        align_1_no_first = tests.base_test.Alignment()
        align_1_no_first.chr_name = chr_1.name
        align_1_no_first.match_to_isoform(isoform_1_no_first)

        isoform_1_no_last = isoform_1.copy()
        isoform_1_no_last.exons.pop()
        align_1_no_last = tests.base_test.Alignment()
        align_1_no_last.chr_name = chr_1.name
        align_1_no_last.match_to_isoform(isoform_1_no_last)

        isoform_2_no_first = isoform_2.copy()
        isoform_2_no_first.exons.pop(0)
        align_2_no_first = tests.base_test.Alignment()
        align_2_no_first.chr_name = chr_1.name
        align_2_no_first.match_to_isoform(isoform_2_no_first)

        isoform_2_no_last = isoform_2.copy()
        isoform_2_no_last.exons.pop()
        align_2_no_last = tests.base_test.Alignment()
        align_2_no_last.chr_name = chr_1.name
        align_2_no_last.match_to_isoform(isoform_2_no_last)

        # gene_2
        gene_2_isoform_1 = chr_1.genes[1].isoforms[0]
        gene_2_align_1 = tests.base_test.Alignment()
        gene_2_align_1.chr_name = chr_1.name
        gene_2_align_1.match_to_isoform(gene_2_isoform_1)

        gene_2_isoform_2 = chr_1.genes[1].isoforms[1]
        gene_2_align_2 = tests.base_test.Alignment()
        gene_2_align_2.chr_name = chr_1.name
        gene_2_align_2.match_to_isoform(gene_2_isoform_2)

        # gene_3
        gene_3_isoform_1 = chr_1.genes[2].isoforms[0]
        gene_3_align_1 = tests.base_test.Alignment()
        gene_3_align_1.chr_name = chr_1.name
        gene_3_align_1.match_to_isoform(gene_3_isoform_1)

        gene_3_isoform_2 = chr_1.genes[2].isoforms[1]
        gene_3_align_2 = tests.base_test.Alignment()
        gene_3_align_2.chr_name = chr_1.name
        gene_3_align_2.match_to_isoform(gene_3_isoform_2)

        # gene_4
        gene_4_isoform_1 = chr_1.genes[3].isoforms[0]
        gene_4_align_1 = tests.base_test.Alignment()
        gene_4_align_1.chr_name = chr_1.name
        gene_4_align_1.match_to_isoform(gene_4_isoform_1)

        gene_4_isoform_2 = chr_1.genes[3].isoforms[1]
        gene_4_align_2 = tests.base_test.Alignment()
        gene_4_align_2.chr_name = chr_1.name
        gene_4_align_2.match_to_isoform(gene_4_isoform_2)

        gene_4_isoform_1_early_start = gene_4_isoform_1.copy()
        gene_4_isoform_1_early_start.exons.pop()
        gene_4_isoform_1_early_start.exons.pop()
        gene_4_isoform_1_early_start.exons.pop()
        gene_4_isoform_1_early_start.exons[0].end -= 150
        gene_4_align_1_early_start = tests.base_test.Alignment()
        gene_4_align_1_early_start.chr_name = chr_1.name
        gene_4_align_1_early_start.match_to_isoform(
            gene_4_isoform_1_early_start)

        gene_4_isoform_1_late_end = gene_4_isoform_1.copy()
        gene_4_isoform_1_late_end.exons.pop(0)
        gene_4_isoform_1_late_end.exons.pop(0)
        gene_4_isoform_1_late_end.exons.pop(0)
        gene_4_isoform_1_late_end.exons[0].start += 150
        gene_4_align_1_late_end = tests.base_test.Alignment()
        gene_4_align_1_late_end.chr_name = chr_1.name
        gene_4_align_1_late_end.match_to_isoform(gene_4_isoform_1_late_end)

        # gene_5
        gene_5_isoform_1 = chr_1.genes[4].isoforms[0]
        gene_5_align_1 = tests.base_test.Alignment()
        gene_5_align_1.chr_name = chr_1.name
        gene_5_align_1.match_to_isoform(gene_5_isoform_1)

        gene_5_isoform_2 = chr_1.genes[4].isoforms[1]
        gene_5_align_2 = tests.base_test.Alignment()
        gene_5_align_2.chr_name = chr_1.name
        gene_5_align_2.match_to_isoform(gene_5_isoform_2)

        gene_5_isoform_2_late_end = gene_5_isoform_2.copy()
        gene_5_isoform_2_late_end.exons[0].start += 150
        gene_5_align_2_late_end = tests.base_test.Alignment()
        gene_5_align_2_late_end.chr_name = chr_1.name
        gene_5_align_2_late_end.match_to_isoform(gene_5_isoform_2_late_end)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(align_1, 11, aligns)
        tests.base_test.append_copies(align_2, 12, aligns)
        tests.base_test.append_copies(align_3, 13, aligns)
        tests.base_test.append_copies(align_1_no_first, 14, aligns)
        tests.base_test.append_copies(align_1_no_last, 15, aligns)
        tests.base_test.append_copies(align_2_no_first, 16, aligns)
        tests.base_test.append_copies(align_2_no_last, 17, aligns)
        tests.base_test.append_copies(gene_2_align_1, 10, aligns)
        tests.base_test.append_copies(gene_2_align_2, 11, aligns)
        tests.base_test.append_copies(gene_3_align_1, 12, aligns)
        tests.base_test.append_copies(gene_3_align_2, 13, aligns)
        tests.base_test.append_copies(gene_4_align_1, 1, aligns)
        tests.base_test.append_copies(gene_4_align_2, 2, aligns)
        tests.base_test.append_copies(gene_4_align_1_early_start, 3, aligns)
        tests.base_test.append_copies(gene_4_align_1_late_end, 4, aligns)
        tests.base_test.append_copies(gene_5_align_1, 1, aligns)
        tests.base_test.append_copies(gene_5_align_2, 2, aligns)
        tests.base_test.append_copies(gene_5_align_2_late_end, 3, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(align_1, 12, aligns)
        tests.base_test.append_copies(align_2, 13, aligns)
        tests.base_test.append_copies(align_3, 14, aligns)
        tests.base_test.append_copies(align_1_no_first, 15, aligns)
        tests.base_test.append_copies(align_1_no_last, 16, aligns)
        tests.base_test.append_copies(align_2_no_first, 17, aligns)
        tests.base_test.append_copies(align_2_no_last, 18, aligns)
        tests.base_test.append_copies(gene_2_align_1, 11, aligns)
        tests.base_test.append_copies(gene_2_align_2, 12, aligns)
        tests.base_test.append_copies(gene_3_align_1, 13, aligns)
        tests.base_test.append_copies(gene_3_align_2, 14, aligns)
        tests.base_test.append_copies(gene_4_align_1, 2, aligns)
        tests.base_test.append_copies(gene_4_align_2, 3, aligns)
        tests.base_test.append_copies(gene_4_align_1_early_start, 4, aligns)
        tests.base_test.append_copies(gene_4_align_1_late_end, 5, aligns)
        tests.base_test.append_copies(gene_5_align_1, 2, aligns)
        tests.base_test.append_copies(gene_5_align_2, 3, aligns)
        tests.base_test.append_copies(gene_5_align_2_late_end, 4, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(align_1, 21, aligns)
        tests.base_test.append_copies(align_2, 10, aligns)
        tests.base_test.append_copies(align_3, 13, aligns)
        tests.base_test.append_copies(align_1_no_first, 24, aligns)
        tests.base_test.append_copies(align_1_no_last, 25, aligns)
        tests.base_test.append_copies(align_2_no_first, 14, aligns)
        tests.base_test.append_copies(align_2_no_last, 15, aligns)
        tests.base_test.append_copies(gene_2_align_1, 30, aligns)
        tests.base_test.append_copies(gene_2_align_2, 9, aligns)
        tests.base_test.append_copies(gene_3_align_1, 9, aligns)
        tests.base_test.append_copies(gene_3_align_2, 30, aligns)
        tests.base_test.append_copies(gene_4_align_1, 4, aligns)
        tests.base_test.append_copies(gene_4_align_2, 6, aligns)
        tests.base_test.append_copies(gene_4_align_1_early_start, 0, aligns)
        tests.base_test.append_copies(gene_4_align_1_late_end, 0, aligns)
        tests.base_test.append_copies(gene_5_align_1, 3, aligns)
        tests.base_test.append_copies(gene_5_align_2, 3, aligns)
        tests.base_test.append_copies(gene_5_align_2_late_end, 0, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(align_1, 22, aligns)
        tests.base_test.append_copies(align_2, 11, aligns)
        tests.base_test.append_copies(align_3, 14, aligns)
        tests.base_test.append_copies(align_1_no_first, 25, aligns)
        tests.base_test.append_copies(align_1_no_last, 26, aligns)
        tests.base_test.append_copies(align_2_no_first, 15, aligns)
        tests.base_test.append_copies(align_2_no_last, 16, aligns)
        tests.base_test.append_copies(gene_2_align_1, 31, aligns)
        tests.base_test.append_copies(gene_2_align_2, 10, aligns)
        tests.base_test.append_copies(gene_3_align_1, 10, aligns)
        tests.base_test.append_copies(gene_3_align_2, 31, aligns)
        tests.base_test.append_copies(gene_4_align_1, 7, aligns)
        tests.base_test.append_copies(gene_4_align_2, 6, aligns)
        tests.base_test.append_copies(gene_4_align_1_early_start, 0, aligns)
        tests.base_test.append_copies(gene_4_align_1_late_end, 0, aligns)
        tests.base_test.append_copies(gene_5_align_1, 4, aligns)
        tests.base_test.append_copies(gene_5_align_2, 5, aligns)
        tests.base_test.append_copies(gene_5_align_2_late_end, 0, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    def _check_output(self):
        summary_txt_path = os.path.join(self._rmats_long_out, 'summary.txt')
        filt_path = os.path.join(self._rmats_long_out,
                                 'differential_isoforms_filtered.tsv')
        chr_id = '0'
        event_i = '0'
        gene_id = 'gene_1'
        by_gene_dir = os.path.join(self._rmats_long_out, 'results_by_gene')
        gene_dir = os.path.join(by_gene_dir, gene_id)
        asm_id = '{}_{}'.format(chr_id, event_i)
        isoform_1_id = '{}_t_1'.format(gene_id)
        isoform_2_id = '{}_t_2'.format(gene_id)
        isoform_3_id = '{}_t_3'.format(gene_id)

        gene_2_event_i = '1'
        gene_2_id = 'gene_2'
        gene_2_dir = os.path.join(by_gene_dir, gene_2_id)
        gene_2_asm_id = '{}_{}'.format(chr_id, gene_2_event_i)
        gene_2_isoform_1_id = '{}_t_1'.format(gene_2_id)
        gene_2_isoform_2_id = '{}_t_2'.format(gene_2_id)

        gene_3_event_i = '2'
        gene_3_id = 'gene_3'
        gene_3_dir = os.path.join(by_gene_dir, gene_3_id)
        gene_3_asm_id = '{}_{}'.format(chr_id, gene_3_event_i)
        gene_3_isoform_1_id = '{}_t_1'.format(gene_3_id)
        gene_3_isoform_2_id = '{}_t_2'.format(gene_3_id)

        gene_4_event_i = '3'
        gene_4_id = 'gene_4'
        gene_4_dir = os.path.join(by_gene_dir, gene_4_id)
        gene_4_asm_id = '{}_{}'.format(chr_id, gene_4_event_i)
        gene_4_isoform_1_id = '{}_t_1'.format(gene_4_id)
        gene_4_isoform_2_id = '{}_t_2'.format(gene_4_id)

        gene_5_event_i = '4'
        gene_5_id = 'gene_5'
        gene_5_dir = os.path.join(by_gene_dir, gene_5_id)
        gene_5_asm_id = '{}_{}'.format(chr_id, gene_5_event_i)
        gene_5_isoform_1_id = '{}_t_1'.format(gene_5_id)
        gene_5_isoform_2_id = '{}_t_2'.format(gene_5_id)

        count_path = os.path.join(self._rmats_long_out, 'count.tsv')
        counts = self.parse_asm_counts(count_path)
        self.assertEqual(sorted(counts.keys()), [
            isoform_1_id, isoform_2_id, isoform_3_id, gene_2_isoform_1_id,
            gene_2_isoform_2_id, gene_3_isoform_1_id, gene_3_isoform_2_id,
            gene_4_isoform_1_id, gene_4_isoform_2_id, gene_5_isoform_1_id,
            gene_5_isoform_2_id
        ])
        self.assertEqual(counts[isoform_1_id]['sample_1_1'], 40.19)
        self.assertEqual(counts[isoform_2_id]['sample_1_1'], 44.81)
        self.assertEqual(counts[isoform_3_id]['sample_1_1'], 13)
        self.assertEqual(counts[gene_2_isoform_1_id]['sample_1_1'], 10)
        self.assertEqual(counts[gene_2_isoform_2_id]['sample_1_1'], 11)
        self.assertEqual(counts[gene_3_isoform_1_id]['sample_1_1'], 12)
        self.assertEqual(counts[gene_3_isoform_2_id]['sample_1_1'], 13)
        self.assertAlmostEqual(counts[gene_4_isoform_1_id]['sample_1_1'],
                               10,
                               delta=0.02)
        self.assertEqual(counts[gene_4_isoform_2_id].get('sample_1_1', 0), 0)
        self.assertEqual(counts[gene_5_isoform_1_id].get('sample_1_1', 0), 0)
        self.assertAlmostEqual(counts[gene_5_isoform_2_id]['sample_1_1'],
                               6,
                               delta=0.02)
        self.assertEqual(counts[isoform_1_id]['sample_1_2'], 43.19)
        self.assertEqual(counts[isoform_2_id]['sample_1_2'], 47.81)
        self.assertEqual(counts[isoform_3_id]['sample_1_2'], 14)
        self.assertEqual(counts[gene_2_isoform_1_id]['sample_1_2'], 11)
        self.assertEqual(counts[gene_2_isoform_2_id]['sample_1_2'], 12)
        self.assertEqual(counts[gene_3_isoform_1_id]['sample_1_2'], 13)
        self.assertEqual(counts[gene_3_isoform_2_id]['sample_1_2'], 14)
        self.assertAlmostEqual(counts[gene_4_isoform_1_id]['sample_1_2'],
                               14,
                               delta=0.02)
        self.assertEqual(counts[gene_4_isoform_2_id].get('sample_1_2', 0), 0)
        self.assertEqual(counts[gene_5_isoform_1_id].get('sample_1_2', 0), 0)
        self.assertAlmostEqual(counts[gene_5_isoform_2_id]['sample_1_2'],
                               9,
                               delta=0.02)
        self.assertEqual(counts[isoform_1_id]['sample_2_1'], 70.59)
        self.assertEqual(counts[isoform_2_id]['sample_2_1'], 38.41)
        self.assertEqual(counts[isoform_3_id]['sample_2_1'], 13)
        self.assertEqual(counts[gene_2_isoform_1_id]['sample_2_1'], 30)
        self.assertEqual(counts[gene_2_isoform_2_id]['sample_2_1'], 9)
        self.assertEqual(counts[gene_3_isoform_1_id]['sample_2_1'], 9)
        self.assertEqual(counts[gene_3_isoform_2_id]['sample_2_1'], 30)
        self.assertEqual(counts[gene_4_isoform_1_id]['sample_2_1'], 5)
        self.assertEqual(counts[gene_4_isoform_2_id]['sample_2_1'], 5)
        self.assertEqual(counts[gene_5_isoform_1_id]['sample_2_1'], 3)
        self.assertEqual(counts[gene_5_isoform_2_id]['sample_2_1'], 3)
        self.assertEqual(counts[isoform_1_id]['sample_2_2'], 73.57)
        self.assertEqual(counts[isoform_2_id]['sample_2_2'], 41.43)
        self.assertEqual(counts[isoform_3_id]['sample_2_2'], 14)
        self.assertEqual(counts[gene_2_isoform_1_id]['sample_2_2'], 31)
        self.assertEqual(counts[gene_2_isoform_2_id]['sample_2_2'], 10)
        self.assertEqual(counts[gene_3_isoform_1_id]['sample_2_2'], 10)
        self.assertEqual(counts[gene_3_isoform_2_id]['sample_2_2'], 31)
        self.assertEqual(counts[gene_4_isoform_1_id]['sample_2_2'], 6.5)
        self.assertEqual(counts[gene_4_isoform_2_id]['sample_2_2'], 6.5)
        self.assertEqual(counts[gene_5_isoform_1_id]['sample_2_2'], 4.5)
        self.assertEqual(counts[gene_5_isoform_2_id]['sample_2_2'], 4.5)
        self.assert_exists(summary_txt_path)
        summary = self.parse_summary_txt(summary_txt_path)
        self.assertEqual(summary['total ASMs with significant isoforms'], '5')
        self.assertEqual(summary["alternative 5'-splice site"], '1')
        self.assertEqual(summary['complex'], '1')
        self.assertEqual(summary['alternative last exon'], '1')
        self.assertEqual(summary['alternative endpoints'], '2')
        self.assert_exists(filt_path)
        diff_transcripts = self.parse_differential_transcripts(filt_path,
                                                               is_asm=True)
        self.assertIn(isoform_1_id, diff_transcripts)
        row = diff_transcripts[isoform_1_id]
        self.assertEqual(row['asm_id'], asm_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), -0.1637)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertGreater(float(row['lr']), 0)
        self.assertGreater(float(row['df']), 0)

        self.assertIn(gene_2_isoform_1_id, diff_transcripts)
        row = diff_transcripts[gene_2_isoform_1_id]
        self.assertEqual(row['asm_id'], gene_2_asm_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), -0.2854)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertGreater(float(row['lr']), 0)
        self.assertGreater(float(row['df']), 0)

        self.assertIn(gene_3_isoform_1_id, diff_transcripts)
        row = diff_transcripts[gene_3_isoform_1_id]
        self.assertEqual(row['asm_id'], gene_3_asm_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), 0.2434)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertGreater(float(row['lr']), 0)
        self.assertGreater(float(row['df']), 0)

        self.assertIn(gene_4_isoform_1_id, diff_transcripts)
        row = diff_transcripts[gene_4_isoform_1_id]
        self.assertEqual(row['asm_id'], gene_4_asm_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), 0.5)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertGreater(float(row['lr']), 0)
        self.assertGreater(float(row['df']), 0)

        self.assertIn(gene_5_isoform_1_id, diff_transcripts)
        row = diff_transcripts[gene_5_isoform_1_id]
        self.assertEqual(row['asm_id'], gene_5_asm_id)
        self.assertAlmostEqual(float(row['delta_isoform_proportion']), -0.5)
        self.assertLess(float(row['adj_pvalue']), 0.05)
        self.assertGreater(float(row['lr']), 0)
        self.assertGreater(float(row['df']), 0)

        isoform_diff_name = '{}_isoform_differences_{}_to_{}.tsv'.format(
            asm_id, isoform_1_id, isoform_2_id)
        diff_path = os.path.join(gene_dir, isoform_diff_name)
        isoform_diffs = self.parse_isoform_diff_tsv(diff_path)
        self.assertEqual(len(isoform_diffs), 1)
        diff_row = isoform_diffs[0]
        self.assertEqual(diff_row['event'], 'A5SS')
        self.assertEqual(diff_row['coordinates'],
                         'chr_1:101:201:+;chr_1:101:301:+')

        isoform_diff_name = '{}_isoform_differences_{}_to_{}.tsv'.format(
            gene_2_asm_id, gene_2_isoform_1_id, gene_2_isoform_2_id)
        diff_path = os.path.join(gene_2_dir, isoform_diff_name)
        isoform_diffs = self.parse_isoform_diff_tsv(diff_path)
        self.assertEqual(len(isoform_diffs), 1)
        diff_row = isoform_diffs[0]
        self.assertEqual(diff_row['event'], 'COMPLEX')
        self.assertEqual(
            diff_row['coordinates'],
            ('chr_1:10301:10301:+;chr_1:10501:10601:+;chr_1:10901:10901:+;'
             'chr_1:10301:10301:+;chr_1:10601:10701:+;chr_1:10901:10901:+'))

        isoform_diff_name = '{}_isoform_differences_{}_to_{}.tsv'.format(
            gene_3_asm_id, gene_3_isoform_2_id, gene_3_isoform_1_id)
        diff_path = os.path.join(gene_3_dir, isoform_diff_name)
        isoform_diffs = self.parse_isoform_diff_tsv(diff_path)
        self.assertEqual(len(isoform_diffs), 1)
        diff_row = isoform_diffs[0]
        self.assertEqual(diff_row['event'], 'ALE')
        self.assertEqual(diff_row['coordinates'],
                         'chr_1:20801:20901:+;chr_1:21101:21201:+')

        isoform_diff_name = '{}_isoform_differences_{}_to_{}.tsv'.format(
            gene_4_asm_id, gene_4_isoform_2_id, gene_4_isoform_1_id)
        diff_path = os.path.join(gene_4_dir, isoform_diff_name)
        isoform_diffs = self.parse_isoform_diff_tsv(diff_path)
        self.assertEqual(len(isoform_diffs), 1)
        diff_row = isoform_diffs[0]
        self.assertEqual(diff_row['event'], 'ENDS')
        self.assertEqual(diff_row['coordinates'],
                         'chr_1:31201:30201:-;chr_1:31301:30101:-')

        isoform_diff_name = '{}_isoform_differences_{}_to_{}.tsv'.format(
            gene_5_asm_id, gene_5_isoform_1_id, gene_5_isoform_2_id)
        diff_path = os.path.join(gene_5_dir, isoform_diff_name)
        isoform_diffs = self.parse_isoform_diff_tsv(diff_path)
        self.assertEqual(len(isoform_diffs), 1)
        diff_row = isoform_diffs[0]
        self.assertEqual(diff_row['event'], 'ENDS')
        self.assertEqual(diff_row['coordinates'],
                         'chr_1:40201:40301:+;chr_1:40201:40401:+')


if __name__ == '__main__':
    unittest.main(verbosity=2)
