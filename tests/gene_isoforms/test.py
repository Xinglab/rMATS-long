import os
import os.path
import unittest

import tests.base_test


class GeneIsoformsTest(tests.base_test.BaseTest):
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
        self._event_dir = self.run_detect_events(
            self._gtf_dir,
            self._out_dir,
            self._log_dir,
            align_dir=self._align_dir,
            min_reads=5,
            output_full_gene_asm=True,
            simplify_gene_isoform_endpoints=True,
            filter_gene_isoforms_by_edge=True)
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

        # transcript 1 has exons 1,2,3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(700, 800))
        gene.isoforms.append(isoform)
        # transcript 2 has an early exon 1 start
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(80, 200))
        isoform.exons.append(tests.base_test.Region(400, 500))
        isoform.exons.append(tests.base_test.Region(700, 800))
        gene.isoforms.append(isoform)
        # transcript 3 skips exons 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(700, 800))
        gene.isoforms.append(isoform)
        # transcript 4 skips exon 2 and has an early exon 3 end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(700, 780))
        gene.isoforms.append(isoform)

        # gene_2
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1 has exons 1,2,3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        isoform.exons.append(tests.base_test.Region(1400, 1500))
        isoform.exons.append(tests.base_test.Region(1700, 1800))
        gene.isoforms.append(isoform)
        # transcript 2 has an early exon 1 start and early exon 3 end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(1080, 1200))
        isoform.exons.append(tests.base_test.Region(1400, 1500))
        isoform.exons.append(tests.base_test.Region(1700, 1780))
        gene.isoforms.append(isoform)
        # transcript 3 skips exons 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        isoform.exons.append(tests.base_test.Region(1700, 1800))
        gene.isoforms.append(isoform)

        # gene_3
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_3'

        # transcript 1 has exons 1,2,3,4,5,6
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 2500))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3000, 3100))
        isoform.exons.append(tests.base_test.Region(3300, 3400))
        isoform.exons.append(tests.base_test.Region(3600, 3700))
        gene.isoforms.append(isoform)
        # transcript 2 has exons 1,3,4,5,6
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3000, 3100))
        isoform.exons.append(tests.base_test.Region(3300, 3400))
        isoform.exons.append(tests.base_test.Region(3600, 3700))
        gene.isoforms.append(isoform)
        # transcript 3 has exons 1,2,3,5,6
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 2500))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3300, 3400))
        isoform.exons.append(tests.base_test.Region(3600, 3700))
        gene.isoforms.append(isoform)
        # transcript 4 has exons 1,3,5,6
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3300, 3400))
        isoform.exons.append(tests.base_test.Region(3600, 3700))
        gene.isoforms.append(isoform)
        # transcript 5 has exons 1,2,3,6
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_5'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2400, 2500))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(3600, 3700))
        gene.isoforms.append(isoform)

        # gene_4: like gene_3 but minus strand
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_4'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(4100, 4200))
        isoform.exons.append(tests.base_test.Region(4400, 4500))
        isoform.exons.append(tests.base_test.Region(4700, 4800))
        isoform.exons.append(tests.base_test.Region(5000, 5100))
        isoform.exons.append(tests.base_test.Region(5300, 5400))
        isoform.exons.append(tests.base_test.Region(5600, 5700))
        gene.isoforms.append(isoform)
        # transcript 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(4100, 4200))
        isoform.exons.append(tests.base_test.Region(4700, 4800))
        isoform.exons.append(tests.base_test.Region(5000, 5100))
        isoform.exons.append(tests.base_test.Region(5300, 5400))
        isoform.exons.append(tests.base_test.Region(5600, 5700))
        gene.isoforms.append(isoform)
        # transcript 3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(4100, 4200))
        isoform.exons.append(tests.base_test.Region(4400, 4500))
        isoform.exons.append(tests.base_test.Region(4700, 4800))
        isoform.exons.append(tests.base_test.Region(5300, 5400))
        isoform.exons.append(tests.base_test.Region(5600, 5700))
        gene.isoforms.append(isoform)
        # transcript 4
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(4100, 4200))
        isoform.exons.append(tests.base_test.Region(4700, 4800))
        isoform.exons.append(tests.base_test.Region(5300, 5400))
        isoform.exons.append(tests.base_test.Region(5600, 5700))
        gene.isoforms.append(isoform)
        # transcript 5
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_5'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(4100, 4200))
        isoform.exons.append(tests.base_test.Region(4400, 4500))
        isoform.exons.append(tests.base_test.Region(4700, 4800))
        isoform.exons.append(tests.base_test.Region(5600, 5700))
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
        # 1,2,3
        inc_isoform = chr_1.genes[0].isoforms[0]
        inc_align = tests.base_test.Alignment()
        inc_align.chr_name = chr_1.name
        inc_align.match_to_isoform(inc_isoform)

        # early_1,2,3
        early_inc_isoform = chr_1.genes[0].isoforms[1]
        early_inc_align = tests.base_test.Alignment()
        early_inc_align.chr_name = chr_1.name
        early_inc_align.match_to_isoform(early_inc_isoform)

        # 1,3
        skip_isoform = chr_1.genes[0].isoforms[2]
        skip_align = tests.base_test.Alignment()
        skip_align.chr_name = chr_1.name
        skip_align.match_to_isoform(skip_isoform)

        # 1,early_3
        skip_early_isoform = chr_1.genes[0].isoforms[3]
        skip_early_align = tests.base_test.Alignment()
        skip_early_align.chr_name = chr_1.name
        skip_early_align.match_to_isoform(skip_early_isoform)

        # gene_2
        # 1,2,3
        g_2_inc_isoform = chr_1.genes[1].isoforms[0]
        g_2_inc_align = tests.base_test.Alignment()
        g_2_inc_align.chr_name = chr_1.name
        g_2_inc_align.match_to_isoform(g_2_inc_isoform)

        # early_1,2,early_3
        g_2_early_inc_isoform = chr_1.genes[1].isoforms[1]
        g_2_early_inc_align = tests.base_test.Alignment()
        g_2_early_inc_align.chr_name = chr_1.name
        g_2_early_inc_align.match_to_isoform(g_2_early_inc_isoform)

        # 1,3
        g_2_skip_isoform = chr_1.genes[1].isoforms[2]
        g_2_skip_align = tests.base_test.Alignment()
        g_2_skip_align.chr_name = chr_1.name
        g_2_skip_align.match_to_isoform(g_2_skip_isoform)

        # gene_3
        # 1-6
        g_3_all_isoform = chr_1.genes[2].isoforms[0]
        g_3_all_align = tests.base_test.Alignment()
        g_3_all_align.chr_name = chr_1.name
        g_3_all_align.match_to_isoform(g_3_all_isoform)

        # skip 2
        g_3_skip_2_isoform = chr_1.genes[2].isoforms[1]
        g_3_skip_2_align = tests.base_test.Alignment()
        g_3_skip_2_align.chr_name = chr_1.name
        g_3_skip_2_align.match_to_isoform(g_3_skip_2_isoform)

        # skip 4
        g_3_skip_4_isoform = chr_1.genes[2].isoforms[2]
        g_3_skip_4_align = tests.base_test.Alignment()
        g_3_skip_4_align.chr_name = chr_1.name
        g_3_skip_4_align.match_to_isoform(g_3_skip_4_isoform)

        # skip 2 and skip 4
        g_3_skip_2_and_4_isoform = chr_1.genes[2].isoforms[3]
        g_3_skip_2_and_4_align = tests.base_test.Alignment()
        g_3_skip_2_and_4_align.chr_name = chr_1.name
        g_3_skip_2_and_4_align.match_to_isoform(g_3_skip_2_and_4_isoform)

        # skip 4 and 5
        g_3_skip_4_and_5_isoform = chr_1.genes[2].isoforms[4]
        g_3_skip_4_and_5_align = tests.base_test.Alignment()
        g_3_skip_4_and_5_align.chr_name = chr_1.name
        g_3_skip_4_and_5_align.match_to_isoform(g_3_skip_4_and_5_isoform)

        # gene_4
        # 1-6
        g_4_all_isoform = chr_1.genes[3].isoforms[0]
        g_4_all_align = tests.base_test.Alignment()
        g_4_all_align.chr_name = chr_1.name
        g_4_all_align.match_to_isoform(g_4_all_isoform)

        # skip 2
        g_4_skip_2_isoform = chr_1.genes[3].isoforms[1]
        g_4_skip_2_align = tests.base_test.Alignment()
        g_4_skip_2_align.chr_name = chr_1.name
        g_4_skip_2_align.match_to_isoform(g_4_skip_2_isoform)

        # skip 4
        g_4_skip_4_isoform = chr_1.genes[3].isoforms[2]
        g_4_skip_4_align = tests.base_test.Alignment()
        g_4_skip_4_align.chr_name = chr_1.name
        g_4_skip_4_align.match_to_isoform(g_4_skip_4_isoform)

        # skip 2 and skip 4
        g_4_skip_2_and_4_isoform = chr_1.genes[3].isoforms[3]
        g_4_skip_2_and_4_align = tests.base_test.Alignment()
        g_4_skip_2_and_4_align.chr_name = chr_1.name
        g_4_skip_2_and_4_align.match_to_isoform(g_4_skip_2_and_4_isoform)

        # skip 4 and 5
        g_4_skip_4_and_5_isoform = chr_1.genes[3].isoforms[4]
        g_4_skip_4_and_5_align = tests.base_test.Alignment()
        g_4_skip_4_and_5_align.chr_name = chr_1.name
        g_4_skip_4_and_5_align.match_to_isoform(g_4_skip_4_and_5_isoform)

        next_i = 0
        aligns = list()
        tests.base_test.append_copies(inc_align, 1, aligns)
        tests.base_test.append_copies(early_inc_align, 1, aligns)
        tests.base_test.append_copies(skip_align, 1, aligns)
        tests.base_test.append_copies(skip_early_align, 1, aligns)
        tests.base_test.append_copies(g_2_inc_align, 1, aligns)
        tests.base_test.append_copies(g_2_early_inc_align, 1, aligns)
        tests.base_test.append_copies(g_2_skip_align, 2, aligns)
        tests.base_test.append_copies(g_3_all_align, 1, aligns)
        tests.base_test.append_copies(g_3_skip_2_align, 1, aligns)
        tests.base_test.append_copies(g_3_skip_4_align, 1, aligns)
        tests.base_test.append_copies(g_3_skip_2_and_4_align, 1, aligns)
        tests.base_test.append_copies(g_3_skip_4_and_5_align, 1, aligns)
        tests.base_test.append_copies(g_4_all_align, 1, aligns)
        tests.base_test.append_copies(g_4_skip_2_align, 1, aligns)
        tests.base_test.append_copies(g_4_skip_4_align, 1, aligns)
        tests.base_test.append_copies(g_4_skip_2_and_4_align, 1, aligns)
        tests.base_test.append_copies(g_4_skip_4_and_5_align, 1, aligns)
        tests.base_test.write_sam(sample_1_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 1, aligns)
        tests.base_test.append_copies(early_inc_align, 1, aligns)
        tests.base_test.append_copies(skip_align, 1, aligns)
        tests.base_test.append_copies(skip_early_align, 1, aligns)
        tests.base_test.append_copies(g_2_inc_align, 1, aligns)
        tests.base_test.append_copies(g_2_early_inc_align, 1, aligns)
        tests.base_test.append_copies(g_2_skip_align, 2, aligns)
        tests.base_test.append_copies(g_3_all_align, 2, aligns)
        tests.base_test.append_copies(g_3_skip_2_align, 0, aligns)
        tests.base_test.append_copies(g_3_skip_4_align, 1, aligns)
        tests.base_test.append_copies(g_3_skip_2_and_4_align, 0, aligns)
        tests.base_test.append_copies(g_3_skip_4_and_5_align, 2, aligns)
        tests.base_test.append_copies(g_4_all_align, 2, aligns)
        tests.base_test.append_copies(g_4_skip_2_align, 0, aligns)
        tests.base_test.append_copies(g_4_skip_4_align, 1, aligns)
        tests.base_test.append_copies(g_4_skip_2_and_4_align, 0, aligns)
        tests.base_test.append_copies(g_4_skip_4_and_5_align, 2, aligns)
        tests.base_test.write_sam(sample_1_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 1, aligns)
        tests.base_test.append_copies(early_inc_align, 1, aligns)
        tests.base_test.append_copies(skip_align, 1, aligns)
        tests.base_test.append_copies(skip_early_align, 1, aligns)
        tests.base_test.append_copies(g_2_inc_align, 1, aligns)
        tests.base_test.append_copies(g_2_early_inc_align, 1, aligns)
        tests.base_test.append_copies(g_2_skip_align, 2, aligns)
        tests.base_test.append_copies(g_3_all_align, 1, aligns)
        tests.base_test.append_copies(g_3_skip_2_align, 2, aligns)
        tests.base_test.append_copies(g_3_skip_4_align, 0, aligns)
        tests.base_test.append_copies(g_3_skip_2_and_4_align, 1, aligns)
        tests.base_test.append_copies(g_3_skip_4_and_5_align, 1, aligns)
        tests.base_test.append_copies(g_4_all_align, 1, aligns)
        tests.base_test.append_copies(g_4_skip_2_align, 2, aligns)
        tests.base_test.append_copies(g_4_skip_4_align, 0, aligns)
        tests.base_test.append_copies(g_4_skip_2_and_4_align, 1, aligns)
        tests.base_test.append_copies(g_4_skip_4_and_5_align, 1, aligns)
        tests.base_test.write_sam(sample_2_1_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

        aligns = list()
        tests.base_test.append_copies(inc_align, 1, aligns)
        tests.base_test.append_copies(early_inc_align, 1, aligns)
        tests.base_test.append_copies(skip_align, 1, aligns)
        tests.base_test.append_copies(skip_early_align, 1, aligns)
        tests.base_test.append_copies(g_2_inc_align, 1, aligns)
        tests.base_test.append_copies(g_2_early_inc_align, 1, aligns)
        tests.base_test.append_copies(g_2_skip_align, 2, aligns)
        tests.base_test.append_copies(g_3_all_align, 2, aligns)
        tests.base_test.append_copies(g_3_skip_2_align, 1, aligns)
        tests.base_test.append_copies(g_3_skip_4_align, 0, aligns)
        tests.base_test.append_copies(g_3_skip_2_and_4_align, 0, aligns)
        tests.base_test.append_copies(g_3_skip_4_and_5_align, 2, aligns)
        tests.base_test.append_copies(g_4_all_align, 2, aligns)
        tests.base_test.append_copies(g_4_skip_2_align, 1, aligns)
        tests.base_test.append_copies(g_4_skip_4_align, 0, aligns)
        tests.base_test.append_copies(g_4_skip_2_and_4_align, 0, aligns)
        tests.base_test.append_copies(g_4_skip_4_and_5_align, 2, aligns)
        tests.base_test.write_sam(sample_2_2_sam, self._chrs, aligns, next_i)
        next_i += len(aligns)

    def _check_output(self):
        asms_by_isoforms = self.parse_asm_definitions(self._event_dir)
        chr_id = '0'
        event_i = 0
        gene_1_asm_id = '{}_{}'.format(chr_id, event_i)
        event_i += 1
        gene_2_asm_id = '{}_{}'.format(chr_id, event_i)
        event_i += 1
        gene_3_asm_id = '{}_{}'.format(chr_id, event_i)
        event_i += 1
        gene_4_asm_id = '{}_{}'.format(chr_id, event_i)

        gene_1 = self._chrs[0].genes[0]
        gene_2 = self._chrs[0].genes[1]
        gene_3 = self._chrs[0].genes[2]
        gene_4 = self._chrs[0].genes[3]
        inc_isoform_id = gene_1.isoforms[1].id
        skip_isoform_id = gene_1.isoforms[2].id
        gene_1_inc_exons = ((81, 201), (401, 501), (701, 801))
        gene_1_skip_exons = ((101, 201), (701, 801))
        gene_1_isoforms = (gene_1_inc_exons, gene_1_skip_exons)

        g_2_inc_isoform_1_id = gene_2.isoforms[0].id
        g_2_inc_isoform_2_id = gene_2.isoforms[1].id
        g_2_inc_isoform_id = '_'.join(
            [g_2_inc_isoform_1_id, g_2_inc_isoform_2_id])
        g_2_skip_isoform_id = gene_2.isoforms[2].id
        gene_2_inc_exons = ((1081, 1201), (1401, 1501), (1701, 1801))
        gene_2_skip_exons = ((1101, 1201), (1701, 1801))
        gene_2_isoforms = (gene_2_inc_exons, gene_2_skip_exons)

        g_3_inc_isoform_id = gene_3.isoforms[0].id
        g_3_skip_2_isoform_id = gene_3.isoforms[1].id
        # g_3_skip_4_isoform_id = gene_3.isoforms[2].id
        # g_3_skip_2_4_isoform_id = gene_3.isoforms[3].id
        g_3_skip_4_5_isoform_id = gene_3.isoforms[4].id
        gene_3_inc_exons = ((2101, 2201), (2401, 2501), (2701, 2801),
                            (3001, 3101), (3301, 3401), (3601, 3701))
        gene_3_skip_2_exons = ((2101, 2201), (2701, 2801), (3001, 3101),
                               (3301, 3401), (3601, 3701))
        gene_3_skip_4_5_exons = ((2101, 2201), (2401, 2501), (2701, 2801),
                                 (3601, 3701))
        gene_3_isoforms = (gene_3_inc_exons, gene_3_skip_4_5_exons,
                           gene_3_skip_2_exons)

        g_4_inc_isoform_id = gene_4.isoforms[0].id
        g_4_skip_2_isoform_id = gene_4.isoforms[1].id
        # g_4_skip_4_isoform_id = gene_4.isoforms[2].id
        # g_4_skip_2_4_isoform_id = gene_4.isoforms[3].id
        g_4_skip_4_5_isoform_id = gene_4.isoforms[4].id
        gene_4_inc_exons = ((4101, 4201), (4401, 4501), (4701, 4801),
                            (5001, 5101), (5301, 5401), (5601, 5701))
        gene_4_skip_2_exons = ((4101, 4201), (4701, 4801), (5001, 5101),
                               (5301, 5401), (5601, 5701))
        gene_4_skip_4_5_exons = ((4101, 4201), (4401, 4501), (4701, 4801),
                                 (5601, 5701))
        gene_4_isoforms = (gene_4_inc_exons, gene_4_skip_4_5_exons,
                           gene_4_skip_2_exons)

        self.assertIn(gene_1_isoforms, asms_by_isoforms)
        gene_1_asm = asms_by_isoforms[gene_1_isoforms]
        self.assertIn(gene_2_isoforms, asms_by_isoforms)
        gene_2_asm = asms_by_isoforms[gene_2_isoforms]
        self.assertIn(gene_3_isoforms, asms_by_isoforms)
        gene_3_asm = asms_by_isoforms[gene_3_isoforms]
        self.assertIn(gene_4_isoforms, asms_by_isoforms)
        gene_4_asm = asms_by_isoforms[gene_4_isoforms]

        self.assertEqual(gene_1_asm['asm_id'], gene_1_asm_id)
        self.assertEqual(gene_1_asm['gene_id'], gene_1.id)
        self.assertEqual(gene_1_asm['chr'], self._chrs[0].name)
        self.assertEqual(gene_1_asm['isoform_ids'],
                         [inc_isoform_id, skip_isoform_id])
        self.assertEqual(gene_1_asm['start'], 'source')
        self.assertEqual(gene_1_asm['end'], 'sink')
        self.assertEqual(gene_1_asm['strand'], '+')

        self.assertEqual(gene_2_asm['asm_id'], gene_2_asm_id)
        self.assertEqual(gene_2_asm['gene_id'], gene_2.id)
        self.assertEqual(gene_2_asm['chr'], self._chrs[0].name)
        self.assertEqual(gene_2_asm['isoform_ids'],
                         [g_2_inc_isoform_id, g_2_skip_isoform_id])
        self.assertEqual(gene_2_asm['start'], 'source')
        self.assertEqual(gene_2_asm['end'], 'sink')
        self.assertEqual(gene_2_asm['strand'], '+')

        self.assertEqual(gene_3_asm['asm_id'], gene_3_asm_id)
        self.assertEqual(gene_3_asm['gene_id'], gene_3.id)
        self.assertEqual(gene_3_asm['chr'], self._chrs[0].name)
        self.assertEqual(gene_3_asm['isoform_ids'], [
            g_3_inc_isoform_id, g_3_skip_4_5_isoform_id, g_3_skip_2_isoform_id
        ])
        self.assertEqual(gene_3_asm['start'], 'source')
        self.assertEqual(gene_3_asm['end'], 'sink')
        self.assertEqual(gene_3_asm['strand'], '+')

        self.assertEqual(gene_4_asm['asm_id'], gene_4_asm_id)
        self.assertEqual(gene_4_asm['gene_id'], gene_4.id)
        self.assertEqual(gene_4_asm['chr'], self._chrs[0].name)
        self.assertEqual(gene_4_asm['isoform_ids'], [
            g_4_inc_isoform_id, g_4_skip_4_5_isoform_id, g_4_skip_2_isoform_id
        ])
        self.assertEqual(gene_4_asm['start'], 'source')
        self.assertEqual(gene_4_asm['end'], 'sink')
        self.assertEqual(gene_4_asm['strand'], '-')


if __name__ == '__main__':
    unittest.main(verbosity=2)
