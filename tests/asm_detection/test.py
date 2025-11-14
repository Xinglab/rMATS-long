import os
import os.path
import unittest

import tests.base_test


class ASMDetectionBaseTest(tests.base_test.BaseTest):
    def setUp(self):
        super().setUp()
        self._test_dir = os.path.dirname(__file__)
        self._input_dir = None
        self._out_dir = None
        self._log_dir = None
        self._gtf = None
        self._gtf_dir = None
        self._event_dir = None
        self._max_nodes = None
        self._max_paths = None

    def _set_input_out_log_with_prefix(self, prefix):
        self._input_dir = os.path.join(self._test_dir,
                                       '{}_input'.format(prefix))
        self._out_dir = os.path.join(self._test_dir, '{}_out'.format(prefix))
        self._log_dir = os.path.join(self._test_dir, '{}_log'.format(prefix))

    def _test(self):
        self._initialize_dirs()
        self._create_gtf()
        cleaned = os.path.join(self._out_dir, 'cleaned.gtf')
        self.run_clean_gtf(self._gtf, cleaned, self._log_dir)
        self._gtf_dir = self.run_organize_gtf(cleaned, self._out_dir,
                                              self._log_dir)
        self._event_dir = self.run_detect_events(self._gtf_dir,
                                                 self._out_dir,
                                                 self._log_dir,
                                                 max_nodes=self._max_nodes,
                                                 max_paths=self._max_paths)
        self._check_output()

    def _initialize_dirs(self):
        dirs = [self._input_dir, self._out_dir, self._log_dir]
        tests.base_test.remove_then_create_directories(dirs)

    def _create_gtf(self):
        self.fail('implement in subclass')

    def _check_output(self):
        self.fail('implement in subclass')


class SkippedExonTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'skipped_exon'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

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
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        gene.isoforms.append(isoform)
        # transcript 2 skips the exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        gene.isoforms.append(isoform)

        # minus strand
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1 includes the exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        gene.isoforms.append(isoform)
        # transcript 2 skips the exon
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 2)
        include_exons = ((301, 401), (501, 601), (701, 801))
        skip_exons = ((301, 401), (701, 801))
        key = tuple(sorted([include_exons, skip_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        include_exons = ((2301, 2401), (2501, 2601), (2701, 2801))
        skip_exons = ((2301, 2401), (2701, 2801))
        key = tuple(sorted([include_exons, skip_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')


class Alt5Test(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'alt_5'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        gene.isoforms.append(isoform)
        # transcript 2 has longer exon 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 605))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        gene.isoforms.append(isoform)
        # transcript 2 has exon 2 start earlier
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2495, 2600))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 2)
        long_exons = ((501, 606), (701, 801))
        short_exons = ((501, 601), (701, 801))
        key = tuple(sorted([long_exons, short_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], "alternative 5'-splice site usage")
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        long_exons = ((2301, 2401), (2496, 2601))
        short_exons = ((2301, 2401), (2501, 2601))
        key = tuple(sorted([long_exons, short_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], "alternative 5'-splice site usage")
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')


class Alt3Test(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'alt_3'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        gene.isoforms.append(isoform)
        # transcript 2: exon 2 starts sooner
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(495, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        gene.isoforms.append(isoform)
        # transcript 2: exon 2 ends later
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2605))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 2)
        long_exons = ((301, 401), (496, 601))
        short_exons = ((301, 401), (501, 601))
        key = tuple(sorted([long_exons, short_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], "alternative 3'-splice site usage")
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        long_exons = ((2501, 2606), (2701, 2801))
        short_exons = ((2501, 2601), (2701, 2801))
        key = tuple(sorted([long_exons, short_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], "alternative 3'-splice site usage")
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')


class MutuallyExclusiveTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'mxe'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1: exon 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)
        # transcript 2: exon 3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1: exon 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        isoform.exons.append(tests.base_test.Region(3100, 3200))
        gene.isoforms.append(isoform)
        # transcript 2: exon 3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        isoform.exons.append(tests.base_test.Region(3100, 3200))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 2)
        first_exons = ((301, 401), (501, 601), (901, 1001))
        second_exons = ((301, 401), (701, 801), (901, 1001))
        key = tuple(sorted([first_exons, second_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'mutually exclusive exons')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        first_exons = ((2301, 2401), (2701, 2801), (2901, 3001))
        second_exons = ((2301, 2401), (2501, 2601), (2901, 3001))
        key = tuple(sorted([first_exons, second_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'mutually exclusive exons')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')


class RetainedIntronTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'ri'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        gene.isoforms.append(isoform)
        # transcript 2: retained
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        gene.isoforms.append(isoform)
        # transcript 2: retained
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2600))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 2)
        spliced_exons = ((301, 401), (501, 601))
        retained_exons = ((301, 601), )
        key = tuple(sorted([spliced_exons, retained_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'intron retention')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        spliced_exons = ((2301, 2401), (2501, 2601))
        retained_exons = ((2301, 2601), )
        key = tuple(sorted([spliced_exons, retained_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'intron retention')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')


class AltFirstTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'alt_first'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1: exon 0
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(500, 600))
        gene.isoforms.append(isoform)
        # transcript 2: exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 600))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1: exon 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        gene.isoforms.append(isoform)
        # transcript 2: exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 2)
        first_exons = ((101, 201), (501, 601))
        second_exons = ((301, 401), (501, 601))
        key = tuple(sorted([first_exons, second_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'alternative first exon usage')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        first_exons = ((2101, 2201), (2501, 2601))
        second_exons = ((2101, 2201), (2301, 2401))
        key = tuple(sorted([first_exons, second_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'alternative first exon usage')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')


class AltLastTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'alt_last'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1: exon 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(500, 600))
        gene.isoforms.append(isoform)
        # transcript 2: exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1: exon 0
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        gene.isoforms.append(isoform)
        # transcript 2: exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 2)
        last_exons = ((101, 201), (501, 601))
        second_exons = ((101, 201), (301, 401))
        key = tuple(sorted([last_exons, second_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'alternative last exon usage')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        last_exons = ((2101, 2201), (2501, 2601))
        second_exons = ((2301, 2401), (2501, 2601))
        key = tuple(sorted([last_exons, second_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'alternative last exon usage')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')


class ComplexTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'complex'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)
        # transcript 2: skip 2 and 3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        isoform.exons.append(tests.base_test.Region(3100, 3200))
        gene.isoforms.append(isoform)
        # transcript 2: skip 2 and 3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        isoform.exons.append(tests.base_test.Region(3100, 3200))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 2)
        include_exons = ((301, 401), (501, 601), (701, 801), (901, 1001))
        skip_exons = ((301, 401), (901, 1001))
        key = tuple(sorted([include_exons, skip_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        include_exons = ((2301, 2401), (2501, 2601), (2701, 2801), (2901,
                                                                    3001))
        skip_exons = ((2301, 2401), (2901, 3001))
        key = tuple(sorted([include_exons, skip_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')


class AltStartEndTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'alt_start_end'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 600))
        gene.isoforms.append(isoform)
        # transcript 2: exon 0 starts later
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(105, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 600))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        isoform.exons.append(tests.base_test.Region(1500, 1600))
        gene.isoforms.append(isoform)
        # transcript 2: exon 2 ends sooner
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(1100, 1200))
        isoform.exons.append(tests.base_test.Region(1300, 1400))
        isoform.exons.append(tests.base_test.Region(1500, 1595))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_3'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        gene.isoforms.append(isoform)
        # transcript 2: exon 2 ends sooner
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2595))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_4'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(3100, 3200))
        isoform.exons.append(tests.base_test.Region(3300, 3400))
        isoform.exons.append(tests.base_test.Region(3500, 3600))
        gene.isoforms.append(isoform)
        # transcript 2: exon 0 starts later
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(3105, 3200))
        isoform.exons.append(tests.base_test.Region(3300, 3400))
        isoform.exons.append(tests.base_test.Region(3500, 3600))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_5'

        # transcript 1: exon 0 starts later
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(4103, 4200))
        isoform.exons.append(tests.base_test.Region(4300, 4400))
        isoform.exons.append(tests.base_test.Region(4500, 4600))
        gene.isoforms.append(isoform)
        # transcript 2: exon 0 starts sooner
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(4097, 4200))
        isoform.exons.append(tests.base_test.Region(4300, 4400))
        isoform.exons.append(tests.base_test.Region(4500, 4600))
        gene.isoforms.append(isoform)
        # transcript 3: exon 2 ends later
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(4100, 4200))
        isoform.exons.append(tests.base_test.Region(4300, 4400))
        isoform.exons.append(tests.base_test.Region(4500, 4603))
        gene.isoforms.append(isoform)
        # transcript 4: exon 2 ends sooner
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(4100, 4200))
        isoform.exons.append(tests.base_test.Region(4300, 4400))
        isoform.exons.append(tests.base_test.Region(4500, 4597))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_6'

        # transcript 1: exon 0 starts later
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(5103, 5200))
        isoform.exons.append(tests.base_test.Region(5300, 5400))
        isoform.exons.append(tests.base_test.Region(5500, 5600))
        gene.isoforms.append(isoform)
        # transcript 2: exon 0 starts sooner
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(5097, 5200))
        isoform.exons.append(tests.base_test.Region(5300, 5400))
        isoform.exons.append(tests.base_test.Region(5500, 5600))
        gene.isoforms.append(isoform)
        # transcript 3: exon 2 ends later
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(5100, 5200))
        isoform.exons.append(tests.base_test.Region(5300, 5400))
        isoform.exons.append(tests.base_test.Region(5500, 5603))
        gene.isoforms.append(isoform)
        # transcript 4: exon 2 ends sooner
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(5100, 5200))
        isoform.exons.append(tests.base_test.Region(5300, 5400))
        isoform.exons.append(tests.base_test.Region(5500, 5597))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_7'

        # transcript 1: includes exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(6100, 6200))
        isoform.exons.append(tests.base_test.Region(6300, 6400))
        isoform.exons.append(tests.base_test.Region(6500, 6600))
        gene.isoforms.append(isoform)
        # transcript 2: +3 start, skips exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(6103, 6200))
        isoform.exons.append(tests.base_test.Region(6500, 6600))
        gene.isoforms.append(isoform)
        # transcript 3: -3 start, skips exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(6097, 6200))
        isoform.exons.append(tests.base_test.Region(6500, 6600))
        gene.isoforms.append(isoform)
        # transcript 4: +6 start, retains 1st intron
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(6106, 6400))
        isoform.exons.append(tests.base_test.Region(6500, 6600))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_8'

        # transcript 1: includes exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(7100, 7200))
        isoform.exons.append(tests.base_test.Region(7300, 7400))
        isoform.exons.append(tests.base_test.Region(7500, 7600))
        gene.isoforms.append(isoform)
        # transcript 2: -3 end, skips exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(7100, 7200))
        isoform.exons.append(tests.base_test.Region(7500, 7597))
        gene.isoforms.append(isoform)
        # transcript 3: +3 end, skips exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(7100, 7200))
        isoform.exons.append(tests.base_test.Region(7500, 7603))
        gene.isoforms.append(isoform)
        # transcript 4: -6 end, retains 2nd intron
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(7100, 7200))
        isoform.exons.append(tests.base_test.Region(7300, 7594))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_9'

        # transcript 1: includes exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(8100, 8200))
        isoform.exons.append(tests.base_test.Region(8300, 8400))
        isoform.exons.append(tests.base_test.Region(8500, 8600))
        gene.isoforms.append(isoform)
        # transcript 2: -3 end, skips exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(8100, 8200))
        isoform.exons.append(tests.base_test.Region(8500, 8597))
        gene.isoforms.append(isoform)
        # transcript 3: +3 end, skips exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(8100, 8200))
        isoform.exons.append(tests.base_test.Region(8500, 8603))
        gene.isoforms.append(isoform)
        # transcript 4: -6 end, retains 2nd intron
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(8100, 8200))
        isoform.exons.append(tests.base_test.Region(8300, 8594))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_10'

        # transcript 1: includes exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(9100, 9200))
        isoform.exons.append(tests.base_test.Region(9300, 9400))
        isoform.exons.append(tests.base_test.Region(9500, 9600))
        gene.isoforms.append(isoform)
        # transcript 2: +3 start, skips exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(9103, 9200))
        isoform.exons.append(tests.base_test.Region(9500, 9600))
        gene.isoforms.append(isoform)
        # transcript 3: -3 start, skips exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(9097, 9200))
        isoform.exons.append(tests.base_test.Region(9500, 9600))
        gene.isoforms.append(isoform)
        # transcript 4: +6 start, retains 1st intron
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(9106, 9400))
        isoform.exons.append(tests.base_test.Region(9500, 9600))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_11'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(10100, 10200))
        isoform.exons.append(tests.base_test.Region(10300, 10400))
        gene.isoforms.append(isoform)
        # transcript 2: exon 0 +3 start
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(10103, 10200))
        isoform.exons.append(tests.base_test.Region(10300, 10400))
        gene.isoforms.append(isoform)
        # transcript 3: exon 0 +3 start, +3 end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(10103, 10203))
        isoform.exons.append(tests.base_test.Region(10300, 10400))
        gene.isoforms.append(isoform)
        # transcript 4: exon 0 -3 start, +3 end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(10097, 10203))
        isoform.exons.append(tests.base_test.Region(10300, 10400))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_12'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(11100, 11200))
        isoform.exons.append(tests.base_test.Region(11300, 11400))
        gene.isoforms.append(isoform)
        # transcript 2: exon 1 -3 end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(11100, 11200))
        isoform.exons.append(tests.base_test.Region(11300, 11397))
        gene.isoforms.append(isoform)
        # transcript 3: exon 1 -3 start, -3 end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(11100, 11200))
        isoform.exons.append(tests.base_test.Region(11297, 11397))
        gene.isoforms.append(isoform)
        # transcript 4: exon 1 -3 start, +3 end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(11100, 11200))
        isoform.exons.append(tests.base_test.Region(11297, 11403))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_13'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(12100, 12200))
        isoform.exons.append(tests.base_test.Region(12300, 12400))
        gene.isoforms.append(isoform)
        # transcript 2: exon 1 -3 end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(12100, 12200))
        isoform.exons.append(tests.base_test.Region(12300, 12397))
        gene.isoforms.append(isoform)
        # transcript 3: exon 1 -3 start, -3 end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(12100, 12200))
        isoform.exons.append(tests.base_test.Region(12297, 12397))
        gene.isoforms.append(isoform)
        # transcript 4: exon 1 -3 start, +3 end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(12100, 12200))
        isoform.exons.append(tests.base_test.Region(12297, 12403))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_14'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(13100, 13200))
        isoform.exons.append(tests.base_test.Region(13300, 13400))
        gene.isoforms.append(isoform)
        # transcript 2: exon 0 +3 start
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(13103, 13200))
        isoform.exons.append(tests.base_test.Region(13300, 13400))
        gene.isoforms.append(isoform)
        # transcript 3: exon 0 +3 start, +3 end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(13103, 13203))
        isoform.exons.append(tests.base_test.Region(13300, 13400))
        gene.isoforms.append(isoform)
        # transcript 4: exon 0 -3 start, +3 end
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(13097, 13203))
        isoform.exons.append(tests.base_test.Region(13300, 13400))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 16)
        include_exons = ((6104, 6201), (6301, 6401), (6501, 6601))
        skip_exons = ((6104, 6201), (6501, 6601))
        key = tuple(sorted([include_exons, skip_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_7')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        start_with_intron = ((6107, 6401), )
        start_with_exon = ((6098, 6201), (6301, 6401))
        key = tuple(sorted([start_with_intron, start_with_exon]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_7')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        start_with_intron = ((6107, 6401), (6501, 6601))
        start_with_exon = ((6098, 6201), (6301, 6401), (6501, 6601))
        start_with_skip = ((6098, 6201), (6501, 6601))
        key = tuple(
            sorted([start_with_intron, start_with_exon, start_with_skip]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_7')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        include_exons = ((7101, 7201), (7301, 7401), (7501, 7598))
        skip_exons = ((7101, 7201), (7501, 7598))
        key = tuple(sorted([include_exons, skip_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_8')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        end_with_intron = ((7301, 7595), )
        end_with_exon = ((7301, 7401), (7501, 7604))
        key = tuple(sorted([end_with_intron, end_with_exon]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_8')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        end_with_intron = ((7101, 7201), (7301, 7595))
        end_with_exon = ((7101, 7201), (7301, 7401), (7501, 7604))
        end_with_skip = ((7101, 7201), (7501, 7604))
        key = tuple(sorted([end_with_intron, end_with_exon, end_with_skip]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_8')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        include_exons = ((8101, 8201), (8301, 8401), (8501, 8598))
        skip_exons = ((8101, 8201), (8501, 8598))
        key = tuple(sorted([include_exons, skip_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_9')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        start_with_intron = ((8301, 8595), )
        start_with_exon = ((8301, 8401), (8501, 8604))
        key = tuple(sorted([start_with_intron, start_with_exon]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_9')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        start_with_intron = ((8101, 8201), (8301, 8595))
        start_with_exon = ((8101, 8201), (8301, 8401), (8501, 8604))
        start_with_skip = ((8101, 8201), (8501, 8604))
        key = tuple(
            sorted([start_with_intron, start_with_exon, start_with_skip]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_9')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        include_exons = ((9104, 9201), (9301, 9401), (9501, 9601))
        skip_exons = ((9104, 9201), (9501, 9601))
        key = tuple(sorted([include_exons, skip_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_10')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        end_with_intron = ((9107, 9401), )
        end_with_exon = ((9098, 9201), (9301, 9401))
        key = tuple(sorted([end_with_intron, end_with_exon]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_10')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        end_with_intron = ((9107, 9401), (9501, 9601))
        end_with_exon = ((9098, 9201), (9301, 9401), (9501, 9601))
        end_with_skip = ((9098, 9201), (9501, 9601))
        key = tuple(sorted([end_with_intron, end_with_exon, end_with_skip]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_10')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        long_isoform = ((10104, 10204), (10301, 10401))
        short_isoform = ((10104, 10201), (10301, 10401))
        key = tuple(sorted([long_isoform, short_isoform]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], "alternative 5'-splice site usage")
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_11')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        long_isoform = ((11101, 11201), (11298, 11398))
        short_isoform = ((11101, 11201), (11301, 11398))
        key = tuple(sorted([long_isoform, short_isoform]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], "alternative 5'-splice site usage")
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_12')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        long_isoform = ((12101, 12201), (12298, 12398))
        short_isoform = ((12101, 12201), (12301, 12398))
        key = tuple(sorted([long_isoform, short_isoform]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], "alternative 3'-splice site usage")
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_13')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        long_isoform = ((13104, 13204), (13301, 13401))
        short_isoform = ((13104, 13201), (13301, 13401))
        key = tuple(sorted([long_isoform, short_isoform]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], "alternative 3'-splice site usage")
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_14')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')


class MixedEventsTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'mixed'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        gene.isoforms.append(isoform)
        # transcript 2: skip exon 1, exon 4 starts sooner
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(895, 1000))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        gene.isoforms.append(isoform)
        # transcript 2: skip exon 3, exon 0 ends later
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2205))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 4)
        include_exons = ((101, 201), (301, 401), (501, 601))
        skip_exons = ((101, 201), (501, 601))
        key = tuple(sorted([include_exons, skip_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        long_exons = ((701, 801), (896, 1001))
        short_exons = ((701, 801), (901, 1001))
        key = tuple(sorted([long_exons, short_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], "alternative 3'-splice site usage")
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        include_exons = ((2501, 2601), (2701, 2801), (2901, 3001))
        skip_exons = ((2501, 2601), (2901, 3001))
        key = tuple(sorted([include_exons, skip_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        long_exons = ((2101, 2206), (2301, 2401))
        short_exons = ((2101, 2201), (2301, 2401))
        key = tuple(sorted([long_exons, short_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], "alternative 3'-splice site usage")
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')


class SeveralIsoformsTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'several'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        gene.isoforms.append(isoform)
        # transcript 2: skip 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        gene.isoforms.append(isoform)
        # transcript 3: skip 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(700, 800))
        gene.isoforms.append(isoform)
        # transcript 4: skip 1 and 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(700, 800))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        gene.isoforms.append(isoform)
        # transcript 2: skip 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        gene.isoforms.append(isoform)
        # transcript 3: skip 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        gene.isoforms.append(isoform)
        # transcript 4: skip 1 and 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_4'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 6)
        inc1_exons = ((101, 201), (301, 401), (501, 601))
        skip1_exons = ((101, 201), (501, 601))
        key = tuple(sorted([inc1_exons, skip1_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        inc2_exons = ((301, 401), (501, 601), (701, 801))
        skip2_exons = ((301, 401), (701, 801))
        key = tuple(sorted([inc2_exons, skip2_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        all_exons = ((101, 201), (301, 401), (501, 601), (701, 801))
        partial_skip1_exons = ((101, 201), (301, 401), (701, 801))
        partial_skip2_exons = ((101, 201), (501, 601), (701, 801))
        skip3_exons = ((101, 201), (701, 801))
        key = tuple(
            sorted([
                all_exons, partial_skip1_exons, partial_skip2_exons,
                skip3_exons
            ]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        inc1_exons = ((2301, 2401), (2501, 2601), (2701, 2801))
        skip1_exons = ((2301, 2401), (2701, 2801))
        key = tuple(sorted([inc1_exons, skip1_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        inc2_exons = ((2101, 2201), (2301, 2401), (2501, 2601))
        skip2_exons = ((2101, 2201), (2501, 2601))
        key = tuple(sorted([inc2_exons, skip2_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        all_exons = ((2101, 2201), (2301, 2401), (2501, 2601), (2701, 2801))
        partial_skip1_exons = ((2101, 2201), (2301, 2401), (2701, 2801))
        partial_skip2_exons = ((2101, 2201), (2501, 2601), (2701, 2801))
        skip3_exons = ((2101, 2201), (2701, 2801))
        key = tuple(
            sorted([
                all_exons, partial_skip1_exons, partial_skip2_exons,
                skip3_exons
            ]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')


class IsoformInMultipleGenesTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'multiple_gene_isoform'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1,gene_2'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1: skip exon 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1: skip exon 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 2)
        include_exons = ((101, 201), (301, 401), (501, 601))
        skip_exons = ((101, 201), (501, 601))
        key = tuple(sorted([include_exons, skip_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        include_exons = ((301, 401), (501, 601), (701, 801))
        skip_exons = ((301, 401), (701, 801))
        key = tuple(sorted([include_exons, skip_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')


class EachExonSkippedTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'each_exon'
        self._set_input_out_log_with_prefix(self._prefix)
        self._max_paths = 2000
        self.num_exons = 15

    def test(self):
        start_time = tests.base_test.start_timer()
        self._test()
        self.assert_seconds_at_most(start_time, 30)

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'
        # Gene with 15 exons.
        # 15 transcripts that each skip 1 exon.
        # ^ -> 1 -> 2 -----> 5 -> 6 -----> 9 -> 10 ...
        #  \         \      /      \      /        ...
        #   --------> 3 -> 4 -----> 7 -> 8 ----->  ...
        exons = list()
        for i in range(self.num_exons):
            start = (i * 200) + 100
            end = start + 100
            exons.append(tests.base_test.Region(start, end))

        for i in range(self.num_exons):
            isoform_exons = exons[:]
            isoform_exons.pop(i)
            isoform = tests.base_test.Isoform()
            isoform.id = '{}_t_{}'.format(gene.id, i)
            isoform.strand = '+'
            isoform.exons = isoform_exons
            gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        # ASMs:
        # n: from source to each start coord (except 1st exon) (also to sink)
        # n-1: from 1st exon end to each later start coord (except 2nd exon) (also to sink)
        # ...
        # 1: from 2nd to last exon end to sink
        # n + n-1 + n-2 + .. + 1 -> n(n+1)/2
        total_asms = self.num_exons * (self.num_exons + 1) / 2
        self.assertEqual(len(asms), total_asms)
        type_counts = dict()
        for asm in asms.values():
            event_type = asm['event_type']
            old_count = type_counts.get(event_type, 0)
            type_counts[event_type] = old_count + 1

        self.assertEqual(len(type_counts), 2)
        # There is a skipping event for each exon except the first and last.
        # (Skipping the first or last exon doesn't fit the skipped exon definition)
        # All other ASMs are complex splicing
        num_skipping = self.num_exons - 2
        self.assertEqual(type_counts['exon skipping'], num_skipping)
        self.assertEqual(type_counts['complex splicing'],
                         total_asms - num_skipping)


class MostExonsSkippedTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'most_exons'
        self._set_input_out_log_with_prefix(self._prefix)
        self._max_paths = 100
        self.num_exons = 1000
        self.share_every_ith = 10
        self.num_shared = self.num_exons / self.share_every_ith

    def test(self):
        start_time = tests.base_test.start_timer()
        self._test()
        self.assert_seconds_at_most(start_time, 30)

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'
        # Gene with 1000 exons.
        # Every 10th exon is shared in all isoforms.
        # Other exons have a skipping isoform.
        # ^ -> 1 -> 2 -----> 5 -> 6 -----> ... --> 19 -> 20 --------> 23 -> ...
        #            \      /      \       ...    /        \         /      ...
        #             3 -> 4 -----> 7 -> 8 ... --/          21 -> 22  --->  ...
        exons = list()
        for i in range(self.num_exons):
            start = (i * 200) + 100
            end = start + 100
            exons.append(tests.base_test.Region(start, end))

        for i in range(self.num_exons):
            if i % self.share_every_ith == 0:
                continue

            isoform_exons = exons[:]
            isoform_exons.pop(i)
            isoform = tests.base_test.Isoform()
            isoform.id = '{}_t_{}'.format(gene.id, i)
            isoform.strand = '+'
            isoform.exons = isoform_exons
            gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        # ASMs:
        # Each block of (i-1) exons in between shared exons:
        #   (i-1)*(i)/2
        # There are (n/i) blocks
        num_blocks = self.num_shared
        total_asms = (num_blocks *
                      (self.share_every_ith - 1) * self.share_every_ith / 2)
        self.assertEqual(len(asms), total_asms)

        type_counts = dict()
        for asm in asms.values():
            event_type = asm['event_type']
            old_count = type_counts.get(event_type, 0)
            type_counts[event_type] = old_count + 1

        self.assertEqual(len(type_counts), 2)
        # There is a skipping event for each exon except:
        # * the n/i which are shared
        # * the last exon which is not shared but doesn't fit the SE definition
        # All other ASMs are complex splicing
        num_skipping = (self.num_exons - self.num_shared) - 1
        self.assertEqual(type_counts['exon skipping'], num_skipping)
        self.assertEqual(type_counts['complex splicing'],
                         total_asms - num_skipping)


class TwoLongIsoformsTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'two_long'
        self._set_input_out_log_with_prefix(self._prefix)
        self._max_nodes = 2000
        self.num_exons = 1000

    def test(self):
        start_time = tests.base_test.start_timer()
        self._test()
        self.assert_seconds_at_most(start_time, 30)

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'
        # Gene with 1000 exons.
        # 2 transcripts which split the exons.
        # The first and last exon are shared
        # ^ -> 1 -> 2 -----> 5 -> 6 -----> 9 -> 10 ... ---
        #            \                             ...    \
        #             3 -> 4 -----> 7 -> 8 ----->  ... --> last_exon -> $
        exons = list()
        for i in range(self.num_exons):
            start = (i * 200) + 100
            end = start + 100
            exons.append(tests.base_test.Region(start, end))

        isoform_1_exons = [exons[0]]
        isoform_2_exons = [exons[0]]
        for i in range(1, self.num_exons - 1):
            if i % 2 == 0:
                isoform_1_exons.append(exons[i])
            else:
                isoform_2_exons.append(exons[i])

        isoform_1_exons.append(exons[self.num_exons - 1])
        isoform_2_exons.append(exons[self.num_exons - 1])

        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons = isoform_1_exons
        gene.isoforms.append(isoform)

        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons = isoform_2_exons
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 1)
        asm = list(asms.values())[0]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(len(asm['isoforms']), 2)
        self.assertEqual(asm['is_strict'], 'True')


class OverlappingSkipsTest(ASMDetectionBaseTest):
    def setUp(self):
        super().setUp()
        self._prefix = 'overlap_skip'
        self._set_input_out_log_with_prefix(self._prefix)

    def test(self):
        self._test()

    def _create_gtf(self):
        chrom = tests.base_test.Chromosome()
        chrom.name = 'chr_1'
        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_1'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(500, 600))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        gene.isoforms.append(isoform)
        # transcript 2: skip 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(300, 400))
        isoform.exons.append(tests.base_test.Region(700, 800))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        gene.isoforms.append(isoform)
        # transcript 3: skip 1,2,3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '+'
        isoform.exons.append(tests.base_test.Region(100, 200))
        isoform.exons.append(tests.base_test.Region(900, 1000))
        gene.isoforms.append(isoform)

        gene = tests.base_test.Gene()
        chrom.genes.append(gene)
        gene.id = 'gene_2'

        # transcript 1
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_1'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2500, 2600))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        gene.isoforms.append(isoform)
        # transcript 2: skip 2
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_2'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2300, 2400))
        isoform.exons.append(tests.base_test.Region(2700, 2800))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        gene.isoforms.append(isoform)
        # transcript 3: skip 1,2,3
        isoform = tests.base_test.Isoform()
        isoform.id = '{}_t_3'.format(gene.id)
        isoform.strand = '-'
        isoform.exons.append(tests.base_test.Region(2100, 2200))
        isoform.exons.append(tests.base_test.Region(2900, 3000))
        gene.isoforms.append(isoform)

        self._gtf = os.path.join(self._input_dir, 'test.gtf')
        chrs = [chrom]
        tests.base_test.write_gtf(self._gtf, chrs)

    def _check_output(self):
        asms = self.parse_asm_definitions(self._event_dir)
        self.assertEqual(len(asms), 4)
        inc2_exons = ((301, 401), (501, 601), (701, 801))
        skip2_exons = ((301, 401), (701, 801))
        key = tuple(sorted([inc2_exons, skip2_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        inc_all_exons = ((101, 201), (301, 401), (501, 601), (701, 801),
                         (901, 1001))
        skip2_all_exons = ((101, 201), (301, 401), (701, 801), (901, 1001))
        skip123_all_exons = ((101, 201), (901, 1001))
        key = tuple(sorted([inc_all_exons, skip2_all_exons,
                            skip123_all_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '+')
        self.assertEqual(asm['gene_id'], 'gene_1')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')

        inc2_exons = ((2301, 2401), (2501, 2601), (2701, 2801))
        skip2_exons = ((2301, 2401), (2701, 2801))
        key = tuple(sorted([inc2_exons, skip2_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'exon skipping')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'False')

        inc_all_exons = ((2101, 2201), (2301, 2401), (2501, 2601),
                         (2701, 2801), (2901, 3001))
        skip2_all_exons = ((2101, 2201), (2301, 2401), (2701, 2801), (2901,
                                                                      3001))
        skip123_all_exons = ((2101, 2201), (2901, 3001))
        key = tuple(sorted([inc_all_exons, skip2_all_exons,
                            skip123_all_exons]))
        self.assertIn(key, asms)
        asm = asms[key]
        self.assertEqual(asm['event_type'], 'complex splicing')
        self.assertEqual(asm['chr'], 'chr_1')
        self.assertEqual(asm['strand'], '-')
        self.assertEqual(asm['gene_id'], 'gene_2')
        self.assertEqual(asm['isoforms'], list(key))
        self.assertEqual(asm['is_strict'], 'True')


if __name__ == '__main__':
    unittest.main(verbosity=2)
