import unittest

import tests.asm_detection.test as asm_detection_test
import tests.asm_read_counts.test as asm_read_counts_test
import tests.filtered_asms.test as filtered_asms_test
import tests.gene_isoforms.test as gene_isoforms_test
import tests.isoform_diffs.test as isoform_diffs_test
import tests.no_counts_replicate.test as no_counts_replicate_test
import tests.se_asm.test as se_asm_test
import tests.se_gene.test as se_gene_test
import tests.significance_filters.test as significance_filters_test


def build_test_suite():
    loader = unittest.defaultTestLoader
    suite = unittest.TestSuite()
    suite.addTest(loader.loadTestsFromModule(asm_detection_test))
    suite.addTest(loader.loadTestsFromModule(asm_read_counts_test))
    suite.addTest(loader.loadTestsFromModule(filtered_asms_test))
    suite.addTest(loader.loadTestsFromModule(gene_isoforms_test))
    suite.addTest(loader.loadTestsFromModule(isoform_diffs_test))
    suite.addTest(loader.loadTestsFromModule(no_counts_replicate_test))
    suite.addTest(loader.loadTestsFromModule(se_asm_test))
    suite.addTest(loader.loadTestsFromModule(se_gene_test))
    suite.addTest(loader.loadTestsFromModule(significance_filters_test))
    return suite


def main():
    runner = unittest.TextTestRunner(verbosity=2)
    suite = build_test_suite()
    runner.run(suite)


if __name__ == '__main__':
    main()
