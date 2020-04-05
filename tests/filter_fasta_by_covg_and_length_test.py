import os
import unittest
import filecmp

from datafunk.filter_fasta_by_covg_and_length import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'filter_fasta_by_covg_and_length')

class TestFilterFastaByCovgAndLength(unittest.TestCase):
    def test_threshold_70(self):
        infile = "%s/test.fasta" %data_dir
        outfile = "%s/tmp.threshold_70.fasta" %data_dir
        threshold = 70
        expected = "%s/expected_threshold_70.fasta" %data_dir
        filter_sequences(infile, outfile, min_covg=threshold)
        self.assertTrue(filecmp.cmp(outfile, expected, shallow=False))
        os.unlink(outfile)

    def test_threshold_80(self):
        infile = "%s/test.fasta" %data_dir
        outfile = "%s/tmp.threshold_80.fasta" %data_dir
        threshold = 80
        expected = "%s/expected_threshold_80.fasta" %data_dir
        filter_sequences(infile, outfile, min_covg=threshold)
        self.assertTrue(filecmp.cmp(outfile, expected, shallow=False))
        os.unlink(outfile)

    def test_threshold_90(self):
        infile = "%s/test.fasta" %data_dir
        outfile = "%s/tmp.threshold_90.fasta" %data_dir
        threshold = 90
        expected = "%s/expected_threshold_90.fasta" %data_dir
        filter_sequences(infile, outfile, min_covg=threshold)
        self.assertTrue(filecmp.cmp(outfile, expected, shallow=False))
        os.unlink(outfile)

    def test_threshold_100(self):
        infile = "%s/test.fasta" %data_dir
        outfile = "%s/tmp.threshold_100.fasta" %data_dir
        threshold = 100
        expected = "%s/expected_threshold_100.fasta" %data_dir
        filter_sequences(infile, outfile, min_covg=threshold)
        self.assertTrue(filecmp.cmp(outfile, expected, shallow=False))
        os.unlink(outfile)

