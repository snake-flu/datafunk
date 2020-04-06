import os
import unittest
import filecmp

from datafunk.remove_fasta import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'remove_fasta')

class TestRemoveFasta(unittest.TestCase):
    def test_run(self):
        input_file = "%s/test_remove.fasta" %data_dir
        filter_file = "%s/test_filter.txt" %data_dir
        output_file = "%s/tmp.filtered.fasta" %data_dir
        expected = "%s/filtered.fasta" %data_dir
        filter_dictionary = filter_list(filter_file)
        remove_fasta(input_file, filter_dictionary, output_file)
        self.assertTrue(filecmp.cmp(output_file, expected, shallow=False))
        os.unlink(output_file)

        self.assertTrue(filecmp.cmp(output_file + ".log", expected + ".log", shallow=False))
        os.unlink(output_file + ".log")