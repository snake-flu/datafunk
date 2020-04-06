import os
import unittest
import filecmp

from datafunk.merge_fasta import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'merge_fasta')

class TestMergeFasta(unittest.TestCase):
    def test_run(self):
        input_file = "%s/metadata.csv" %data_dir
        output_file = "%s/tmp.merged.fasta" %data_dir
        expected = "%s/merged.fasta" %data_dir
        merge_fasta(data_dir, input_file, output_file)
        self.assertTrue(filecmp.cmp(output_file, expected, shallow=False))
        os.unlink(output_file)
        os.unlink(output_file + ".log")