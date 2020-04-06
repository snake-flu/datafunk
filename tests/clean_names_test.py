import os
import unittest
import filecmp

from datafunk.clean_names import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'clean_names')

class TestCleanNames(unittest.TestCase):
    def test_run(self):
        input_file = "%s/caseData.csv" %data_dir
        input_trait = "country/region"
        output_file = "%s/tmp.filtered.csv" %data_dir
        expected = "%s/filtered.csv" %data_dir
        clean_name(input_file, input_trait, output_file)
        self.assertTrue(filecmp.cmp(output_file, expected, shallow=False))
        os.unlink(output_file)
        os.unlink(output_file + ".log")