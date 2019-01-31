import unittest
from MixtureAnalysis import ParseBam
from MixtureAnalysis.CalculateBinCoverage import CalculateCompleteBins
import os
import pandas as pd


user_home = os.path.expanduser("~")
test_data_location = os.path.join(user_home, '.MixtureAnalysis') # TODO change this once name is final


def download_data():
    print("Download data called")
    pass


if not os.path.exists(test_data_location):
    print("Creating test data location...")
    os.makedirs(test_data_location)
    download_data()


def check_data_exists(required_data):
    for data in required_data:
        if data not in os.listdir(test_data_location):
            download_data()

bamA = "TEST_DATA_A.bam"
test_bin = "chr1_910700"
test_bin_bad = ""

class TestParseBam(unittest.TestCase):

    def setUp(self):
        self.required_data = [bamA, 'TEST_DATA_B.bam', 'TEST_DATA_A.bam.bai', 'TEST_DATA_B.bam.bai']
        check_data_exists(self.required_data)
        self.parserA = ParseBam.BamFileReadParser(os.path.join(test_data_location, bamA), 20)

    def test_cases_works(self):
        self.assertTrue(True, "Truth isn't truth -Rudy Giuliani")

    def test_parser_loaded(self):
        self.assertIsInstance(self.parserA, ParseBam.BamFileReadParser, "Failed to load Parser object")

    def test_can_extract_reads(self):
        reads = self.parserA.parse_reads("chr1", 910600, 910700)
        self.assertEqual(len(reads), 118, "Reads could not be read from BAM file")
        matrix = self.parserA.create_matrix(reads)
        self.assertIsInstance(matrix, pd.DataFrame, "Reads failing to convert to data frame")
        self.assertEqual(matrix.shape, (118 ,4), "Dataframe fails to be expected shape")


class TestCoverageCalculation(unittest.TestCase):

    def setUp(self):
        self.required_data = [bamA]
        check_data_exists(self.required_data)
        self.calc = CalculateCompleteBins(bam_file=os.path.join(test_data_location, bamA), bin_size=100, output_directory=os.path.join(test_data_location, 'TestCoverageCalculation'))

    def testCoverageCalcLoaded(self):
        self.assertIsInstance(self.calc, CalculateCompleteBins, "Failed to load module")

    def testCoverageCalc(self):
        b, matrix = self.calc.calculate_bin_coverage(test_bin)
        self.assertEqual(b, test_bin, "Bin failed to return ID correctly")
        self.assertIsInstance(matrix, pd.DataFrame, "matrix failed to return correctly")


if __name__ == "__main__":
    unittest.main()
