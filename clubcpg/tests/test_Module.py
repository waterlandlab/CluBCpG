import unittest
from clubcpg import ParseBam
from clubcpg.CalculateBinCoverage import CalculateCompleteBins
from clubcpg.ClusterReads import ClusterReads
import os
import pandas as pd
import numpy as np
from urllib.request import urlretrieve
from sklearn.cluster import DBSCAN


user_home = os.path.expanduser("~")
test_data_location = os.path.join(user_home, '.CluBCpG')


def download_data():
    print("Downloading test data...")
    test_data = ["TEST_DATA_A.bam", "TEST_DATA_B.bam", "TEST_DATA_A.bam.bai", "TEST_DATA_B.bam.bai", "TEST_BINS.csv"]
    for data in test_data:
        urlretrieve("https://s3.amazonaws.com/canthonyscott-mixture-analysis/{}".format(data),
                    os.path.join(test_data_location, data)) # todo should update this S3 bucketname
    return


if not os.path.exists(test_data_location):
    print("Creating test data location...")
    os.makedirs(test_data_location)
    download_data()


def check_data_exists(required_data):
    for data in required_data:
        if data not in os.listdir(test_data_location):
            download_data()

bamA = "TEST_DATA_A.bam"
bamB = "TEST_DATA_B.bam"
TEST_BINS="TEST_BINS.csv"
test_bin = "chr1_910700"
test_bin_bad = "chr1_500"

class TestParseBam(unittest.TestCase):
    """
    Test that the PraseBam class can function
    """

    # Make sure the test data is available
    def setUp(self):
        self.required_data = [bamA, 'TEST_DATA_B.bam', 'TEST_DATA_A.bam.bai', 'TEST_DATA_B.bam.bai']
        check_data_exists(self.required_data)
        self.parserA = ParseBam.BamFileReadParser(os.path.join(test_data_location, bamA), 20)

    # Make sure the unit tests work at all
    def test_cases_works(self):
        self.assertTrue(True, "Truth isn't truth -Rudy Giuliani, https://www.youtube.com/watch?v=CljsZ7lgbtw")

    def test_parser_loaded(self):
        self.assertIsInstance(self.parserA, ParseBam.BamFileReadParser, "Failed to load Parser object")

    def test_can_extract_reads(self):
        reads = self.parserA.parse_reads("chr1", 910600, 910700)
        self.assertEqual(len(reads), 118, "Reads could not be read from BAM file")
        matrix = self.parserA.create_matrix(reads)
        self.assertIsInstance(matrix, pd.DataFrame, "Reads failing to convert to data frame")
        self.assertEqual(matrix.shape, (118, 4), "Dataframe fails to be expected shape")


class TestCoverageCalculation(unittest.TestCase):
    """
    Test the features within the CoverageCalculation classes
    """

    def setUp(self):
        self.required_data = [bamA]
        check_data_exists(self.required_data)
        self.calc = CalculateCompleteBins(bam_file=os.path.join(test_data_location, bamA),
                                          bin_size=100,
                                          output_directory=os.path.join(test_data_location, 'TestCoverageCalculation'))

    def testCoverageCalcLoaded(self):
        self.assertIsInstance(self.calc, CalculateCompleteBins, "Failed to load module")

    def testCoverageCalc(self):
        b, matrix = self.calc.calculate_bin_coverage(test_bin)
        self.assertEqual(b, test_bin, "Bin failed to return ID correctly")
        self.assertIsInstance(matrix, pd.DataFrame, "matrix failed to return correctly")
        bad_result = self.calc.calculate_bin_coverage(test_bin_bad)
        self.assertIsNone(bad_result, "empty bin shouldn't have returned data")


class TestClustering(unittest.TestCase):
    """
    Test clustering functions
    """

    def setUp(self):
        self.required_data = [bamA, bamB, TEST_BINS]
        check_data_exists(self.required_data)
        self.parserA = ParseBam.BamFileReadParser(os.path.join(test_data_location, bamA), 20)
        self.parserB = ParseBam.BamFileReadParser(os.path.join(test_data_location, bamB), 20)
        self.reads_A = self.parserA.parse_reads("chr1", 910600, 910700)
        self.reads_B = self.parserB.parse_reads("chr1", 910600, 910700)
        self.matrix_A = self.parserA.create_matrix(self.reads_A).dropna()
        self.matrix_B = self.parserA.create_matrix(self.reads_B).dropna()
        self.cluster = ClusterReads(bamA, bamB, bins_file=TEST_BINS)

        self.matrix_A = self.matrix_A.copy()
        self.matrix_B = self.matrix_B.copy()
        labels_A = ['A'] * len(self.matrix_A)
        labels_B = ['B'] * len(self.matrix_B)
        self.matrix_A['input'] = labels_A
        self.matrix_B['input'] = labels_B
        self.full_matrix = pd.concat([self.matrix_A, self.matrix_B])
        self.data_to_cluster = np.matrix(self.full_matrix)[:, :-1]
        clf = DBSCAN(min_samples=2)
        labels = clf.fit_predict(self.data_to_cluster)
        self.full_matrix['class'] = labels
        self.filtered = self.cluster.filter_data_frame(self.full_matrix)


    def testExtractedReads(self):
        self.assertIsInstance(self.matrix_A, pd.DataFrame, "Failed to convert reads into matrix")
        self.assertIsInstance(self.matrix_B, pd.DataFrame, "Failed to convert reads into matrix")

    def testClustering(self):

        self.assertEqual(self.matrix_A.shape, (21, 5), "Failed to add labels to matrix")
        self.assertEqual(self.full_matrix.shape, (61, 6), "Failed to add cluster labels to matrix")
        self.assertIsInstance(self.cluster, ClusterReads, "Failed to load clustering class")


    def testPostClusterFiltering(self):
        self.assertEqual(self.filtered.shape, (56, 6), "Matrix did not filter correctly")
        self.assertEqual(len(self.cluster.get_common_matrices(self.filtered)), 4, "Failed to get common matrices")
        self.assertEqual(len(self.cluster.get_unique_matrices(self.filtered)), 2, "Failed to get unique matrices")


if __name__ == "__main__":
    unittest.main()