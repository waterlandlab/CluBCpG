from unittest import TestCase
from MixtureAnalysis import ParseBam


class TestParseBam(TestCase):

    def can_retrieve_reads(self):
        parser = ParseBam.BamFileReadParser("test_data/test_A.bam", 20)
        reads = parser.parse_reads("chr18", 37442711, 37450806)
        self.assertEquals(len(reads), 1793, "Failed to extract reads, is samtools installed?")
