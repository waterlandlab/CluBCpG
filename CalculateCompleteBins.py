from ParseBam import BamFileReadParser
import sys
import os
import logging
from multiprocessing import Pool
import numpy as np


class CalculateCompleteBins:

    def __init__(self, bam_file: str, bin_size: int, output_directory: str, number_of_processors=1):
        """
        This class is initialized with a path to a bam file and a bin size
        :param bam_file: One of the BAM files for analysis to be performed
        :param bin_size: Size of the bins for the analysis, integer
        :number_of_processors: How many CPUs to use for parallel computation, default=1
        """
        self.input_bam_file = bam_file
        self.bin_size = int(bin_size)
        self.number_of_processors = int(number_of_processors)
        self.output_directory = output_directory

    def calculate_bin_coverage(self, bin):
        """
        Take a single bin, return a matrix
        :param bin: Bin should be passed as "Chr19_4343343"
        :return: pd.DataFrame with rows containing NaNs dropped
        """

        # Get reads from bam file
        parser = BamFileReadParser(self.input_bam_file, 20)
        # Split bin into parts
        chromosome, bin_location = bin.split("_")
        bin_location = int(bin_location)
        try:
            reads = parser.parse_reads(chromosome, (bin_location-self.bin_size), bin_location)
            matrix = parser.create_matrix(reads)

        except ValueError:
            # No reads are within this window, do nothing
            logging.info("No reads found for bin {}".format(bin))
            return None

        # convert to data_frame of 1s and 0s, drop rows with NaN
        matrix = matrix.dropna()

        return bin, matrix

    def get_chromosome_lengths(self):
        """
        Get dictionary containing lengths of the chromosomes. Uses bam file for reference
        :return: Dictionary of chromosome lengths, ex: {"chrX": 222222}
        """
        parser = BamFileReadParser(self.input_bam_file, 20)
        return dict(zip(parser.OpenBamFile.references, parser.OpenBamFile.lengths))

    def remove_scaffolds(self, chromosome_len_dict):
        """
        Return a dict containing only the standard chromosomes starting with "chr"
        :param chromosome_len_dict: A dict generated by get_chromosome_lenghts()
        :return: a dict containing only chromosomes starting with "chr"
        """
        new_dict = dict(chromosome_len_dict)
        for key in chromosome_len_dict.keys():
            if not key.startswith("chr"):
                new_dict.pop(key)

        return new_dict

    def generate_bins_list(self, chromosome_len_dict: dict):
        """
        Get a list of all bins according to desired bin size for all chromosomes in the passed dict
        :param chromosome_len_dict: A dict of chromosome length sizes from get_chromosome_lenghts, cleaned up by remove_scaffolds() if desired
        :return: list of all bins
        """
        all_bins = []
        for key, value in chromosome_len_dict.items():
            bins = list(np.arange(self.bin_size, value + self.bin_size, self.bin_size))
            bins = ["_".join([key, str(x)]) for x in bins]
            all_bins.extend(bins)

        return all_bins

    def analyze_bins(self, individual_chrom=None):
        # Get and clean dict of chromosome lenghts, convert to list of bins
        print("Getting Chromosome lengths from bam files...")
        chromosome_lengths = self.get_chromosome_lengths()
        chromosome_lengths = self.remove_scaffolds(chromosome_lengths)

        # If one chromosome was specified use only that chromosome
        if individual_chrom:
            new = dict()
            new[individual_chrom] = chromosome_lengths[individual_chrom]
            chromosome_lengths = new

        print("Generating bins for the entire genome...")
        bins_to_analyze = self.generate_bins_list(chromosome_lengths)

        # Set up for multiprocessing
        print("Beginning analysis of bins using {} processors".format(self.number_of_processors))
        pool = Pool(processes=self.number_of_processors)
        results = pool.map(self.calculate_bin_coverage, bins_to_analyze)
        logging.info("Analysis complete")

        # Write to output file
        output_file = os.path.join(self.output_directory, "CalculateCompleteBins_{}.csv".format(os.path.basename(self.input_bam_file)))

        with open(output_file, "w") as out:
            for result in results:
                out.write(result[0] + ",")
                out.write(str(result[1].shape[0]) + ",")
                out.write(str(result[1].shape[1]) + "\n")

        logging.info("Full read coverage analysis complete!")


if __name__ == "__main__":

    # Input params
    # todo add these as command line args using argparse
    input_bam_file = sys.argv[1]
    num_of_processors = sys.argv[2]
    if not num_of_processors:
        num_of_processors = 1

    # temp chromosome override for testing todo set this as an input arg
    chrom_of_interest = "chr19"

    log_file = "CalculateCompleteBins.{}.log".format(os.path.basename(input_bam_file))
    BASE_DIR = os.path.dirname(input_bam_file)

    logging.basicConfig(filename=os.path.join(BASE_DIR, log_file), level=logging.DEBUG)

    calc = CalculateCompleteBins(input_bam_file, 100, BASE_DIR, num_of_processors)

    calc.analyze_bins(chrom_of_interest)

