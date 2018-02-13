from ParseBam import BamFileReadParser
import sys
import os
import logging
from multiprocessing import Pool
import numpy as np
import argparse
from collections import defaultdict


class CalculateCompleteBins:

    def __init__(self, bam_file: str, bin_size: int, output_directory: str, number_of_processors=1, mbias_read1_5=None, mbias_read1_3=None,
                 mbias_read2_5= None, mbias_read2_3=None):
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
        self.bins_no_reads = 0
        self.bins_yes_reads = 0

        self.mbias_read1_5 = mbias_read1_5
        self.mbias_read1_3 = mbias_read1_3
        self.mbias_read2_5 = mbias_read2_5
        self.mbias_read2_3 = mbias_read2_3

    def calculate_bin_coverage(self, bin):
        """
        Take a single bin, return a matrix
        :param bin: Bin should be passed as "Chr19_4343343"
        :return: pd.DataFrame with rows containing NaNs dropped
        """

        # Get reads from bam file
        parser = BamFileReadParser(self.input_bam_file, 20, self.mbias_read1_5, self.mbias_read1_3,
                                   self.mbias_read2_5, self.mbias_read2_3, no_overlap)
        # Split bin into parts
        chromosome, bin_location = bin.split("_")
        bin_location = int(bin_location)
        try:
            reads = parser.parse_reads(chromosome, bin_location-self.bin_size, bin_location)
            matrix = parser.create_matrix(reads)

        except ValueError:
            # No reads are within this window, do nothing
            # logging.info("No reads found for bin {}".format(bin))
            self.bins_no_reads += 1
            return None
        except:
            logging.error("Unknown error: {}".format(bin))
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

    @staticmethod
    def remove_scaffolds(chromosome_len_dict):
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
        Get a dict of lists of all bins according to desired bin size for all chromosomes in the passed dict
        :param chromosome_len_dict: A dict of chromosome length sizes from get_chromosome_lenghts, cleaned up by remove_scaffolds() if desired
        :return: dict with each key being a chromosome. ex: chr1
        """
        all_bins = defaultdict(list)
        for key, value in chromosome_len_dict.items():
            bins = list(np.arange(self.bin_size, value + self.bin_size, self.bin_size))
            bins = ["_".join([key, str(x)]) for x in bins]
            all_bins[key].extend(bins)

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
        print("Beginning analysis of bins using {} processors.".format(self.number_of_processors))
        print("This will take awhile.....")

        # Loop over bin dict and pool.map them individually
        final_results = []
        for key in bins_to_analyze.keys():
            print("Analyzing chromosome {}".format(key))
            pool = Pool(processes=self.number_of_processors)
            results = pool.map(self.calculate_bin_coverage, bins_to_analyze[key])
            final_results.extend(results)
            print("Finished chromosome {}".format(key))

        logging.info("Analysis complete")

        print("Complete.")
        print("Found {} bins without reads".format(self.bins_no_reads)) #todo idk if this works
        print("Found {} bins with reads. Writing these to a file.".format(len(final_results)))

        # Write to output file
        output_file = os.path.join(self.output_directory, "CalculateCompleteBins_{}.csv".format(os.path.basename(self.input_bam_file)))

        with open(output_file, "w") as out:
            for result in final_results:
                if result:
                    out.write(result[0] + ",")
                    out.write(str(result[1].shape[0]) + ",")
                    out.write(str(result[1].shape[1]) + "\n")

        logging.info("Full read coverage analysis complete!")


if __name__ == "__main__":

    def str2bool(v):
        if v.lower() == 'true':
            return True
        elif v.lower() == 'false':
            return False
        else:
            raise argparse.ArgumentTypeError("Boolean value expected.")

    # Input params
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-a", "--input_bam_A",
                            help="First Input bam file, coordinate sorted with index present")
    arg_parser.add_argument("-o", "--output_dir",
                            help="Output directory to save figures, defaults to bam file loaction")
    arg_parser.add_argument("-bin_size", help="Size of bins to extract and analyze, default=100", default=100)
    arg_parser.add_argument("-m", "--cluster_member_minimum",
                            help="Minimum number of members a cluster should have for it to be considered, default=4",
                            default=4)
    arg_parser.add_argument("-r", "--read_depth",
                            help="Minium number of reads covering all CpGs that the bins should have to analyze, default=10",
                            default=10)
    arg_parser.add_argument("-n", "--num_processors",
                            help="Number of processors to use for analysis, default=1",
                            default=1)
    arg_parser.add_argument("-chr", "--chromosome",
                            help="Chromosome to analyze, example: 'chr19', not required but encouraged, default=all chromosomes")

    arg_parser.add_argument("--read1_5", help="integer, read1 5' m-bias ignore bp, default=0", default=0)
    arg_parser.add_argument("--read1_3", help="integer, read1 3' m-bias ignore bp, default=0", default=0)
    arg_parser.add_argument("--read2_5", help="integer, read2 5' m-bias ignore bp, default=0", default=0)
    arg_parser.add_argument("--read2_3", help="integer, read2 3' m-bias ignore bp, default=0", default=0)
    arg_parser.add_argument("--no_overlap", help="bool, remove any overlap between paired reads and stitch"
                                                 " reads together when possible, default=True",
                            type=str2bool, const=True, default='True', nargs='?')

    args = arg_parser.parse_args()

    input_bam_file = args.input_bam_A
    num_of_processors = int(args.num_processors)
    bin_size = int(args.bin_size)
    min_cluster_members = int(args.cluster_member_minimum)
    read_depth_req = int(args.read_depth)
    no_overlap = args.no_overlap

    # Get the mbias inputs and adjust to work correctly, 0s should be converted to None
    mbias_read1_5 = int(args.read1_5)
    mbias_read1_3 = int(args.read1_3)
    mbias_read2_5 = int(args.read2_5)
    mbias_read2_3 = int(args.read2_3)


    if args.chromosome:
        chrom_of_interest = args.chromosome
    else:
        chrom_of_interest = None

    if args.output_dir:
        BASE_DIR = args.output_dir
    else:
        BASE_DIR = os.path.dirname(input_bam_file)

    # Create output dir if it doesnt exist
    if not os.path.exists(BASE_DIR):
        os.makedirs(BASE_DIR)

    log_file = os.path.join(BASE_DIR, "CalculateCompleteBins.{}.log".format(os.path.basename(input_bam_file)))

    logging.basicConfig(filename=os.path.join(BASE_DIR, log_file), level=logging.DEBUG)

    calc = CalculateCompleteBins(input_bam_file, 100, BASE_DIR, num_of_processors,
                                 mbias_read1_5, mbias_read1_3, mbias_read2_5, mbias_read2_3)

    logging.info("M bias inputs ignoring the following:\nread 1 5': {}bp\n"
                 "read1 3': {}bp\nread2 5: {}bp\nread2 3': {}bp".format(mbias_read1_5, mbias_read1_3, mbias_read2_5, mbias_read2_3))

    calc.analyze_bins(chrom_of_interest)

