from ParseBam import BamFileReadParser
import sys
import os
import logging
from multiprocessing import Pool
import numpy as np

#########################
#
# Use this script to read in a coordinate sorted bam file (with its bai index present)
# and output a CSV file identifying the number of CpGs present
# and how many reads cover ALL present CpGs.
# Edit the input params to output the correct chromosome and bin size desired.
#
# CalculateCompleteBins.py input_file.bam
#
##########################

class CalculateComplteBins:

    def __init__(self, input_bam_file, bin_size: int, number_of_processors=1):
        """
        This class is initialized with a path to a bam file and a bin size
        :param input_bam_file: One of the BAM files for analysis to be performed
        :param bin_size: Size of the bins for the analysis, integer
        :number_of_processors: How many CPUs to use for parallel computation, default=1
        """
        self.input_bam_file = input_bam_file
        self.bin_size = bin_size
        self.number_of_processors = number_of_processors

    def CalculateBinCoverate(self, bin):
        """
        Take a single bin, return a matrix
        :param bin: Bin should be passed as "Chr19_4343343"
        :return: pd.DataFrame with rows containing NaNs dropped
        """

        # Get reads from bam file
        parser = BamFileReadParser(self.input_bam_file, 20)
        # Split bin into parts
        chromosome, bin_location = bin.split("_")
        try:
            reads = parser.parse_reads(chromosome, int(bin_location)-self.bin_size, int(bin_location))
        except ValueError:
            # No reads are within this window, do nothing
            logging.info("No reads found in the window of {}:{}-{}".format(chromosome, str(int(bin_location-100)), bin_location))
            return None

        # convert to data_frame of 1s and 0s, drop rows with NaN
        matrix = parser.create_matrix(reads)
        matrix = matrix.dropna()

        return matrix

    def GetChromosomeLenghts(self):
        """
        Get dictionary containing lengths of the chromosomes. Uses bam file for reference
        :return: Dictionary of chromosome lengths, ex: {"chrX": 222222}
        """
        parser = BamFileReadParser(self.input_bam_file, 20)
        return dict(zip(parser.OpenBamFile.references, parser.OpenBamFile.lengths))




def CalculateBinCoverage(bin):
    parser = BamFileReadParser(input_bam_file, 20)
    try:
        reads = parser.parse_reads(chromosome, start_pos, stop_pos)
        matrix = parser.create_matrix(reads)
    except ValueError:
        # No reads in this window
        logging.info("No reads found in the window of {} to {}".format(start_pos, stop_pos))

    matrix = matrix.dropna()

    return matrix




if __name__ == "__main__":

    # Input params
    # todo add these as command line args using argparse
    input_bam_file = sys.argv[1]
    bin_size = 100
    chromosome = 'chr19'
    log_file = "CalculateCompleteBins.{}.log".format(os.path.basename(input_bam_file))
    output_filename = "CalculateCompleteBins.{}.csv".format(os.path.basename(input_bam_file))
    BASE_DIR = os.path.dirname(input_bam_file)

    logging.basicConfig(filename=os.path.join(BASE_DIR, log_file), level=logging.DEBUG)


    logging.info("Input file is {}".format(input_bam_file))
    logging.info("Input params, bin size: {}".format(bin_size))
    logging.info("Chromosome: {}".format(chromosome))

    # Setup parser object on the bam file
    chrom_lengths = dict(zip(parser.OpenBamFile.references, parser.OpenBamFile.lengths))

    # Open output file for writing
    output_file = open(os.path.join(BASE_DIR, output_filename), 'w')
    #output_file.write("chromosome,start,stop,full_reads,CpGs\n")


    # generate all bins for chromosome
    bins = list(np.arange(bin_size, chrom_lengths[chromosome] + bin_size, bin_size))
    # convert to input format
    bins = ["{}_".format(chromosome) + str(x) for x in bins]

    # Start looping over the bam file
    pool = Pool(processes=4)
    results = pool.map(CalculateBinCoverage, bins)

    for result in results:
        # todo bin must be linked with matrix output maybe use temp files and merge
        output_file.write(chromosome + ",")
        output_file.write(str(start_pos) + ",")
        output_file.write(str(stop_pos) + ",")
        output_file.write(str(result.shape[0]) + ",")
        output_file.write(str(result.shape[1]) + "\n")


    output_file.close()
    logging.info("Complete")

