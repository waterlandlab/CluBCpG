from ParseBam import BamFileReadParser
import sys
import os
import logging
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

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
parser = BamFileReadParser(input_bam_file, 20)
chrom_lengths = dict(zip(parser.OpenBamFile.references, parser.OpenBamFile.lengths))

# Open output file for writing
output_file = open(os.path.join(BASE_DIR, output_filename), 'w')
output_file.write("chromosome,start,stop,full_reads,CpGs\n")

# Start looping over the bam file
current_bin = 0
while current_bin <= chrom_lengths['chr19']:
    current_bin = current_bin + bin_size
    stop_pos = current_bin
    start_pos = stop_pos - bin_size
    try:
        reads = parser.parse_reads(chromosome, start_pos, stop_pos)
        matrix = parser.create_matrix(reads)
    except ValueError:
        # No reads in this window
        logging.info("No reads found in the window of {} to {}".format(start_pos, stop_pos))
        continue

    # Drop reads with empty values
    matrix = matrix.dropna()

    output_file.write(chromosome + ",")
    output_file.write(str(start_pos) + ",")
    output_file.write(str(stop_pos) + ",")
    output_file.write(str(matrix.shape[0]) + ",")
    output_file.write(str(matrix.shape[1]) + "\n")


output_file.close()
logging.info("Complete")
