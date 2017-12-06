from ParseBam import BamFileReadParser
import sys
import os
import logging
from multiprocessing import Pool
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score


# Input params
input_bam_file = sys.argv[1]
bin_size = 100
chromosome = 'chr19'
log_file = "CalculateCompleteBinsMulti.{}.log".format(os.path.basename(input_bam_file))
output_filename = "CalculateCompleteBinsMulti.{}.csv".format(os.path.basename(input_bam_file))
BASE_DIR = os.path.dirname(input_bam_file)

logging.basicConfig(filename=os.path.join(BASE_DIR, log_file), level=logging.DEBUG)


logging.info("Input file is {}".format(input_bam_file))
logging.info("Input params, bin size: {}".format(bin_size))
logging.info("Chromosome: {}".format(chromosome))

# Setup parser object on the bam file
# parser = BamFileReadParser(input_bam_file, 20)
# chrom_lengths = dict(zip(parser.OpenBamFile.references, parser.OpenBamFile.lengths))


def get_matrix_size(bin_loc):
    parser_b = BamFileReadParser(input_bam_file, 20)
    start_pos = bin_loc - 100
    stop_pos = bin_loc
    try:
        reads = parser_b.parse_reads(chromosome, start_pos, stop_pos)
        matrix = parser_b.create_matrix(reads)
    except ValueError:
        # No reads present
        return None

    matrix = matrix.dropna()

    return matrix.shape

pool = Pool(processes=2)

bins = pd.read_table("chr19_bins.txt")
bins = np.array(bins)

print("Getting results...")

results = pool.map(get_matrix_size, bins)

print("Got results")

output = open("test.txt", 'w')

for result in results:
    output.write(str(result))
    output.write("\n")

output.close()