from sklearn.cluster import DBSCAN
from collections import Counter
import pandas as pd
import numpy as np
import sys
import logging
import os
from ParseBam import BamFileReadParser

# Input params
# todo set this up as arg parser
input_bam_file = sys.argv[1]
bins_file = sys.argv[2]
bin_size = 100
chromosome = 'chr19'
log_file = "CalcDBSCAN.{}.log".format(os.path.basename(input_bam_file))
output_filename = "CalcDBSCAN.{}.csv".format(os.path.basename(input_bam_file))
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
output_file.write("chromosome,stop_loc,DBSCAN_clusters\n")

bins = []
with open(bins_file, 'r') as f:
    for line in f:
        data = line.strip().split(",")
        bins.append(data[2])

bins = list(map(int, bins))

for bin in bins:
    reads = parser.parse_reads(chromosome, bin-bin_size, bin)
    matrix = parser.create_matrix(reads)
    matrix = matrix.dropna()

    scn = DBSCAN(min_samples=2)
    cluster_labels = scn.fit_predict(matrix)

    # Count memebers of each cluster
    counter = Counter(cluster_labels)
    counts = list(counter.values())
    counts = [str(x) for x in counts]

    clusters = len(set(cluster_labels))

    output_file.write(str(chromosome) + ",")
    output_file.write(str(bin) + ",")
    output_file.write(str(clusters) + ",")
    output_file.write(":".join(counts))
    output_file.write("\n")
    output_file.flush()

output_file.close()

print("Done")
