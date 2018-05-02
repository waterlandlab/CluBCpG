from MixtureAnalysis.ParseBam import BamFileReadParser
import sys
import os
import logging
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


# Adapted from Stackoverflow
# https://stackoverflow.com/questions/35611465/python-scikit-learn-clustering-with-missing-data


def kmeans_missing(X, n_clusters, max_iter=10):
    """Perform K-Means clustering on data with missing values.

    Args:
      X: An [n_samples, n_features] array of data to cluster.
      n_clusters: Number of clusters to form.
      max_iter: Maximum number of EM iterations to perform.

    Returns:
      labels: An [n_samples] vector of integer labels.
      centroids: An [n_clusters, n_features] array of cluster centroids.
      X_hat: Copy of X with the missing values filled in.
    """

    # Initialize missing values to their column means
    missing = ~np.isfinite(X)
    mu = np.nanmean(X, 0, keepdims=1)
    X_hat = np.where(missing, mu, X)

    for i in range(max_iter):
        if i > 0:
            # initialize KMeans with the previous set of centroids. this is much
            # faster and makes it easier to check convergence (since labels
            # won't be permuted on every iteration), but might be more prone to
            # getting stuck in local minima.
            cls = KMeans(n_clusters, init=prev_centroids)
        else:
            # do multiple random initializations in parallel
            cls = KMeans(n_clusters, n_jobs=-1)

        # perform clustering on the filled-in data
        labels = cls.fit_predict(X_hat)
        centroids = cls.cluster_centers_

        # fill in the missing values based on their cluster centroids
        X_hat[missing] = centroids[labels][missing]

        # when the labels have stopped changing then we have converged
        if i > 0 and np.all(labels == prev_labels):
            break

        prev_labels = labels
        prev_centroids = cls.cluster_centers_

    return labels, centroids, X_hat, missing


# Input params
input_bam_file = sys.argv[1]
bin_size = 100
chromosome = 'chr19'
matrix_read_req = 5
matrix_cpg_req = 3
log_file = "AnalyzeBamBins.{}.log".format(os.path.basename(input_bam_file))
output_filename = "AnalyzeBamBins.{}.csv".format(os.path.basename(input_bam_file))
BASE_DIR = os.path.dirname(input_bam_file)

logging.basicConfig(filename=os.path.join(BASE_DIR, log_file), level=logging.DEBUG)


logging.info("Input file is {}".format(input_bam_file))
logging.info("Input params, bin size: {}".format(bin_size))
logging.info("Chromosome: {}".format(chromosome))
logging.info("Matrix shape reqs is {}, {}".format(matrix_read_req, matrix_cpg_req))

# Setup parser object on the bam file
parser = BamFileReadParser(input_bam_file, 20)
chrom_lengths = dict(zip(parser.OpenBamFile.references, parser.OpenBamFile.lengths))

# Open output file for writing
output_file = open(os.path.join(BASE_DIR, output_filename), 'w')
output_file.write("chrom,start,stop,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,reads,cpgs\n")

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

    # Check the matrix shape for size requirements
    if matrix.shape[0] < matrix_read_req:
        logging.info("Less than the required number of reads in region {} to {}".format(start_pos, stop_pos))
        continue
    if matrix.shape[1] < matrix_read_req:
        logging.info("Less than the required number of CpGs in the region {} to {}".format(start_pos, stop_pos))
        continue

    # Passes checks, do the analysis for the bin
    # Determine number of clusters as explained above
    n_clusters = min(len(matrix), 24)
    labels, centroids, X_hat, missing = kmeans_missing(matrix, n_clusters)

    # Silouette analysis on the estimated values
    X = X_hat.round()
    range_n_clusters = np.arange(2, n_clusters, 1)

    sil_scores = []

    for n_clusters in range_n_clusters:
        # Cluster the data and get the cluster labels
        clusterer = KMeans(n_clusters=n_clusters, random_state=7)
        cluster_labels = clusterer.fit_predict(X)
        # If there is only one predicted label append a nan and continue to prevent crashing
        if len(set(cluster_labels)) < 2:
            logging.info("Only one predicted cluster for {} to {}".format(start_pos, stop_pos))
            sil_scores.append(np.nan)
            continue
        silhouette_avg = silhouette_score(X, cluster_labels)
        sil_scores.append(silhouette_avg)

    # Write data
    output_file.write("{},{},{},".format(chromosome, start_pos, stop_pos))
    output_file.write(",".join(str(s) for s in sil_scores) + ",")
    output_file.write(str(matrix.shape[0]) + ",")
    output_file.write(str(matrix.shape[1]) + "\n")


output_file.close()
logging.info("Complete")
