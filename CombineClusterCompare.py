from sklearn.cluster import DBSCAN
from collections import Counter
import pandas as pd
import numpy as np
import sys
import logging
import os
from ParseBam import BamFileReadParser
import argparse


def filter_data_frame(matrix: pd.DataFrame, cluster_memeber_min):
    output = matrix.copy()
    for cluster in output['class'].unique():
        df = output[output['class'] == cluster]
        if len(df) < cluster_memeber_min:
            indexes = df.index
            output.drop(indexes, inplace=True)

    return output



if __name__ == "__main__":

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-a", "--input_bam_A",
                            help="First Input bam file, coordinate sorted with index present")
    arg_parser.add_argument("-b", "--input_bam_B",
                            help="Second Input bam file, coordinate sorted with index present")
    arg_parser.add_argument("-o", "--output_dir",
                            help="Output directory to save figures, defaults to bam file loaction")
    arg_parser.add_argument("-bins",
                            help="File with each line being one bin to extract and analyze, "
                                 "generated by CalculateCompleteBins")
    arg_parser.add_argument("-bin_size", help="Size of bins to extract and analyze, default=100", default=100)
    arg_parser.add_argument("-m", "--cluster_member_minimum",
                            help="Minimum number of members a cluster should have for it to be considered, default=4",
                            default=4)
    arg_parser.add_argument("-r", "--read_depth",
                            help="Minium number of reads covering all CpGs that the bins should have to analyze, default=20",
                            default=20)

    args = arg_parser.parse_args()

    input_bam_a = args.input_bam_A
    input_bam_b = args.input_bam_B
    bins_file = args.bins
    bin_size = int(args.bin_size)
    cluster_min = int(args.cluster_member_minimum)
    read_depth_req = int(args.read_depth)


    # Check all inputs are supplied
    if not input_bam_a or not input_bam_b or not bins_file:
        raise FileNotFoundError("Please make sure all required input files are supplied")


    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = os.path.dirname(input_bam_a)

    # Create output dir if it doesnt exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    bam_parser_A = BamFileReadParser(input_bam_a, 20)
    bam_parser_B = BamFileReadParser(input_bam_b, 20)

    # Read in bins
    bins=[]
    with open(bins_file, 'r') as f:
        for line in f:
            data = line.strip().split(",")
            bins.append("_".join([data[0], data[2]]))

    # a list to hold the clusters we identify as interesting
    bins_with_unique_clusters = []

    # Open output file for writing all data and unique data
    output_all = open(os.path.join(output_dir, "CombinedClusterCompare_output_all.csv"), 'w')
    output_unique = open(os.path.join(output_dir, "CombinedClusterCompare_output_unique.csv"), 'w')

    # Go through and cluster reads from each bin
    for bin in bins:
        chromosome, bin_loc = bin.split("_")
        bin_loc = int(bin_loc)

        reads_A = bam_parser_A.parse_reads(chromosome, bin_loc-bin_size, bin_loc)
        reads_B = bam_parser_B.parse_reads(chromosome, bin_loc-bin_size, bin_loc)
        matrix_A = bam_parser_A.create_matrix(reads_A)
        matrix_B = bam_parser_B.create_matrix(reads_B)

        # drop reads without full coverage of CpGs
        matrix_A = matrix_A.dropna()
        matrix_B = matrix_B.dropna()

        # if read depths are still not a minimum, skip
        if matrix_A.shape[0] < read_depth_req or matrix_B.shape[0] < read_depth_req:
            print("{}: Failed read req with out {} reads in one file".format(bin, str(matrix_B.shape[0])))
            continue

        # create labels and add to dataframe
        labels_A = ['A'] * len(matrix_A)
        labels_B = ['B'] * len(matrix_B)
        matrix_A['input'] = labels_A
        matrix_B['input'] = labels_B

        full_matrix = pd.concat([matrix_A, matrix_B])
        data_to_cluster = np.matrix(full_matrix)[:, :-1]

        # Create DBSCAN classifier and cluster add cluster classes to df
        clf = DBSCAN(min_samples=2)
        labels = clf.fit_predict(data_to_cluster)
        full_matrix['class'] = labels

        # Filter out any clusters with less than a minimum
        full_matrix = filter_data_frame(full_matrix, cluster_min)
        total_clusters = len(full_matrix['class'].unique())  # for output

        # Calculate clusters for A and B
        A_clusters = len(full_matrix[full_matrix['input'] == 'A']['class'].unique())  # for output
        B_clusters = len(full_matrix[full_matrix['input'] == 'B']['class'].unique())  # for output

        # todo fix a bug here where unique always equals total
        # Calculate how many clusters are unique to A or B
        num_unique_classes = 0 # for output
        # print(full_matrix.sort_values('class'))
        for label in full_matrix['class'].unique():
            df = full_matrix[full_matrix['class'] == label]
            # This cluster is unique to only one input
            if len(df['input'].unique()) == 1:
                num_unique_classes += 1

        # Write this data for an output
        output_line = ",".join([bin, str(total_clusters), str(A_clusters), str(B_clusters), str(num_unique_classes)])
        output_all.write(output_line + "\n")
        output_all.flush()
        if num_unique_classes > 0:
            output_unique.write(output_line + "\n")
            output_unique.flush()

    output_all.close()
    output_unique.close()
    print("Done")