from sklearn.cluster import DBSCAN
from collections import Counter
import pandas as pd
import pickle
import numpy as np
import sys
import logging
import os
from ParseBam import BamFileReadParser
from OutputComparisonResults import OutputComparisonResults, OutputIndividualMatrixData
import argparse
import datetime
from multiprocessing import Pool

# Remove clusters with less than n members
def filter_data_frame(matrix: pd.DataFrame, cluster_memeber_min):
    output = matrix.copy()
    for cluster in output['class'].unique():
        df = output[output['class'] == cluster]
        if len(df) < cluster_memeber_min:
            indexes = df.index
            output.drop(indexes, inplace=True)

    return output


# Get only the matrices made up of reads from A or B
def get_unique_matrices(filtered_matrix):
    unique_dfs = []
    for label in filtered_matrix['class'].unique():
        df = filtered_matrix[filtered_matrix['class'] == label]
        if len(df['input'].unique()) == 1:
            unique_dfs.append(df)

    return unique_dfs


# Get matrices with reads made up of A and B
def get_common_matrices(filtered_matrix):
    shared_dfs = []
    for label in filtered_matrix['class'].unique():
        df = filtered_matrix[filtered_matrix['class'] == label]
        if len(df['input'].unique()) > 1:
            shared_dfs.append(df)

    return shared_dfs

# Get the means for all unique matrices
def get_unique_means(filtered_matrix):
    output = []
    for matrix in get_unique_matrices(filtered_matrix):
        input_file = matrix['input'].unique()[0]
        matrix_mean = np.array(matrix)[:, :-2].mean()
        output.append((input_file, matrix_mean))

    return output

# Get the means for all common matrices
def get_common_means(filtered_matrix):
    output = []
    for matrix in get_common_matrices(filtered_matrix):
        matrix_mean = np.array(matrix)[:, :-2].mean()
        output.append(matrix_mean)

    return output

# Generate a string label for each bin
def make_bin_label(chromosome, stop_loc):
    return "_".join([chromosome, str(stop_loc)])

#
def generate_output_data(filtered_matrix, chromosome, bin_loc):
    # Individual comparisons data
    lines = []
    unique_groups = get_unique_means(filtered_matrix)
    common_groups = get_common_means(filtered_matrix)
    bin_label = make_bin_label(chromosome, bin_loc)

    for group in unique_groups:
        file_input = group[0]
        mean = group[1]
        for common_group in common_groups:
            diff = mean - common_group
            line = ",".join([bin_label, file_input, str(mean), str(common_group), str(diff)])
            lines.append(line)

    # Bin summary data:
    num_unique = len(unique_groups)
    num_common = len(common_groups)
    num_total = num_unique + num_common
    summary_line = ",".join([bin_label, str(num_unique), str(num_common), str(num_total)])

    return summary_line, lines

# Takes the output of process_bins() and converts it into list of lines of data for output
def generate_individual_matrix_data(filtered_matrix, chromosome, bin_loc):
    # Individual comparisons data
    lines = []
    unique_groups = get_unique_matrices(filtered_matrix)
    common_groups = get_common_matrices(filtered_matrix)
    bin_label = make_bin_label(chromosome, bin_loc)

    for matrix in unique_groups:
        cpg_matrix = np.matrix(matrix.drop(['class', 'input'], axis=1))
        m_mean = cpg_matrix.mean()
        num_cpgs = cpg_matrix.shape[1]
        read_number = len(matrix)
        input_label = matrix['input'].unique()[0]
        class_label = matrix['class'].unique()[0]
        out_line = ",".join([bin_label, input_label, str(m_mean), str(class_label), str(read_number), str(num_cpgs)])
        lines.append(out_line)

    for matrix in common_groups:
        cpg_matrix = np.matrix(matrix.drop(['class', 'input'], axis=1))
        m_mean = cpg_matrix.mean()
        num_cpgs = cpg_matrix.shape[1]
        read_number = len(matrix)
        input_label = 'AB'
        class_label = matrix['class'].unique()[0]
        out_line = ",".join([bin_label, input_label, str(m_mean), str(class_label), str(read_number), str(num_cpgs)])
        lines.append(out_line)

    return lines


# Function to execute in parallel using Pool
# bin should be passed as "chr19_33444"
def process_bins(bin):

    bam_parser_A = BamFileReadParser(input_bam_a, 20, read1_5=mbias_read1_5, read1_3=mbias_read1_3,
                                     read2_5=mbias_read2_5, read2_3=mbias_read2_3)
    bam_parser_B = BamFileReadParser(input_bam_b, 20, read1_5=mbias_read1_5, read1_3=mbias_read1_3,
                                     read2_5=mbias_read2_5, read2_3=mbias_read2_3)

    chromosome, bin_loc = bin.split("_")
    bin_loc = int(bin_loc)

    reads_A = bam_parser_A.parse_reads(chromosome, bin_loc - bin_size, bin_loc)
    reads_B = bam_parser_B.parse_reads(chromosome, bin_loc - bin_size, bin_loc)
    try:
        matrix_A = bam_parser_A.create_matrix(reads_A)
        matrix_B = bam_parser_B.create_matrix(reads_B)
    except ValueError as e:
        logging.error("ValueError when creating matrix at bin {}. Stack trace will be below if log level=DEBUG".format(bin))
        logging.debug(str(e))
        return None

    # drop reads without full coverage of CpGs
    matrix_A = matrix_A.dropna()
    matrix_B = matrix_B.dropna()

    # if read depths are still not a minimum, skip
    if matrix_A.shape[0] < read_depth_req or matrix_B.shape[0] < read_depth_req:
        # print("{}: Failed read req with out {} reads in one file".format(bin, str(matrix_B.shape[0])))
        # sys.stdout.flush()
        return None

    # create labels and add to dataframe
    labels_A = ['A'] * len(matrix_A)
    labels_B = ['B'] * len(matrix_B)
    matrix_A['input'] = labels_A
    matrix_B['input'] = labels_B

    full_matrix = pd.concat([matrix_A, matrix_B])
    data_to_cluster = np.matrix(full_matrix)[:, :-1]

    # Create DBSCAN classifier and cluster add cluster classes to df
    clf = DBSCAN(min_samples=2)
    try:
        labels = clf.fit_predict(data_to_cluster)
    except ValueError as e:
        # log error
        logging.error("ValueError when trying to cluster bin {}".format(bin))
        logging.debug(str(e))
        return None

    full_matrix['class'] = labels

    # Filter out any clusters with less than a minimum
    filtered_matrix = filter_data_frame(full_matrix, cluster_min)

    # return generate_output_data(filtered_matrix, chromosome, bin_loc)
    return generate_individual_matrix_data(filtered_matrix, chromosome, bin_loc)


if __name__ == "__main__":

    # Set command line arguments
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
                            help="Minium number of reads covering all CpGs that the bins should have to analyze, default=10",
                            default=10)
    arg_parser.add_argument("-n", "--num_processors",
                            help="Number of processors to use for analysis, default=1",
                            default=1)
    arg_parser.add_argument("--read1_5", help="integer, read1 5' m-bias ignore bp, default=0", default=0)
    arg_parser.add_argument("--read1_3", help="integer, read1 3' m-bias ignore bp, default=0", default=0)
    arg_parser.add_argument("--read2_5", help="integer, read2 5' m-bias ignore bp, default=0", default=0)
    arg_parser.add_argument("--read2_3", help="integer, read2 3' m-bias ignore bp, default=0", default=0)

    args = arg_parser.parse_args()

    # Assign arg parser vars to new variables
    input_bam_a = args.input_bam_A
    input_bam_b = args.input_bam_B
    bins_file = args.bins
    bin_size = int(args.bin_size)
    cluster_min = int(args.cluster_member_minimum)
    read_depth_req = int(args.read_depth)
    num_processors = int(args.num_processors)

    # Get the mbias inputs and adjust to work correctly, 0s should be converted to None
    mbias_read1_5 = int(args.read1_5)
    mbias_read1_3 = int(args.read1_3)
    mbias_read2_5 = int(args.read2_5)
    mbias_read2_3 = int(args.read2_3)

    if mbias_read1_5 == 0:
        mbias_read1_5 = None
    if mbias_read1_3 == 0:
        mbias_read1_3 = None
    if mbias_read2_5 == 0:
        mbias_read2_5 = None
    if mbias_read2_3 == 0:
        mbias_read2_3 = None

    if mbias_read1_5 or mbias_read1_3 or mbias_read2_5 or mbias_read2_3:
        mbias_on = True
    else:
        mbias_on = False

    # todo write input params to log_file for record keeping

    # Check all inputs are supplied
    if not input_bam_a or not input_bam_b or not bins_file:
        print("Please make sure all required input files are supplied")
        logging.error("All required input files were not supplied. Exiting error code 1.")
        sys.exit(1)


    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = os.path.dirname(input_bam_a)

    # Create output dir if it doesnt exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Set up logging
    start_time = datetime.datetime.now().strftime("%y-%m-%d_%H-%M-%S")
    log_file = os.path.join(output_dir, "{}.log".format(start_time))
    logging.basicConfig(filename=os.path.join(output_dir, log_file), level=logging.DEBUG) #todo adjust this with a -v imput param

    # Read in bins
    bins=[]
    with open(bins_file, 'r') as f:
        for line in f:
            data = line.strip().split(",")
            # bins.append("_".join([data[0], data[2]]))
            bins.append(data[0])

    pool = Pool(processes=num_processors)
    logging.info("Starting workers pool, using {} processors".format(num_processors))
    if mbias_on:
        logging.info("M bias inputs received, ignoring the following:\nread 1 5': {}bp\n"
                     "read1 3': {}bp\nread2 5: {}bp\nread2 3': {}bp".format(mbias_read1_5, mbias_read1_3, mbias_read2_5, mbias_read2_3))

    # Results is a list of lists
    results = pool.map(process_bins, bins)

    # Convert the results into two output csv files for human analysis
    # output = OutputComparisonResults(results)
    output = OutputIndividualMatrixData(results)
    output.write_to_output(output_dir, start_time)

    logging.info("Done")
