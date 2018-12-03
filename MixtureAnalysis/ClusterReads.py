from sklearn.cluster import DBSCAN
import pandas as pd
from pandas.core.indexes.base import InvalidIndexError
import numpy as np
import sys
import logging
import os
from MixtureAnalysis.ParseBam import BamFileReadParser
from MixtureAnalysis.OutputComparisonResults import OutputIndividualMatrixData
from MixtureAnalysis.Imputation import Imputation
import keras.backend as K
import argparse
import datetime
from multiprocessing import Pool
from collections import defaultdict
import time
import tempfile
import pickle

class ClusterReads:

    def __init__(self, bam_a: str, bam_b=None, bin_size=100, bins_file=None, output_directory=None, num_processors=1,
        cluster_member_min=4, read_depth_req=10, remove_noise=True, mbias_read1_5=None, 
        mbias_read1_3=None, mbias_read2_5=None, mbias_read2_3=None, suffix="", no_overlap=True):

        self.bam_a = bam_a
        self.bam_b = bam_b
        self.bin_size = int(bin_size)
        self.bins_file = bins_file
        self.output_directory = output_directory
        self.num_processors = num_processors
        self.cluster_member_min = cluster_member_min
        self.read_depth_req = read_depth_req
        self.remove_noise = remove_noise
        self.mbias_read1_5 = mbias_read1_5
        self.mbias_read1_3 = mbias_read1_3
        self.mbias_read2_5 = mbias_read2_5
        self.mbias_read2_3 = mbias_read2_3
        self.suffix = suffix
        self.no_overlap = no_overlap
        
        if bam_b:
            self.single_file_mode = False
        else:
            self.single_file_mode = True

    # Remove clusters with less than n members
    def filter_data_frame(self, matrix: pd.DataFrame):
        output = matrix.copy()
        # duplicate indexes exist from concatenation, reset in the index to prevent dropping unintended rows
        output.reset_index(drop=True, inplace=True)
        for cluster in output['class'].unique():
            df = output[output['class'] == cluster]
            if len(df) < self.cluster_member_min:
                indexes = df.index
                output.drop(indexes, inplace=True)

        return output

    # Get only the matrices made up of reads from A OR B
    @staticmethod
    def get_unique_matrices(filtered_matrix):
        unique_dfs = []
        for label in filtered_matrix['class'].unique():
            df = filtered_matrix[filtered_matrix['class'] == label]
            if len(df['input'].unique()) == 1:
                unique_dfs.append(df)

        return unique_dfs

    # Get matrices with reads made up of A AND B
    @staticmethod
    def get_common_matrices(filtered_matrix):
        shared_dfs = []
        for label in filtered_matrix['class'].unique():
            df = filtered_matrix[filtered_matrix['class'] == label]
            if len(df['input'].unique()) > 1:
                shared_dfs.append(df)

        return shared_dfs

    # Get the means for all unique matrices
    def get_unique_means(self, filtered_matrix):
        output = []
        for matrix in self.get_unique_matrices(filtered_matrix):
            input_file = matrix['input'].unique()[0]
            matrix_mean = np.array(matrix)[:, :-2].mean()
            output.append((input_file, matrix_mean))

        return output

    # Get the means for all common matrices
    def get_common_means(self, filtered_matrix):
        output = []
        for matrix in self.get_common_matrices(filtered_matrix):
            matrix_mean = np.array(matrix)[:, :-2].mean()
            output.append(matrix_mean)

        return output

    # Generate a string label for each bin
    @staticmethod
    def make_bin_label(chromosome, stop_loc):
        return "_".join([chromosome, str(stop_loc)])

    # Takes the output of process_bins() and converts it into list of lines of data for output
    def generate_individual_matrix_data(self, filtered_matrix, chromosome, bin_loc):
        # Individual comparisons data
        lines = []
        unique_groups = self.get_unique_matrices(filtered_matrix)
        common_groups = self.get_common_matrices(filtered_matrix)
        bin_label = self.make_bin_label(chromosome, bin_loc)

        for matrix in unique_groups:
            cpg_matrix = np.array(matrix.drop(['class', 'input'], axis=1))
            # get a semi-colon separated string of 1s and 0s representing the CpG pattern
            cpg_pattern = ";".join([str(int(x)) for x in list(cpg_matrix[0])])
            m_mean = cpg_matrix.mean()
            num_cpgs = cpg_matrix.shape[1]
            read_number = len(matrix)
            input_label = matrix['input'].unique()[0]
            class_label = matrix['class'].unique()[0]
            out_line = ",".join([bin_label, input_label, str(m_mean), str(class_label), str(read_number),
                                str(num_cpgs), cpg_pattern])
            lines.append(out_line)

        for matrix in common_groups:
            cpg_matrix = np.array(matrix.drop(['class', 'input'], axis=1))
            # get a semi-colon separated string of 1s and 0s representing the CpG pattern
            cpg_pattern = ";".join([str(int(x)) for x in list(cpg_matrix[0])])
            m_mean = cpg_matrix.mean()
            num_cpgs = cpg_matrix.shape[1]
            read_number = len(matrix)
            input_label = "".join(list(matrix['input'].unique()))
            class_label = matrix['class'].unique()[0]
            out_line = ",".join([bin_label, input_label, str(m_mean), str(class_label), str(read_number),
                                str(num_cpgs), cpg_pattern])
            lines.append(out_line)

        return lines


    # MAIN METHOD
    def process_bins(self, bin):
        """
        This is the main method and should be called using Pool.map It takes one bin location and uses the other helper
        functions to get the reads, form the matrix, cluster it with DBSCAN, and output the cluster data as text lines
        ready to writing to a file.
        :param bin: string in this format: "chr19_55555"
        :return: a list of lines representing the cluster data from that bin
        """
        chromosome, bin_loc = bin.split("_")
        bin_loc = int(bin_loc)

        # Create bam parser and parse reads
        bam_parser_A = BamFileReadParser(self.bam_a, 20, read1_5=self.mbias_read1_5, read1_3=self.mbias_read1_3,
                                        read2_5=self.mbias_read2_5, read2_3=self.mbias_read2_3, no_overlap=self.no_overlap)
        reads_A = bam_parser_A.parse_reads(chromosome, bin_loc - self.bin_size, bin_loc)

        if not self.single_file_mode:
            bam_parser_B = BamFileReadParser(self.bam_b, 20, read1_5=self.mbias_read1_5, read1_3=self.mbias_read1_3,
                                            read2_5=self.mbias_read2_5, read2_3=self.mbias_read2_3, no_overlap=self.no_overlap)
            reads_B = bam_parser_B.parse_reads(chromosome, bin_loc - self.bin_size, bin_loc)

        # This try/catch block returns None for a bin if any discrepancies in the data format of the bins are detected.
        # The Nones are filtered out during the output of the data
        try:
            #create matrix  drop NA
            # This matrix is actually a pandas dataframe
            matrix_A = bam_parser_A.create_matrix(reads_A).dropna()
            if not self.single_file_mode:
                matrix_B = bam_parser_B.create_matrix(reads_B).dropna()

        except ValueError as e:
            logging.error("ValueError when creating matrix at bin {}. Stack trace will be below if log level=DEBUG".format(bin))
            logging.debug(str(e))
            return None
        except InvalidIndexError as e:
            logging.error("Invalid Index error when creating matrices at bin {}".format(bin))
            logging.debug(str(e))
            return None

        # if read depths are still not a minimum, skip
        if matrix_A.shape[0] < self.read_depth_req:
            return None
        if not self.single_file_mode:
            if matrix_B.shape[0] < self.read_depth_req:
                return None

        # create labels and add to dataframe

        # If two files label each A and B, otherwise use file_name as label
        if not self.single_file_mode:
            labels_A = ['A'] * len(matrix_A)
            matrix_A['input'] = labels_A
            labels_B = ['B'] * len(matrix_B)
            matrix_B['input'] = labels_B
        else:
            labels_A = [os.path.basename(self.bam_a)] * len(matrix_A)
            matrix_A['input'] = labels_A

        if not self.single_file_mode:
            try:
                full_matrix = pd.concat([matrix_A, matrix_B])
            except ValueError as e:
                logging.error("Matrix concat error in bin {}".format(bin))
                # logging.debug(str(e))
                return None
        else:
            full_matrix = matrix_A

        # Get data without labels for clustering
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
        filtered_matrix = self.filter_data_frame(full_matrix)
        if self.remove_noise:
            filtered_matrix = filtered_matrix[filtered_matrix['class'] != -1]

        # return generate_output_data(filtered_matrix, chromosome, bin_loc)
        return self.generate_individual_matrix_data(filtered_matrix, chromosome, bin_loc)

    def execute(self):
        start_time = datetime.datetime.now().strftime("%y-%m-%d")
        def track_progress(job, update_interval=60):
            while job._number_left > 0:
                logging.info("Tasks remaining = {0}".format(
                    job._number_left * job._chunksize
                ))
                time.sleep(update_interval)

        bins = []
        with open(self.bins_file, "r") as f:
            for line in f:
                data = line.strip().split(",")
                bins.append(data[0])

        pool = Pool(processes=self.num_processors)
        results = pool.map_async(self.process_bins, bins)
        track_progress(results)
        results = results.get()

        output = OutputIndividualMatrixData(results)
        output.write_to_output(self.output_directory, "Clustering.{}{}.{}".format(os.path.basename(self.bam_a), self.suffix, start_time))


class ClusterReadsWithImputation(ClusterReads):
    
    def __init__(self, bam_a: str, bam_b=None, bin_size=100, bins_file=None, output_directory=None, num_processors=1,
        cluster_member_min=4, read_depth_req=10, remove_noise=True, mbias_read1_5=None, 
        mbias_read1_3=None, mbias_read2_5=None, mbias_read2_3=None, suffix="", no_overlap=True, models_A=None, models_B=None):
        self.models_A = models_A
        self.models_B = models_B

        super().__init__(bam_a, bam_b, bin_size, bins_file, output_directory, 
        num_processors, cluster_member_min, read_depth_req, remove_noise, 
        mbias_read1_5, mbias_read1_3, mbias_read2_5, mbias_read2_3, suffix, no_overlap)



    def get_coverage_data(self, cpg_density=None):
        coverage_data = pd.read_csv(self.bins_file, header=None)
        coverage_data.columns = ['bin', 'reads', 'cpgs']

        return coverage_data

    # def create_dictionary(self, bins_A, matrices_A, bins_B=None, matrices_B=None):
    #     output = defaultdict(list)
    #     for b, m in zip(bins_A, matrices_A):
    #         output[b].append(m)
    #     if bins_B is not None and matrices_B is not None:
    #         for b, m in zip(bins_B, matrices_B):
    #             output[b].append(m)
    #     return output

    def create_dictionary(self, bins, matrices):
        output = dict()
        for b, m in zip(bins, matrices):
            output[b] = m

        return output


    @staticmethod
    def filter_coverage_data(coverage_data, cpg_density):
        return coverage_data[coverage_data['cpgs'] == cpg_density]

    
    def execute(self):
        start_time = datetime.datetime.now().strftime("%y-%m-%d")
        def track_progress(job, update_interval=60):
            while job._number_left > 0:
                logging.info("Tasks remaining = {0}".format(
                    job._number_left * job._chunksize
                ))
                time.sleep(update_interval)
        
        coverage_data = self.get_coverage_data()

        # start the main loop for imputation of these CpGs
        final_results_tf = tempfile.TemporaryFile(mode="w+t")
        for i in range(2,7):
            print("Starting Imputation of CpG density {}...".format(i))
            imputer_A = Imputation(cpg_density=i, 
            bam_file=self.bam_a, 
            mbias_read1_5=self.mbias_read1_5,
            mbias_read1_3=self.mbias_read1_3, 
            mbias_read2_5=self.mbias_read2_5, 
            mbias_read2_3=self.mbias_read2_3,
            processes=self.num_processors,
            )
            if self.bam_b:
                imputer_B = Imputation(cpg_density=i, 
                bam_file=self.bam_b, 
                mbias_read1_5=self.mbias_read1_5,
                mbias_read1_3=self.mbias_read1_3, 
                mbias_read2_5=self.mbias_read2_5, 
                mbias_read2_3=self.mbias_read2_3,
                processes=self.num_processors
                )

            # Subset for CpG density
            sub_coverage_data = self.filter_coverage_data(coverage_data, i)

            # Split into chunks for memory management
            n = 10000
            chunks = [sub_coverage_data[i * n:(i + 1) * n] for i in range((len(sub_coverage_data) + n - 1) // n )]

            for j, chunk in enumerate(chunks):
                print("Divided into {} chunks for processing...".format(len(chunks)), flush=True)
                print("Extracting matrices from chunk {}...".format(j+1))
                bins_A, matrices_A = imputer_A.extract_matrices(chunk, return_bins=True)

                if self.bam_b:
                    print("Extracting from second imput bam...", flush=True)
                    bins_B, matrices_B = imputer_B.extract_matrices(chunk, return_bins=True)
                else:
                    bins_B = None
                    matrices_B = None

                # Create dictionary with key as bin and values as a list like [matrix_A, matrix_B]
                data_A_dict = self.create_dictionary(bins_A, matrices_A)
                if self.bam_b:
                    data_B_dict = self.create_dictionary(bins_B, matrices_B)

                # ### TESTING ###
                # with open("unimputed_matrix_dict_A_cpgs{}_chunk{}.pickle".format(i,j), "w+b") as f:
                #     pickle.dump(data_A_dict, f)
                # if self.bam_b:
                #     with open("unimputed_matrix_dict_B_cpgs{}_chunk{}.pickle".format(i,j), "w+b") as f:
                #         pickle.dump(data_B_dict, f)
                # ### TESTING ###

                # Attempt to impute
                print("Imputing chunk {}...".format(j), flush=True)
                imputed_matrices_A = imputer_A.impute_from_model(self.models_A, list(data_A_dict.values()))
                if self.bam_b:
                    imputed_matrices_B = imputer_B.impute_from_model(self.models_B, data_B_dict.values())
                else:
                    imputed_matrices_B = None

                data_imputed_A_dict = self.create_dictionary(data_A_dict.keys(), imputed_matrices_A)
                if self.bam_b:
                    data_imputed_B_dict = self.create_dictionary(data_B_dict.keys(), imputed_matrices_B)

                # ### TESTING ###
                # with open("IMPUTED_matrix_dict_A_cpgs{}_chunk{}.pickle".format(i,j), "w+b") as f:
                #     pickle.dump(data_imputed_A_dict, f)
                # if self.bam_b:
                #     with open("IMPUTED_matrix_dict_B_cpgs{}_chunk{}.pickle".format(i,j), "w+b") as f:
                #         pickle.dump(data_imputed_B_dict, f)
                # ### TESTING ###

                # Combine and cluster as normal
                # Maybe do this in the future
                def _multiprocess_cluster_work(bin_):
                    pass

                for bin_ in data_imputed_A_dict.keys():
                    matrix_A = data_imputed_A_dict[bin_]
                    matrix_A = pd.DataFrame(matrix_A)
                    if matrix_A.shape[0] < self.read_depth_req:
                        # TODO OUTPUT SOME EMTPY RESULT
                        continue
                    if self.bam_b:
                        matrix_B = data_imputed_B_dict[bin_]
                        matrix_B = pd.DataFrame(matrix_B)
                        if matrix_B.shape[0] < self.read_depth_req:
                            # TODO SAME AS ABOVE
                            continue
                    
                    if self.bam_b:
                        labels_A = ['A'] * len(matrix_A)
                        matrix_A['input'] = labels_A
                        labels_B = ['B'] * len(matrix_B)
                        matrix_B['input'] = labels_B
                    else:
                        labels_A = [os.path.basename(self.bam_a)] * len(matrix_A)
                        matrix_A['input'] = labels_A

                    if self.bam_b:
                        try:
                            full_matrix = pd.concat([matrix_A, matrix_B])
                        except ValueError as e:
                            logging.error("Matrix concat eror in bin {}".format(bin_))
                            continue
                    else:
                        full_matrix = matrix_A

                    # get data to cluster
                    data_to_cluster = np.matrix(full_matrix)[:,:-1]
                    clf = DBSCAN(min_samples=2)
                    try:
                        labels = clf.fit_predict(data_to_cluster)
                    except ValueError as e:
                        logging.error("ValueError when trying to cluster bin {}".format(bin_))
                        continue

                    full_matrix['class'] = labels

                    filtered_matrix = super().filter_data_frame(full_matrix)
                    if self.remove_noise:
                        filtered_matrix = filtered_matrix[filtered_matrix['class'] != -1]
                    
                    # generate output lines
                    chromosome, bin_loc = bin_.split("_")
                    output_lines = super().generate_individual_matrix_data(filtered_matrix, chromosome, bin_loc)
                    for line in output_lines:
                        final_results_tf.write(line+"\n")

        # TODO CLUSTER ALL OTHER BINS LIKE NORMAL WITHOUT IMPUTATION
        
        final_results_tf.seek(0)
        with open("final_results.csv", "w") as final:
            for line in final_results_tf:
                final.write(line)
                