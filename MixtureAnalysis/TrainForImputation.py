import pandas as pd
import numpy as np
from tqdm import tqdm
import time
from multiprocessing import Pool
from MixtureAnalysis.ConnectToCpGNet import TrainWithCpGNet
from MixtureAnalysis.ParseBam import BamFileReadParser


class TrainForImputation:

    def __init__(self, cpg_density: int, bam_file: str, mbias_read1_5=None, 
        mbias_read1_3=None, mbias_read2_5= None, mbias_read2_3=None, processes=-1):
        self.cpg_density = cpg_density
        self.bam_file = bam_file
        self.mbias_read1_5 = mbias_read1_5
        self.mbias_read1_3 = mbias_read1_3
        self.mbias_read2_5 = mbias_read2_5
        self.mbias_read2_3 = mbias_read2_3
        self.processes = processes



    def extract_matrices(self, coverage_data_frame: pd.DataFrame):

        def track_progress(job, update_interval=30):
            while job._number_left > 0:
                print("Tasks remaining = {0}".format(job._number_left * job._chunksize), flush=True)
                time.sleep(update_interval)

        subset = coverage_data_frame[coverage_data_frame['cpgs'] == self.cpg_density]
        bins_of_interest = subset['bin'].unique()

        pool = Pool(processes=self.processes)
        matrices = pool.map_async(self.multiprocess_extract, bins_of_interest)

        track_progress(matrices)
        matrices = matrices.get()

        clean_matrices = []
        for matrix in matrices:
            if matrix.shape[1] == self.cpg_density:
                clean_matrices.append(matrix)

        return np.array(clean_matrices)


    def multiprocess_extract(self, one_bin: str):
        read_parser = BamFileReadParser(self.bam_file, 20, read1_5=self.mbias_read1_5, read1_3=self.mbias_read1_3, read2_5=self.mbias_read2_5, read2_3=self.mbias_read2_3)
        chrom, loc = one_bin.split("_")
        loc = int(loc)
        reads = read_parser.parse_reads(chrom, loc-100, loc)
        matrix = read_parser.create_matrix(reads)
        matrix = matrix.dropna(how="all")
        matrix = matrix.fillna(-1)
        matrix = np.array(matrix)

        return matrix


    def train_model(self, output_folder: str, matrices: iter):
        train_net = TrainWithCpGNet(cpg_density=self.cpg_density, save_path=output_folder)
        model = train_net.train_model(matrices)

        return model
