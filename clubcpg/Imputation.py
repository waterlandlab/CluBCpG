import pandas as pd
import numpy as np
import logging
import os
from clubcpg.ConnectToCpGNet import TrainWithPReLIM
from clubcpg.ParseBam import BamFileReadParser
from clubcpg_prelim import PReLIM
from pebble import ProcessPool
from joblib import load


class Imputation:
    """The class providing convienent APIs to train models and impute from models using PReLIM
    """

    def __init__(self, cpg_density: int, bam_file: str, mbias_read1_5=None, 
        mbias_read1_3=None, mbias_read2_5= None, mbias_read2_3=None, processes=-1):
        """[summary]
        
        Arguments:
            cpg_density {int} -- Number of CpGs this class instance will be used for
            bam_file {str} -- path to the bam file
        
        Keyword Arguments:
            mbias_read1_5 {[type]} -- [description] (default: {None})
            mbias_read1_3 {[type]} -- [description] (default: {None})
            mbias_read2_5 {[type]} -- [description] (default: {None})
            mbias_read2_3 {[type]} -- [description] (default: {None})
            processes {int} -- number or CPUs to use when parallelization can be utilized, default= All available (default: {-1})
        """

        self.cpg_density = cpg_density
        self.bam_file = bam_file
        self.mbias_read1_5 = mbias_read1_5
        self.mbias_read1_3 = mbias_read1_3
        self.mbias_read2_5 = mbias_read2_5
        self.mbias_read2_3 = mbias_read2_3
        self.processes = processes

    def extract_matrices(self, coverage_data_frame: pd.DataFrame, sample_limit: int = None, return_bins=False):
        """Extract CpG matrices from bam file.
        
        Arguments:
            coverage_data_frame {pd.DataFrame} -- Output of clubcpg-coverage read in as a csv file
        
        Keyword Arguments:
            return_bins {bool} -- Return the bin location along with the matrix (default: {False})
        
        Returns:
            [tuple] -- Returns tuple of (bin, np.array) if returns_bins = True else returns only np.array
        """

        subset = coverage_data_frame[coverage_data_frame['cpgs'] == self.cpg_density]
        bins_of_interest = subset['bin'].unique()

        # Downsample the training bins if requested and necessary
        if sample_limit and len(bins_of_interest) > sample_limit:
            bins_of_interest = np.random.choice(bins_of_interest, size=sample_limit)

        # Use the pebbel ProcessPool because it can handle hanging processes with a timeout
        complete_results = []
        with ProcessPool(max_workers=self.processes) as pool:
            future = pool.map(self._multiprocess_extract, bins_of_interest, timeout=5)

            iterator = future.result()

            while True:
                try:
                    result = next(iterator)
                    complete_results.append(result)
                except StopIteration:
                    break
                except TimeoutError as error:
                    print("Timeout caught - {}".format(error.args[1]))
                except Exception as error:
                    print("Unknown exception = {}".format(error))

        bins, matrices = zip(*complete_results)

        # destroy the pool
        pool.close()

        # Remove any potential bad data
        clean_matrices = []
        clean_bins = []
        for matrix, bin_ in zip(matrices, bins):
            try:
                if matrix.shape[1] == self.cpg_density:
                    clean_matrices.append(matrix)
                    clean_bins.append(bin_)
            except IndexError as e:
                logging.info("Index error at bin {}".format(bin_))
                logging.error(str(e))
                continue

        # if len(clean_matrices) > 0:
        #     clean_matrices = np.array(clean_matrices)

        clean_matrices = np.array(clean_matrices)

        if return_bins:
            return clean_bins, clean_matrices
        else:
            return clean_matrices


    def _multiprocess_extract(self, one_bin: str):
        """Function to be used for multiprocessing
        
        Arguments:
            one_bin {str} -- bin id as "chr7_222222"
        
        Returns:
            [tuple] -- bin, matrix
        """
        try:
            read_parser = BamFileReadParser(self.bam_file, 20, read1_5=self.mbias_read1_5, read1_3=self.mbias_read1_3, read2_5=self.mbias_read2_5, read2_3=self.mbias_read2_3)
            chrom, loc = one_bin.split("_")
            loc = int(loc)
            reads = read_parser.parse_reads(chrom, loc-100, loc) # TODO unhardcode bin size
            matrix = read_parser.create_matrix(reads)
            matrix = matrix.dropna(how="all")
            # if matrix.shape[0] == 0:
            #     return None
            matrix = matrix.fillna(-1)
            matrix = np.array(matrix)
            matrix = matrix.astype('int8')
        except: # BAD EXCEPTION
            return (one_bin, np.array([]))

        return (one_bin, matrix)


    def train_model(self, output_folder: str, matrices: iter):
        """Train a CpGNet model using :class:`.TrainWithCpGNet`
        
        Arguments:
            output_folder {str} -- Folder to save trained models
            matrices {iter} -- An iterable of CpGMatrices - ideally obtained through Imputation.extract_matrices()
        
        Returns:
            [keras model] -- Returns the trained CpGNet model
        """

        train_model = TrainWithPReLIM(cpg_density=self.cpg_density, save_path=output_folder)
        model = train_model.train_model(matrices)

        return model

    @staticmethod
    def postprocess_predictions(predicted_matrix):
        """Takes array with predicted values and rounds them to 0 or 1 if threshold is exceeded
        
        Arguments:
            predicted_matrix {[type]} -- matrix generated by imputation
        
        Returns:
            [type] -- predicted matrix predictions as 1, 0, or NaN
        """

        processed_array = []
        for array in predicted_matrix:
            new_array = []
            for item in array:
                if item != 1 and item != 0:
                    if item <= 0.2: #TODO un-hardcode this
                        new_array.append(0.0)
                    elif item >= 0.8: #TODO un-hardcode this
                        new_array.append(1.0)
                    else:
                        new_array.append(np.nan)
                else:
                    new_array.append(item)

            processed_array.append(new_array)

        return np.array(processed_array)




    def impute_from_model(self, models_folder: str, matrices: iter, postprocess=True):
        """Generator to provide imputed matrices on-the-fly
        
        Arguments:
            models_folder {str} -- Path to directory containing trained CpGNet models
            matrices {iter} -- An iterable containging n x m matrices with n=cpgs and m=reads
        
        Keyword Arguments:
            postprocess {bool} -- Round imputed values to 1s and 0s  (default: {True})
        """

        model_path = os.path.join(models_folder, "saved_model_{}_cpgs.prelim".format(self.cpg_density))

        trained_model = PReLIM(cpgDensity=self.cpg_density)
        print("Successfully loaded model: {}".format(model_path), flush=True)
        trained_model.model = load(model_path)

        for m in matrices:
            # only impute if there is an unknown
            if -1 in m:
                m = m.astype(float)
                pm = trained_model.impute(m)
                if postprocess:
                    pm = self.postprocess_predictions(pm)
            # Nothing to impute, passback original matrix to keep list in order
            else:
                pm = m.copy()
                
            # K.clear_session()
            yield pm   
