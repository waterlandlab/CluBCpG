from CpGNet import CpGNet, TRAINED_CPG_MODELS
import numpy as np
from keras.models import load_model
import keras.backend as K

class ImputeWithCpGNet:

    def __init__(self, cpg_density: int, bin: str, bin_size=100):
        self.cpg_density = cpg_density
        self.bin = bin
        self.bin_size = bin_size
        self.model = self.get_model()

    def get_model(self):
        net = CpGNet(cpgDensity=self.cpg_density)
        net.model = load_model(TRAINED_CPG_MODELS[self.cpg_density])
        return net

    # Take the predicted matrix, convert to 1 and 0 when confident, otherwise nan
    def postprocess(self, predicted_matrix):
        processed_array = []
        for array in predicted_matrix:
            new_array = []
            for item in array:
                if item != 1 and item != 0:
                    if item <= 0.2:
                        new_array.append(0.0)
                    elif item >= 0.8:
                        new_array.append(1.0)
                    else:
                        new_array.append(np.nan)
                else:
                    new_array.append(item)

            processed_array.append(new_array)

        return np.array(processed_array)

    def impute_matrix(self, matrix, positions, postprocess=True):
        """
        :param matrix: 2D np array with unknowns as -1
        :param positions: 1D np array of chromosomal positions (ints)
        :param postprocess: Imply classes from imputation, if False, returns probabilities instead of predictions
        :return: A matrix imputed with highly confident predictions, less confident regions are assigned np.nan
        """
        bin_stop = int(self.bin.split("_")[1])
        bin_start = bin_stop - self.bin_size

        predicted_matrix = self.model.impute(matrix, positions, bin_start, bin_stop)
        if postprocess:
            predicted_matrix = self.postprocess(predicted_matrix)

        # This call forces tf to release memory, otherwise a large memory leak happens when executing in parallel
        K.clear_session()
        return predicted_matrix


class Train