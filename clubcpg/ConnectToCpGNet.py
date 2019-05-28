from clubcpg_prelim import PReLIM
import os
from joblib import dump


class TrainWithPReLIM:
    """
    Used to train models using CpGnet
    """

    def __init__(self, cpg_density=None, save_path=None):
        """
        Class to train a CpGNet model from input data

        :param cpg_density: Number of CpGs
        :type cpg_density: int
        :param save_path: Location of folder to save the resulting model files. One per cpg density
        """
        if not cpg_density:
            raise AttributeError("CpG density must be specified")
        if not save_path:
            raise AttributeError("Folder to save trained model must be specified")
        self.save_path = save_path
        self.cpg_density = cpg_density
        self.model = PReLIM(cpgDensity=cpg_density)

    def save_net(self, model):
        """
        Save the network to a file

        :param model: The trained PReLIM model. Located at PReLIM.model
        :type model: :class:`clubcpg_prelim.PReLIM`
        :return: Path to the saved model
        """
        file_name = "saved_model_{}_cpgs.prelim".format(self.cpg_density)
        output = os.path.join(self.save_path, file_name)
        dump(model, output)
        print("Saved {} cpg model to {}".format(self.cpg_density, output))

        return output


    def train_model(self, bins: iter):
        """
        Train the CpGNet model on a list of provided bins

        :param bins: iterable containing CpG matrices of 1 (methylated), 0 (unmethylated), and -1 (unknown)
        :return: Path to the saved model file
        """
        self.model.train(bins, model_file="no")
        output = self.save_net(self.model.model)

        return output

