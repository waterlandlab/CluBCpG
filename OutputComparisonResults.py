import os


# This class will be called by CombineClusterCompare.py to output its multiprocessing results output into csv files

class OutputIndividualMatrixData:

    def __init__(self, results):
        """
        :param results: Results objected generated by CompareClusterCompare.py's Pool.map of process_bins
        """
        self.results = results

    def write_to_output(self, filepath=None, prefix=None):
        """
        :param prefix: Prefix to name the output files
        :param filepath: Path to save the output files
        :return: None
        """

        if prefix and filepath:
            output_comparisons = open("{}_matrix_data.csv".format(os.path.join(filepath, prefix)), "w")

            output_comparisons.write("bin,input_label,methylation,class_label,read_number,cpg_number, cpg_pattern\n")

            for result in self.results:
                if result:
                    for line in result:
                        output_comparisons.write(line + "\n")
                else:
                    continue

            output_comparisons.close()

            return

