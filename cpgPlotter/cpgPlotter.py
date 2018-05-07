import matplotlib.pyplot as plt
import numpy as np


class CpGMatrixPlotter:

    def __init__(self):
        pass

    @staticmethod
    def buffer_spacings(spacings, radius):
        for i in range(len(spacings)):
            if i == 0:
                continue
            distance = spacings[i] - spacings[i - 1]
            if distance < radius:
                spacings[i] = spacings[i] + (radius - distance)

        return np.array(spacings)

    @staticmethod
    def get_range_from_bin(bin_str: str, bin_size=100):
        chrom, loc = bin_str.split("_")
        loc = int(loc)
        start = loc - bin_size
        return start, loc

    @staticmethod
    def prep_clustered_data_frame(clustered_data_frame):
        cluster_labels = clustered_data_frame['class']
        input_labels = clustered_data_frame['input']

        working_df = clustered_data_frame.drop(['class', 'input'], axis=1)

        cpgMatrix = np.array(working_df)
        cpgPositions = np.array([int(x) for x in working_df.columns])

        return cpgMatrix, cpgPositions, cluster_labels, input_labels

    @staticmethod
    def get_color(cpg_value: int):
        if cpg_value == 1:
            return "black"
        elif cpg_value == 0:
            return "white"

        else:
            NotImplementedError("I cannot yet accept unknown values. But I will soon.")

    def plotCpGMatrix(self, cpgMatrix, cpgPositions, title=None, figsize=(8, 8)):
        fig, ax = plt.subplots(figsize=figsize)
        v_steps = 1 / cpgMatrix.shape[0]
        v_spacings = np.arange(0, 1, v_steps)
        h_spacings = (cpgPositions - min(cpgPositions)) * 0.01
        ax.set_ylim(-.1, 1.1)
        ax.set_xlim(-.1, 1.1)
        ax.set_xticks(h_spacings)
        ax.set_xticklabels(cpgPositions, rotation=90)
        ax.set_yticks([])
        if title:
            ax.set_title(title)
        radius = min(v_steps / 2.5, 0.05)
        h_spacings = self.buffer_spacings(h_spacings, radius * 2)

        for read, vspace in zip(cpgMatrix, v_spacings):
            ax.axhline(vspace, color="black", zorder=1)
            for cpg, hspace in zip(read, h_spacings):
                x = hspace
                y = vspace
                circle = plt.Circle((x, y), radius=radius, facecolor=self.get_color(cpg), edgecolor="black")
                ax.add_artist(circle)

        return fig

    def plotCpGReadDataFrame(self, clustered_data_frame, title=None, figsize=(8, 8)):
        """
        Wrapper method to combine prepping a clustered data frame and plotting the matrix
        :param clustered_data_frame: #todo describe this
        :param figsize: tuple passed to matplotlib figsize (8,8)
        :return: matplotlib figure
        """
        prepped_data = self.prep_clustered_data_frame(clustered_data_frame)

        return self.plotCpGMatrix(prepped_data[0], prepped_data[1], title=title, figsize=figsize)
