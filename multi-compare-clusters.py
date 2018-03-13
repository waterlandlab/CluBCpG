import pandas as pd
from multiprocessing import Pool
import pickle
import sys


class Analysis:
    """
    A class to allow multiprocessing of the cluster comparions. Initializes with all needed data.
    """

    def __init__(self, overlapping_A, overlapping_B, unique_bins):
        self.overlapping_clusters_a = overlapping_A
        self.overlapping_clusters_b = overlapping_B
        self.unique_bins = unique_bins

    def get_comparisons(self, bin: str):
        As_temp = []
        Bs_temp = []
        bins_temp = []
        A = self.overlapping_clusters_a[self.overlapping_clusters_a['bin'] == bin]
        B = self.overlapping_clusters_b[self.overlapping_clusters_b['bin'] == bin]
        for index_a, cluster_a in A.iterrows():
            for index_b, cluster_b in B.iterrows():
                As_temp.append(cluster_a['methylation'])
                Bs_temp.append(cluster_b['methylation'])
                bins_temp.append(bin)

        return bins_temp, As_temp, Bs_temp

    def process(self):
        """
        Call this method to start multiprocessing across multiple cores
        :return: List of tuples which will need unpacked
        """
        pool = Pool(processes=24) # 24 for Sphere
        results = pool.map(self.get_comparisons, self.unique_bins)

        return results


if __name__ == "__main__":

    df = pd.read_csv(sys.argv[1], header=0, index_col=None)
    print("Input df len: ", len(df))

    # Filter for at least 4 reads per cluster
    df_prefilter = df[df['read_number'] >= 4]
    df_prefilter = df[df['class_label'] != -1]
    
    print("Pre-filetred df len: ", len(df_prefilter))

    all_results = []

    for i in range(2, 14):
        print("Cpg Number: {}".format(i))
        d = df_prefilter[df_prefilter['cpg_number'] == i]
        a_clusters = d[d.input_label == 'A']
        b_clusters = d[d.input_label == 'B']

        # for each label, find clusters common to both,
        overlapping_clusters_a = a_clusters[a_clusters['bin'].isin(b_clusters['bin'])]
        overlapping_clusters_b = b_clusters[b_clusters['bin'].isin(a_clusters['bin'])]

        # Get bins present in both
        unique_bins = overlapping_clusters_a['bin'].unique()

        # Sanity check
        if 'chr1_59119400' in unique_bins:
            print("WHY IS THIS BIN HERE????")
            sys.exit(13)

        print("Starting pool...")
        analyze = Analysis(overlapping_clusters_a, overlapping_clusters_b, unique_bins)
        results = analyze.process()
        print("Pool complete...")

        bins = []
        As = []
        Bs = []

        for result in results:
            for item in result[0]:
                bins.append(item)
            for item in result[1]:
                As.append(item)
            for item in result[2]:
                Bs.append(item)


        df_out = pd.DataFrame({
            "bin": bins,
            "A": As,
            "B": Bs,
        })

        df_out.to_csv("cluster_comparisons_{}cpgs.csv.gz".format(i), compression='gzip', index=False)

    print("Done")


