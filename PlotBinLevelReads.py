from ParseBam import BamFileReadParser
import  matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import argparse
import logging
import sys
import os

sns.set_context("talk")
sns.set_style("darkgrid")


def plot_complete_bin_reads(matrix, chromosome: str, start_pos: int, stop_pos: int, output_loc: str):
    # Take a martix of CpG status and plot
    # Visualze the reads from the bam file, 1=methylated, 0=unmethylated
    plt.figure(figsize=(10, 6))
    g = sns.heatmap(matrix, vmax=1, vmin=0, cmap='coolwarm', linewidths=0.1)
    g.set_title("{}: {}-{}".format(chromosome, start_pos, stop_pos))
    g.set_ylabel("reads")
    g.set_xlabel("CpG site")
    file_name = "{}_{}.png".format(chromosome, str(stop_pos))
    output_file = os.path.join(output_loc, file_name)
    plt.savefig(output_file)


if __name__ == "__main__":
    # Get arguments from the command line
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("bins_to_plot",
                            help="File with each line being one bin to extract and plot, "
                                 "bins should be in format: chr19:3333333")
    arg_parser.add_argument("input_bam",
                            help="Input bam file, coordinate sorted with index present")
    arg_parser.add_argument("-o", "--output_dir",
                            help="Output directory to save figures, defaults to bam file loaction",)

    args = arg_parser.parse_args()

    # Set output directory to user input or default
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = os.path.dirname(args.input_bam)

    print("Output directory set to {}".format(output_dir))
    file_spec_path = os.path.join(output_dir, "plots_{}".format(os.path.basename(os.path.splitext(args.input_bam)[0])))

    # Create subdirectory for plots
    if not os.path.exists(file_spec_path):
        print("Creating file specific plots folder in output directory")
        os.makedirs(file_spec_path)
    else:
        print("plots folder already exists, saving there...")

    # write out inputs for reference
    print("Input bam file is: {}".format(args.input_bam))
    print("List of bins is: {}".format(args.bins_to_plot))
    print("Plots being saved to {}".format(file_spec_path))
    sys.stdout.flush()

    # Load bins into memory
    bins = []
    with open(args.bins_to_plot, 'r') as f:
        for line in f:
            bins.append(line.strip())

    # Create bam file parser object
    bam_parser = BamFileReadParser(args.input_bam, 20)

    # loop over bins generating a matrix for each
    for bin in bins:
        bin_chr = bin.split("_")
        chromosome = bin_chr[0]
        stop_pos = int(bin_chr[1])

        reads = bam_parser.parse_reads(chromosome, stop_pos-100, stop_pos)
        matrix = bam_parser.create_matrix(reads)
        matrix = matrix.dropna()

        # plot the matrix
        plot_complete_bin_reads(matrix, chromosome, stop_pos-100, stop_pos, file_spec_path)

