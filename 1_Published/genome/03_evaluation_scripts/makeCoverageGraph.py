"""Script for making a coverage graph for every chromosome.

This script makes a coverage graph for every chromosome of the coverage of centromere regions on
every chromosome.

Author: Evy Kuenen
Date: 11-1-2024

"""
import matplotlib.pyplot as plt
import numpy as np

def read_coverage_data(file_path):
    """
    Read coverage data from a TSV file.
    Returns a list of lists containing the data.
    """
    coverage_data = []

    with open(file_path, "r") as tsv_file:
        for line in tsv_file:
            fields = line.strip().split("\t")
            coverage_data.append(fields)

    return coverage_data

def plot_coverage_per_chromosome(coverage_data, output_plot_path):
    """
    Plot coverage data per chromosome and save the plot to the specified file.
    """
    unique_chromosomes = set(item[0] for item in coverage_data)

    num_rows = len(unique_chromosomes)
    num_cols = 1

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(8, 4 * num_rows))

    for i, chromosome in enumerate(unique_chromosomes):
        x_values = []
        y_values = []

        for item in coverage_data:
            if item[0] == chromosome:
                x_values.append((int(item[1]) + int(item[2])) / 2)
                y_values.append(float(item[6]) * 100)

        axs[i].plot(x_values, y_values, label=f"{chromosome}")
        axs[i].set_title(f"{chromosome}")
        axs[i].set_xlabel('Position on chromosome in bp')
        axs[i].set_ylabel('Coverage in percentage')
        axs[i].set_ylim(0, 100)
        axs[i].get_xaxis().get_major_formatter().set_scientific(False)

    plt.tight_layout()
    plt.savefig(output_file_path)
    plt.show()

def main():
    input_file_path = "../../../data/genome/03_evaluation_genome/trFinder/coverage.txt"
    output_plot_path = "../../../data/genome/03_evaluation_genome/trFinder/coverage.png"

    coverage_data = read_coverage_data(input_file_path)
    plot_coverage_per_chromosome(coverage_data, output_plot_path)

main()