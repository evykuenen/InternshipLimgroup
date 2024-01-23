"""Script for making a graph of centromere coverage in chromosomes

This script makes for every chromosome a graph of the coverage of centromere repeats on the chromosomes.
Author: Evy Kuenen
Date: 11-1-2024

"""
import matplotlib.pyplot as plt
import numpy as np

def read_tsv(input_file_path, line_skip=1):
    """
    Read TSV data from a file, skipping the specified number of lines.
    Returns a list of lists containing the data.
    """
    tsv_data = []

    with open(file_path, "r") as tsv_file:
        for _ in range(line_skip):
            next(tsv_file)

        for line in tsv_file:
            fields = line.strip().split("\t")
            tsv_data.append(fields)

    return tsv_data

def filter_and_plot_tsv(tsv_data, output_plot_path):
    """
    Filter and plot TSV data and save the plot to the specified file.
    """
    unique_item0_values = set(item[0] for item in tsv_data)

    num_rows = len(unique_item0_values)
    num_cols = 1

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(8, 4 * num_rows))

    for i, item0_value in enumerate(unique_item0_values):
        x_values = []
        y_values = []

        for item in tsv_data:
            if item[0] == item0_value:
                total_value = int(item[2]) + int(item[3])

                if total_value < 75:
                    continue

                x_values.append(int(item[1]))
                y_values.append(total_value)

        axs[i].plot(x_values, y_values, label=item0_value)
        axs[i].set_title(item0_value)
        axs[i].set_xlabel('positie in bp')
        axs[i].set_ylabel('hoeveelheid TTAGGG')
        axs[i].set_ylim(0, 1200)

    plt.tight_layout()
    plt.savefig(output_file_path)
    plt.show()

def main():
    input_file_path = "../../../data/genome/03_evaluation_genome/tidk/asparagusV2_telomeric_repeat_windows.tsv"
    output_plot_path = "../../../data/genome/03_evaluation_genome/tidk/tidk_filtered_subplots.png"

    # Read TSV data
    tsv_data = read_tsv(input_file_path)

    # Filter and plot TSV data
    filter_and_plot_tsv(tsv_data, output_plot_path)
main()