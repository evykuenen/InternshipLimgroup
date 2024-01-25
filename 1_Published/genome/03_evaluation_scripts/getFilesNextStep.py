"""Script for extracting onemapbed and weights file 

This script makes from the output files of extractAndMapMarkers a new bed file and a weights file
in the bed file are the markers and their new location on the new genome where they map and the 
weights file is now for every marker a weight of 1 because there is only one type of marker.
Author: Evy Kuenen
Date: 11-1-2024

"""
import csv


def read_data(file_path):
    """
    Read data from the part 1 file and return a list of lists.

    Args:
    - file_path (str): The path to the file to be read.

    Returns:
    - list: A list of lists representing the data read from the file.
    """
    data = []
    with open(file_path, "r") as file:
        for line in file:
            fields = line.strip().split("\t")
            data.append(fields)
    return data


def process_data(part1_data, part2_data):
    """
    Process the input data and create a new list based on specific conditions.

    Args:
    - part1_data (list): List of data from the first file.
    - part2_data (list): List of data from the second file.

    Returns:
    - list: A processed list based on specific conditions.
    """
    onemapbed = []
    for item in part1_data:
        for line in part2_data:
            if item[0] == line[0]:
                data = [
                    item[1],
                    int(item[2]) - 1,
                    item[2],
                    item[0].split("_")[0] + "-" + line[1] + ":" + line[2],
                    item[1] + ":" + item[2]
                ]
                onemapbed.append(data)
    return onemapbed


def write_to_file(data, file_path, delimiter='\t'):
    """
    Write the provided data to a file.

    Args:
    - data (list): List of data to be written to the file.
    - file_path (str): The path to the file to be written.
    - delimiter (str): The delimiter to be used in the CSV file. Default is '\t'.
    """
    with open(file_path, "w") as output_file:
        writer = csv.writer(output_file, delimiter=delimiter)
        writer.writerows(data)


def main():
    part1_file_path = '../../../data/genome/allmapsGetContigOrderCurrentContig/1_onemap/scaffolds_adjusted_markers_part1.txt'
    part2_file_path = '../../../data/genome/allmapsGetContigOrderCurrentContig/1_onemap/scaffolds_adjusted_markers_part2.txt'

    part1_data = read_data(part1_file_path)
    part2_data = read_data(part2_file_path)
    onemapbed = process_data(part1_data, part2_data)
    onemapbed_output_path = '../../../data/genome/allmapsGetContigOrderCurrentContig/3_allmaps/scaffolds_mapped_onemap-cM.bed'
    weights_output_path = '../../../data/genome/allmapsGetContigOrderCurrentContig/3_allmaps/weights.txt'

    write_to_file(onemapbed, onemapbed_output_path)
    weights = [['AoV1', '1']]
    write_to_file(weights, weights_output_path)


main()
