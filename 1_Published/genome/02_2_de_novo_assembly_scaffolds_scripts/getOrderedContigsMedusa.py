"""Script for ordering medusa scaffolds

This script orders the medusa scaffolds, finds the right orientation and 
outputs this in an agp file
Author: Evy Kuenen
Date: 8-1-2024

"""
import pandas as pd
from collections import defaultdict
import csv

def read_mapinfo_file(file_path, skip_lines=5):
    """
    Read a mapinfo file and extract all scaffold information.

    Parameters:
    - file_path: The path to the mapinfo file.
    - skip_lines: The number of lines to skip at the beginning of the file.

    Returns:
    - A list of lists containing scaffold information.
    """
    scaffolds = []
    with open(file_path, "r") as mapinfo_file:
        for _ in range(skip_lines):
            next(mapinfo_file)
        
        for line in mapinfo_file:
            data = []
            fields = line.strip().split("\t")
            fields_final = [field.split(" | ") for field in fields]
            flat_fields_final = sum(fields_final, [])
            formatted_fields = [item.strip() for sublist in flat_fields_final for item in sublist.split()]
            orientation = "-" if int(formatted_fields[2]) > int(formatted_fields[3]) else "+"
            data.extend([formatted_fields[0], formatted_fields[1], formatted_fields[3], formatted_fields[2],
                         formatted_fields[4], formatted_fields[5], formatted_fields[6], formatted_fields[7],
                         formatted_fields[8], formatted_fields[9], formatted_fields[10], formatted_fields[11],
                         formatted_fields[12], orientation])
            scaffolds.append(data)
    return scaffolds

def get_best_scaffolds(scaffolds):
    """
    Get the best scaffold information based on COVQ and %IDY.

    Parameters:
    - scaffolds: A list of scaffold information.

    Returns:
    - A list containing the best scaffold information.
    """
    consolidated_data = defaultdict(list)
    for row in scaffolds:
        tag = row[12]  # Last column is the scaffold
        covq = float(row[10])  # COVQ column
        idy = float(row[6])  # %IDY column
        consolidated_data[tag].append((covq, idy, row))

    best_scaffolds = []
    for tag, data_list in consolidated_data.items():
        # Sort based on COVQ in reverse order (highest first), and then on %IDY in reverse order
        sorted_data = sorted(data_list, key=lambda x: (x[0], x[1]), reverse=True)
        best_scaffold = sorted_data[0][2]  # Take the row with the highest COVQ and %IDY
        best_scaffolds.append(best_scaffold)

    return best_scaffolds

def sort_scaffolds(best_scaffolds):
    """
    Sort scaffolds based on chromosome and start location.

    Parameters:
    - best_scaffolds: A list of best scaffold information.

    Returns:
    - A sorted list of scaffold information.
    """
    return sorted(best_scaffolds, key=lambda x: (x[11].split("_")[1], int(x[0])))

def create_agp_file(sorted_scaffolds, output_file_path):
    """
    Create an AGP file from sorted scaffold information and save it.

    Parameters:
    - sorted_scaffolds: A sorted list of scaffold information.
    - output_file_path: The path to save the AGP file.
    """
    agp = []
    for line in sorted_scaffolds:
        agp_line = [line[12], 1, line[8], 1, "W", line[12], 1, line[8], line[13]]
        agp.append(agp_line)

    with open(output_file_path, "w") as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        writer.writerows(agp)

def main():
    mapinfo_file_path = "../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/medusaToCurrentRefseq/medusaToOld.coord"
    output_agp_file_path = "../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsCurrentGenome/orderedOrientedScaffolds.agp"

    scaffolds = read_mapinfo_file(mapinfo_file_path)
    best_scaffolds = get_best_scaffolds(scaffolds)
    sorted_scaffolds = sort_scaffolds(best_scaffolds)
    create_agp_file(sorted_scaffolds, output_agp_file_path)

main()
