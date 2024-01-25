"""Script for making agp

This script makes an AGP file with the right orientation and order of the medusa scaffolds.

Author: Evy Kuenen
Date: 8-1-2024

"""

from collections import defaultdict
import csv


def read_coordinate_file(file_path, skip_lines=5):
    """
    Read the coordinate file and return scaffold information.

    Parameters:
    - file_path: The path to the coordinate file.
    - skip_lines: The number of lines to skip at the beginning of the file.

    Returns:
    - A list of scaffold information.
    """
    data = []
    with open(file_path, "r") as coord_file:
        for _ in range(skip_lines):
            next(coord_file)
        for line in coord_file:
            fields = line.strip().split("\t")
            fields_final = [field.split(" | ") for field in fields]
            flat_fields_final = sum(fields_final, [])
            formatted_fields = [item.strip() for sublist in flat_fields_final for item in sublist.split()]
            data.append(formatted_fields)
    return data


def process_coordinate_data(data, output_orientation=True):
    """
    Process scaffold information and return a defaultdict grouped by a scaffold name.

    Parameters:
    - data: The coordinate data.
    - output_orientation: A flag indicating which orientation scaffold has.

    Returns:
    - A defaultdict grouped by scaffold
    """
    consolidated_data = defaultdict(list)
    for row in data:
        tag = row[12]  # scaffold
        covq = float(row[10])  # COVQ column
        idy = float(row[6])  # %IDY column
        if output_orientation:
            orientation = "-" if row[2] > row[3] else orientation = "+"
            row.append(orientation)
        consolidated_data[tag].append((covq, idy, row))
    return consolidated_data


def get_best_data(consolidated_data):
    """
    Get the best scaffold based on COVQ and %IDY.

    Parameters:
    - consolidated_data: The data grouped by scaffold

    Returns:
    - A list of the best scaffolds.
    """
    best_data = []
    for tag, data_list in consolidated_data.items():
        sorted_data = sorted(data_list, key=lambda x: (x[0], x[1]), reverse=True)
        best_entry = sorted_data[0][2]  # Get the row with the highest COVQ and %IDY
        best_data.append(best_entry)
    return best_data


def create_order_dict(order_scaffolds):
    """
    Create a dictionary to store the order of scaffolds.

    Parameters:
    - order_scaffolds: A list of scaffolds.

    Returns:
    - A dictionary with scaffold names as keys and their order as values.
    """
    return {name: index for index, name in enumerate(order_scaffolds)}


def sort_data_by_order(best_data, order_dict):
    """
    Sort data based on order in which contigs are in coord file.

    Parameters:
    - best_data: The data to be sorted.
    - order_dict: A dictionary with order information.

    Returns:
    - A list of sorted data.
    """
    return sorted(best_data, key=lambda item: (order_dict.get(item[12], float('inf')), item[2]))


def generate_agp_file(sorted_data, agp_file_path):
    """
    Generate an AGP file from sorted data and save it.

    Parameters:
    - sorted_data: The sorted data.
    - agp_file_path: The path to save the AGP file.
    """
    agp = []
    for item in sorted_data:
        data = [item[12], 1, item[8], 1, "W", item[11], 1, item[7], item[13]]
        agp.append(data)

    with open(agp_file_path, "w", newline='') as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        writer.writerows(agp)


def main():
    coord_reference_path = "../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsCurrentGenome/medusaOrderedToCurrentRefseq.coord"
    coord_contigs_path = "../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/medusaToContigs/medusaToContigs.coord"
    output_agp_path = "../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/makeAGP/medusa.agp"

    data_coord_reference = read_coordinate_file(coord_reference_path)
    consolidated_data_reference = process_coordinate_data(data_coord_reference)
    best_mpdis = get_best_data(consolidated_data_reference)

    order_scaffolds = [entry[12] for entry in best_mpdis]
    order_dict = create_order_dict(order_scaffolds)

    data_coord_contigs = read_coordinate_file(coord_contigs_path)
    consolidated_data_contigs = process_coordinate_data(data_coord_contigs, output_orientation=True)
    best_contigs = get_best_data(consolidated_data_contigs)

    sorted_coord_info = sort_data_by_order(best_contigs, order_dict)
    generate_agp_file(sorted_coord_info, output_agp_path)


main()
