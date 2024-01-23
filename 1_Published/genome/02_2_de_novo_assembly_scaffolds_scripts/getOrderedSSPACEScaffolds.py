"""Script for ordering scaffolds

This script orders all scaffolds from SSPACE in order of where they map on the reference genome by using
a coord file and outputting a fasta file that is ordered.

Author: Evy Kuenen
Date: 19-12-2023

"""

import pandas as pd
from collections import defaultdict

def read_mapping_info(mapping_info_file_path):
    """
    Read coord file and return scaffold information.

    Args:
        file_path (str): Path to the mapping information file.

    Returns:
        list: Formatted data from the mapping information file.
    """
    data = []

    with open(file_path, "r") as mapinfo_file:
        for _ in range(6):
            next(mapinfo_file)

        for line in mapinfo_file:
            fields = line.strip().split("\t")
            fields_final = [field.split(" | ") for field in fields]
            flat_fields_final = sum(fields_final, [])
            formatted_fields = [item.strip() for sublist in flat_fields_final for item in sublist.split()]
            data.append(formatted_fields)

    return mapping_info_data

def consolidate_data(mapping_info_data):
    """
    Consolidate mapping information data based on the TAGS column (MPDI).

    Args:
        data (list): Formatted data from the mapping information file.

    Returns:
        defaultdict: Consolidated data with TAGS as keys.
    """
    consolidated_data = defaultdict(list)

    for row in data:
        tag = row[12]  # Last column is the TAGS column (MPDI)
        covq = float(row[10])  # COVQ column
        idy = float(row[6])  # %IDY column

        consolidated_data[tag].append((covq, idy, row))

    return consolidated_data_result

def get_best_scaffolds(consolidated_data_result):
    """
    Get the best scaffold for each TAGS (MPDI) based on COVQ and %IDY.

    Args:
        consolidated_data (defaultdict): Consolidated data with TAGS as keys.

    Returns:
        list: Best scaffolds based on COVQ and %IDY.
    """
    best_scaffolds = []

    for tag, data_list in consolidated_data.items():
        sorted_data = sorted(data_list, key=lambda x: (x[0], x[1]), reverse=True)
        best_scaffold = sorted_data[0][2]  # Select the row with the highest COVQ and %IDY
        best_scaffolds.append(best_scaffold)

    return best_scaffolds

def write_best_scaffolds_to_file(best_scaffolds, output_file_path, fasta_file_path):
    """
    Write the best scaffolds to an output file using information from a FASTA file.

    Args:
        best_scaffolds (list): Best scaffolds based on COVQ and %IDY.
        output_file_path (str): Path to the output file.
        fasta_file_path (str): Path to the FASTA file containing scaffold information.
    """
    scaffolds = set()

    with open(output_file_path, "w") as output_file:
        for scaffold in best_scaffolds:
            scaffold_name = scaffold[11].split()[0]
            if scaffold_name not in scaffolds:
                output_file.write(f">{scaffold[11]}\n")
                scaffolds.add(scaffold_name)

                with open(fasta_file_path, "r") as fasta_file:
                    record = False
                    for line in fasta_file:
                        if line.startswith(f">{scaffold_name}"):
                            record = True
                            continue
                        if record:
                            if line.startswith(">"):
                                break
                            output_file.write(line)

def main():
    mapping_info_file_path = "../../../data/genome/02_deNovoAssembly/superScaffolds/newScaffoldOrContigsAgianstRefSeq/mappingInfoNewContigsV2standard_output.txt"
    output_file_path = "../../../data/genome/02_deNovoAssembly/superScaffolds/OrderedScaffoldsPythonScript/fastascaffoldsordened.fa"
    fasta_file_path = "../../../data/genome/02_deNovoAssembly/scaffolds/SSPACE/v3GoodOutput/standard_output.final.scaffolds_andere_nummering.fasta"
    mapping_info_data = read_mapping_info(mapping_info_file_path)
    consolidated_data_result = consolidate_data(mapping_info_data)
    best_scaffolds_result = get_best_scaffolds(consolidated_data_result)
    write_best_scaffolds_to_file(best_scaffolds_result, output_file_path, fasta_file_path)

main()