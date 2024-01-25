"""Script for getting order and orientation of mitochondrial contigs

This script gets from a coord file the orientation, the contig with the best coverage and identity on the mitochondrial genome.
These contigs are then filtered on more than 1 percent match on the reference genome and if two contigs map to the same place the
one with the longest coverge stays.

Author: Evy Kuenen
Date: 29-11-2023

"""
import csv
from collections import defaultdict


def read_mapinfo_file(file_path, skip_lines):
    """
    Reads the mapinfo file and returns all data from coord file.

    Args:
        file_path (str): Path to the mapinfo file.
        skip_lines (int): Number of lines to skip at the beginning of the file.

    Returns:
        list: List of data from the mapinfo file.
    """
    scaffoldorder = []
    with open(file_path, "r") as mapinfo_file:
        for _ in range(skip_lines):
            next(mapinfo_file)

        for line in mapinfo_file:
            data = []
            fields = line.strip().split("\t")
            fields_final = [field.split(" | ") for field in fields]
            flat_fields_final = sum(fields_final, [])
            formatted_fields = [item.strip() for sublist in flat_fields_final for item in sublist.split()]
            data.append(formatted_fields)
            scaffoldorder.append(data)

    return scaffoldorder


def get_oriented_data(formatted_fields):
    """
    Gets oriented based on start and end position in coord file.

    Args:
        formatted_fields (list): List of formatted fields from the mapinfo file.

    Returns:
        list: Oriented data.
    """
    data = []
    if int(formatted_fields[2]) > int(formatted_fields[3]):
        data.extend([formatted_fields[0], formatted_fields[1], formatted_fields[3], formatted_fields[2]])
        data.extend(formatted_fields[4:13])
        data.append("-")
    else:
        data.extend([formatted_fields[0], formatted_fields[1], formatted_fields[2], formatted_fields[3]])
        data.extend(formatted_fields[4:13])
        data.append("+")
    return data


def consolidate_data(scaffoldorder):
    """
    Consolidates data based on coverage and percent identity.

    Args:
        scaffoldorder (list): List of data from the mapinfo file.

    Returns:
        dict: Consolidated data.
    """
    consolidated_data = defaultdict(list)
    for row in scaffoldorder:
        if float(row[10]):
            tag = row[12]  # contigscaffold
            covq = float(row[10])  # COVQ-column
            idy = float(row[6])  # %IDY-column
            consolidated_data[tag].append((covq, idy, row))

    return consolidated_data


def get_best_scaffolds(consolidated_data):
    """
    Gets scaffolds with the highest coverage and identity percentage.

    Args:
        consolidated_data (dict): Consolidated data.

    Returns:
        list: Best scaffolds.
    """
    best_scaffolds = []
    processed_headers = set()
    for tag, data_list in consolidated_data.items():
        if tag in processed_headers:
            continue
        sorted_data = sorted(data_list, key=lambda x: (x[0], x[1]), reverse=True)
        best_scaffold = sorted_data[0][2]
        best_scaffolds.append(best_scaffold)
        processed_headers.add(tag)

    return best_scaffolds


def get_grouped_data(best_scaffolds):
    """
    Groups data based on scaffold name.

    Args:
        best_scaffolds (list): Best scaffolds.

    Returns:
        dict: Grouped data.
    """
    grouped_data = {}
    for item in best_scaffolds:
        scaffold_name = item[12]
        if scaffold_name not in grouped_data:
            grouped_data[scaffold_name] = []
        grouped_data[scaffold_name].append(item)
    return grouped_data


def get_row_highest_cov_and_id(grouped_data):
    """
    Gets the row with the highest coverage and percent id for every contig.

    Args:
        grouped_data (dict): Grouped data.

    Returns:
        list: Selected rows.
    """
    selected_rows = []
    for scaffold_name, scaffold_data in grouped_data.items():
        max_coverage_row = max(scaffold_data, key=lambda x: float(x[10]))
        selected_rows.append(max_coverage_row)

    selected_rows = sorted(selected_rows, key=lambda l: int(l[0]))

    print(len(selected_rows))
    return selected_rows


def filter_longest_contigs(selected_rows):
    """
    Filter overlapping contigs and remove contigs that are inside other contigs.

    Args:
        selected_rows (list): Selected rows.

    Returns:
        list: Filtered contigs.
    """
    longest_contigs = []

    for data in selected_rows:
        contig_info = {
            'name': data[12],
            'start': int(data[0]),
            'end': int(data[1]),
            'length': int(data[8])
        }
        longest_contigs.append(contig_info)

    # Filter overlapping contigs and remove contigs that are inside other contigs
    filtered_contigs = []
    for i, current_contig in enumerate(longest_contigs):
        overlapping = False
        for j, other_contig in enumerate(longest_contigs):
            if i != j:  # Skip current contig
                if (
                    current_contig['start'] >= other_contig['start']
                    and current_contig['end'] <= other_contig['end']
                ):
                    overlapping = True
                    break  # If there is an overlap, stop the loop
        if not overlapping:
            filtered_contigs.append(current_contig)

    return filtered_contigs


def not_double_contigs(selected_rows, filtered_contigs):
    """
    Filter out double contigs based on their name.

    Args:
        selected_rows (list): Selected rows.
        filtered_contigs (list): Filtered contigs.

    Returns:
        list: Contigs without duplicates.
    """
    notDoubleContigs = []
    for item in selected_rows:
        for filtered_contig in filtered_contigs:
            if filtered_contig['name'] == item[12]:
                notDoubleContigs.append(item)
    return notDoubleContigs


def filter_on_match_percentage(notDoubleContigs):
    """
    Filter contigs based on more then 1 percent match percentage of the contig on the reference genome.

    Args:
        notDoubleContigs (list): Contigs without duplicates.

    Returns:
        list: Filtered contigs based on match percentage.
    """
    filteredOnMatchPercentage = []
    for line in notDoubleContigs:
        if float(line[4]) > float(int(line[8]) / 100):
            filteredOnMatchPercentage.append(line)

    return filteredOnMatchPercentage


def read_flye_agp(file_path):
    """
    Read Flye AGP file.

    Args:
        file_path (str): Path to the Flye AGP file.

    Returns:
        list: List of data from the Flye AGP file.
    """
    flyeAgplist = []
    with open(file_path, "r") as flyeSspace_agp:
        for line in flyeSspace_agp:
            linelist = line.strip().split("\t")
            flyeAgplist.append(linelist)

    return flyeAgplist


def compare_and_append_agp_info(filteredOnMatchPercentage, flyeAgplist):
    """
    Get contig information from AGP to make new AGP.

    Args:
        filteredOnMatchPercentage (list): Filtered contigs based on match percentage.
        flyeAgplist (list): List of data from the Flye AGP file.

    Returns:
        list: Final list with appended Flye AGP information.
    """
    finallist = []
    idlistPlacedScaffolds = []
    for item in filteredOnMatchPercentage:
        for lijst in flyeAgplist:
            if lijst[5] == item[12] and lijst[4] != 'N':
                finallinelist = [
                    item[12], lijst[1], lijst[2], lijst[3], lijst[4],
                    lijst[5], lijst[6], lijst[7], item[13]
                ]
                idlistPlacedScaffolds.append(item[12])
                finallist.append(finallinelist)
    return finallist


def write_to_output_file(output_file_path, data):
    """
    Write data to the output file.

    Args:
        output_file_path (str): Path to the output file.
        data (list): Data to be written to the output file.
    """
    with open(output_file_path, "w") as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        writer.writerows(data)


def main():
    mapinfo_file_path = "../../../data/mitochondrien/contigsToMitochondrien/allContigstomitochondrien.coord"
    skip_lines = 5
    scaffoldorder = read_mapinfo_file(mapinfo_file_path, skip_lines)
    consolidated_data = consolidate_data(scaffoldorder)
    best_scaffolds = get_best_scaffolds(consolidated_data)
    grouped_data = get_grouped_data(best_scaffolds)
    selected_rows = get_row_highest_cov_and_id(grouped_data)
    filtered_contigs = filter_longest_contigs(selected_rows)
    notDoubleContigs = not_double_contigs(selected_rows, filtered_contigs)
    filteredOnMatchPercentage = filter_on_match_percentage(notDoubleContigs)
    flye_agp_path = "../../../data/genome/02_deNovoAssembly/scaffolds/makingAGPfile/flyeSspace.agp"
    flyeAgplist = read_flye_agp(flye_agp_path)
    finallist = compare_and_append_agp_info(filteredOnMatchPercentage, flyeAgplist)
    output_file_path = "../../../data/mitochondrien/contigsToMitochondrien/contigsVsOldRef/orderedContigsWithOrientationMitochondrienAllContigsNoDoubles.agp"
    write_to_output_file(output_file_path, finallist)


main()
