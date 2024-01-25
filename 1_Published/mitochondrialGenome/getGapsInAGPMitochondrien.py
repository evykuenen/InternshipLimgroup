"""Script for getting order and orientation of mitochondrial contigs and putting gaps in agp

This script gets from a coord file the orientation, the contig with the best coverage and identity on the mitochondrial genome.
These contigs are then filtered on more then 1 percent match on the reference genome and if two contigs map to the same place the
one with the longest coverge stays. Then in the agp that will be made gaps are put in between the contigs.

Author: Evy Kuenen
Date: 18-12-2023

"""
import csv
from collections import defaultdict


def read_coord_file(coord_file):
    """
    Reads a coordinate file and extracts scaffold information.

    Parameters:
    - coord_file (str): Path to the coordinate file.

    Returns:
    - scaffoldorder (list): List containing scaffold information.
    """
    scaffoldorder = []
    skip_lines = 5
    with open(coord_file, "r") as mapinfo_file:
        for _ in range(skip_lines):
            next(mapinfo_file) 
            
        for line in mapinfo_file:
            data = []
            fields = line.strip().split("\t")
            fields_final = [field.split(" | ") for field in fields]
            flat_fields_final = sum(fields_final, [])
            formatted_fields = [item.strip() for sublist in flat_fields_final for item in sublist.split()]
            if int(formatted_fields[2]) > int(formatted_fields[3]):  # get orientation
                data.append(formatted_fields[0])
                data.append(formatted_fields[1])
                data.append(formatted_fields[3])
                data.append(formatted_fields[2])
                data.append(formatted_fields[4])
                data.append(formatted_fields[5])
                data.append(formatted_fields[6])
                data.append(formatted_fields[7])
                data.append(formatted_fields[8])
                data.append(formatted_fields[9])
                data.append(formatted_fields[10])
                data.append(formatted_fields[11])
                data.append(formatted_fields[12])
                data.append("-")
            else:
                data.append(formatted_fields[0])
                data.append(formatted_fields[1])
                data.append(formatted_fields[2])
                data.append(formatted_fields[3])
                data.append(formatted_fields[4])
                data.append(formatted_fields[5])
                data.append(formatted_fields[6])
                data.append(formatted_fields[7])
                data.append(formatted_fields[8])
                data.append(formatted_fields[9])
                data.append(formatted_fields[10])
                data.append(formatted_fields[11])
                data.append(formatted_fields[12])
                data.append("+")
            scaffoldorder.append(data)
    return scaffoldorder


def get_best_contigs(scaffoldorder):
    """
    Consolidates scaffold data and selects scaffolds with the highest coverage and identity percentage.

    Parameters:
    - scaffoldorder (list): List containing scaffold information.

    Returns:
    - best_scaffolds (list): List of selected scaffolds.
    """
    consolidated_data = defaultdict(list)

    # loop through every line and consolidate
    for row in scaffoldorder:
        if float(row[10]): 
            tag = row[12]  # scaffold
            covq = float(row[10])  # COVQ-column
            idy = float(row[6])  # %IDY-column
            
            consolidated_data[tag].append((covq, idy, row))

    # get scaffolds with highest coverage and identity percentage
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


def group_data(best_scaffolds):
    """
    Groups selected scaffolds by scaffold name.

    Parameters:
    - best_scaffolds (list): List of selected scaffolds.

    Returns:
    - grouped_data (dict): Dictionary containing grouped scaffolds.
    """
    grouped_data = {}
    for item in best_scaffolds:
        scaffold_name = item[12]
        if scaffold_name not in grouped_data:
            grouped_data[scaffold_name] = []
        grouped_data[scaffold_name].append(item)
    return grouped_data


def select_rows(grouped_data):
    """
    Selects rows with the highest coverage and identity percentage for each scaffold.

    Parameters:
    - grouped_data (dict): Dictionary containing grouped scaffolds.

    Returns:
    - selected_rows (list): List of selected rows.
    """
    # select row with highest coverage and id percentage for every contig
    selected_rows = []
    for scaffold_name, scaffold_data in grouped_data.items():
        max_coverage_row = max(scaffold_data, key=lambda x: float(x[10]))
        selected_rows.append(max_coverage_row)
    selected_rows = sorted(selected_rows, key=lambda l: int(l[0]))
    return selected_rows


def get_longest_contigs(selected_rows):
    """
    Extracts information about the longest contigs.

    Parameters:
    - selected_rows (list): List of selected rows.

    Returns:
    - longest_contigs (list): List containing information about the longest contigs.
    """
    print(selected_rows)
    longest_contigs = []
    for data in selected_rows:
        contig_info = {
            'name': data[12],
            'start': int(data[0]),
            'end': int(data[1]),
            'length': int(data[8])
        }
        longest_contigs.append(contig_info)
    return longest_contigs


def filter_contigs(longest_contigs, selected_rows):
    """
    Filters overlapping contigs and removes contigs that are within other contigs from the list.

    Parameters:
    - longest_contigs (list): List containing information about the longest contigs.

    Returns:
    - filteredOnMatchPercentage (list): List of filtered contigs based on match percentage.
    """
    filtered_contigs = []
    for i, current_contig in enumerate(longest_contigs):
        overlapping = False
        for j, other_contig in enumerate(longest_contigs):
            if i != j:  # skip current contig
                if (
                    current_contig['start'] >= other_contig['start']
                    and current_contig['end'] <= other_contig['end']
                ):
                    overlapping = True
                    break  # if there is overlap stop loop
        if not overlapping:
            filtered_contigs.append(current_contig)

    notDoubleContigs = []
    for item in selected_rows:
        for filtered_contig in filtered_contigs:
            if filtered_contig['name'] == item[12]:
                notDoubleContigs.append(item)

    filteredOnMatchPercentage = []
    for line in notDoubleContigs:
        if float(line[4]) > float(int(line[8]) / 100):
            filteredOnMatchPercentage.append(line)
    return filteredOnMatchPercentage


def get_info_from_agp(agp_file):
    """
    Reads information from an AGP file.

    Parameters:
    - agp_file (str): Path to the AGP file.

    Returns:
    - flyeAgplist (list): List containing information from the AGP file.
    """
    flyeAgplist = []
    with open(agp_file, "r") as flyeSspace_agp:
        for line in flyeSspace_agp:
            linelist = []
            fields = line.strip().split("\t")
            for items in fields:
                linelist.append(items)
            flyeAgplist.append(linelist)
    return flyeAgplist


def compare_to_agp(flyeAgplist, filteredOnMatchPercentage):
    """
    Compares scaffold names and appends information from AGP file.

    Parameters:
    - flyeAgplist (list): List containing information from the AGP file.
    - filteredOnMatchPercentage (list): List of filtered contigs based on match percentage.

    Returns:
    - finallist (list): List containing final information.
    """
    finallist = []
    idlistPlacedScaffolds = []
    for item in filteredOnMatchPercentage:
        for lijst in flyeAgplist:
            if lijst[5] == item[12] and lijst[4] != 'N':
                print("items filterOnmatch:", item)
                finallinelist = []
                finallinelist.append('mitochondrial')
                idlistPlacedScaffolds.append(item[12])
                finallinelist.append(item[0])
                finallinelist.append(item[1])
                finallinelist.append(lijst[3])
                finallinelist.append(lijst[4])
                finallinelist.append(lijst[5])
                finallinelist.append(lijst[6])
                finallinelist.append(lijst[7])
                finallinelist.append(item[13])
                finallist.append(finallinelist)
    return finallist


def put_gaps_in_agp(finallist):
    """
    Inserts gap information into the AGP (A Golden Path) based on the provided list of scaffold information.

    Parameters:
    - finallist (list): List containing scaffold information.

    Returns:
    - output_with_gaps (list): List with gap information included in the AGP.
    """
    output_with_gaps = []
    previous_end = None

    for line in finallist:
        if previous_end is not None:
            current_start = int(line[1]) 

            gap_length = current_start - previous_end - 1
            if gap_length > 0:
                gap_line = ['mitochondrial', str(previous_end + 1), str(current_start - 1), '1', 'N', str(gap_length), 'scaffold', 'yes', 'align_genus', 'no']
                output_with_gaps.append(gap_line)
            if gap_length == 0 or gap_length < 0:
                gap_line = ['mitochondrial', str(previous_end + 1), str(current_start - 1), '1', 'N', str(1), 'scaffold', 'yes', 'align_genus', 'no']
                output_with_gaps.append(gap_line)
            
        output_with_gaps.append(line)
        previous_end = int(line[2])
    return output_with_gaps


def write_to_file(output_with_gaps, output_file):
    """
    Writes the final AGP information with gaps to an output file.

    Parameters:
    - output_with_gaps (list): List containing scaffold information with gaps.
    - output_file (str): Path to the output file.

    Returns:
    - None
    """
    with open(output_file, "w") as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        writer.writerows(output_with_gaps)


def main():
    """
    Main function that orchestrates the execution of the entire process.

    Returns:
    - None
    """
    coord_file = "../../../data/mitochondrien/contigsToMitochondrien/allContigstomitochondrien.coord"
    scaffoldorder = read_coord_file(coord_file)
    best_scaffolds = get_best_contigs(scaffoldorder)
    grouped_data = group_data(best_scaffolds)
    selected_rows = select_rows(grouped_data)
    longest_contigs = get_longest_contigs(selected_rows)
    filteredOnMatchPercentage = filter_contigs(longest_contigs, selected_rows)
    agp_file = "../../../data/genome/02_deNovoAssembly/scaffolds/makingAGPfile/flyeSspace.agp"
    flyeAgplist = get_info_from_agp(agp_file)
    finallist = compare_to_agp(flyeAgplist, filteredOnMatchPercentage)
    output_with_gaps = put_gaps_in_agp(finallist)
    output_file = "../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/orderedContigsWithOrientationMitochondrienRefSeqNoDoubles.agp"
    write_to_file(output_with_gaps, output_file)


main()
