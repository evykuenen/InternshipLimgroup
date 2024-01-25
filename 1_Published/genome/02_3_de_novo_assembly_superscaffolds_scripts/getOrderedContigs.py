"""Script for ordering Contigs

This script orders all Flye contigs in order of the reference genome of Harkess, 
orients them and filters them. This script outputs an AGP file.

Author: Evy Kuenen
Date: 19-12-2023

"""

import csv
from collections import defaultdict
from itertools import groupby


def read_contigorder(file_path):
    """
    Read contig order information from coord file and return contig information.

    Args:
        file_path (str): Path to the contig order information file.

    Returns:
        list: Formatted data of contig order information.
    """
    contigorder = [] 
    skip_lines = 5

    with open(file_path, "r") as mapinfo_file:
        for _ in range(skip_lines):
            next(mapinfo_file)

        for line in mapinfo_file:
            data = []
            fields = line.strip().split("\t")
            fields_final = [field.split(" | ") for field in fields]
            flat_fields_final = sum(fields_final, [])
            formatted_fields = [item.strip() for sublist in flat_fields_final for item in sublist.split()]
            if int(formatted_fields[2]) > int(formatted_fields[3]):  # get right orientation
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
            contigorder.append(data)
    return contigorder 


def get_right_orientatien(contigorder):
    """
    Determine the correct orientation for contigs based on length and orientation.

    Args:
        contigorder (list): Formatted contig order information.

    Returns:
        list: Contig order information with correct orientation.
    """
    orientation_mapping = {}

    contigorder.sort(key=lambda x: x[11])
    grouped_contigs = {key: list(group) for key, group in groupby(contigorder, key=lambda x: x[11])}

    for key, group in grouped_contigs.items():
        values_list = [(item[4], item[13]) for item in group]
        if key not in orientation_mapping:
            orientation_mapping[key] = values_list
        else:
            orientation_mapping[key].extend(values_list)
    total_lengths_per_contig = {}

    for contig, orientations in orientation_mapping.items():
        total_length_minus = sum(int(length) for length, orientation in orientations if orientation == '-')
        total_length_plus = sum(int(length) for length, orientation in orientations if orientation == '+')
        
        total_lengths_per_contig[contig] = {
            'minus': total_length_minus,
            'plus': total_length_plus
        }

    for contig, lengths in total_lengths_per_contig.items():
        percentage_plus = (lengths["plus"] / (lengths["minus"] + lengths["plus"])) * 100
        if percentage_plus > 50:
            for item in contigorder:
                if item[11] == contig:
                    item[13] = '+'  # verander de oriëntatie naar '+'
        if percentage_plus < 50:
            for item in contigorder:
                if item[11] == contig:
                    item[13] = '-'  # verander de oriëntatie naar '-'

    orientated_contig_order = contigorder  # moet erin blijven
    return orientated_contig_order


def consolidate_data(orientated_contig_order):
    """
    Consolidate contig order data based on contig identifiers.

    Args:
        orientated_contig_order (list): Contig order information with correct orientation.

    Returns:
        defaultdict: Consolidated data with contig identifiers as keys.
    """
    consolidated_data = defaultdict(list)
    for row in orientated_contig_order:  # juiste orientatie
        tag = row[11]  # contig
        covq = float(row[10]) 
        idy = float(row[6])
        consolidated_data[tag].append((covq, idy, row))
    return consolidated_data


def get_best_contigs(consolidated_data):
    """
    Get the best contig for each contig identifier based on coverage and identity percentage.

    Args:
        consolidated_data (defaultdict): Consolidated data with contig identifiers as keys.

    Returns:
        list: Best contigs based on coverage and identity percentage.
    """
    best_contigs = []
    for tag, data_list in consolidated_data.items():  
        sorted_data = sorted(data_list, key=lambda x: (x[0], x[1]), reverse=True)  # eerst coverage dan id%
        best_contig = sorted_data[0][2]
        best_contigs.append(best_contig)
    best_contigs = sorted(best_contigs, key=lambda x: float(x[12][4:]))
    return best_contigs


def filter_on_match_percentage(best_contigs):
    """
    Filter contigs based on the match percentage criterion.

    Args:
        best_contigs (list): Best contigs based on coverage and identity percentage.

    Returns:
        list: Contigs filtered based on the match percentage criterion.
    """
    filteredOnMatchPercentage = []
    for line in best_contigs:
        if float(line[4]) > float(int(line[8]) / 100):
            filteredOnMatchPercentage.append(line)
    return filteredOnMatchPercentage


def read_flye_agp(file_path):
    """
    Read AGP (A Golden Path) information from a file.

    Args:
        file_path (str): Path to the AGP file.

    Returns:
        list: Formatted AGP information.
    """
    flyeAgplist = []
    with open(file_path, "r") as flyeSspace_agp:
        for line in flyeSspace_agp:
            linelist = []
            fields = line.strip().split("\t")
            for items in fields:
                linelist.append(items)
            flyeAgplist.append(linelist)
    return flyeAgplist


def process_finallist(flyeAgplist, filteredOnMatchPercentage):
    """
    Process the final list of contigs based on AGP information and match percentage filtering.

    Args:
        flyeAgplist (list): Formatted AGP information.
        filteredOnMatchPercentage (list): Contigs filtered based on match percentage.

    Returns:
        list: Final list of contigs.
    """
    finallist = []
    idlistPlacedScaffolds = []
    for item in filteredOnMatchPercentage:
        for lijst in flyeAgplist:
            if lijst[5] == item[11] and item[11] not in idlistPlacedScaffolds:
                print(item)
                print(idlistPlacedScaffolds)
                finallinelist = []
                finallinelist.append(item[11])
                idlistPlacedScaffolds.append(item[11])
                finallinelist.append(1)
                finallinelist.append(item[7])
                finallinelist.append(1)
                finallinelist.append('W')
                finallinelist.append(item[11])
                finallinelist.append(1)
                finallinelist.append(item[7])
                finallinelist.append(item[13])
                finallist.append(finallinelist)
    unplaced = []
    for lijst in flyeAgplist:
        if lijst[4] != "N":
            if lijst[5] not in idlistPlacedScaffolds:
                finallinelist = []
                finallinelist.append(lijst[5])
                idlistPlacedScaffolds.append(lijst[5])
                finallinelist.append(1)
                finallinelist.append(lijst[7])
                finallinelist.append(1)
                finallinelist.append('W')
                finallinelist.append(lijst[5])
                finallinelist.append(1)
                finallinelist.append(lijst[7])
                finallinelist.append(lijst[8])
                unplaced.append(finallinelist)
    finallist = finallist + unplaced
    return finallist


def write_agp_to_file(finallist, output_file_path):
    """
    Write the final list of contigs to an AGP file.

    Args:
        finallist (list): Final list of contigs.
        output_file_path (str): Path to the output AGP file.
    """
    with open(output_file_path, "w") as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        writer.writerows(finallist)
        

def main():
    contigorder = read_contigorder("../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq.coord")
    orientated_contig_order = get_right_orientatien(contigorder)
    consolidated_data = consolidate_data(orientated_contig_order)
    ordered_best = get_best_contigs(consolidated_data)
    filteredOnMatchPercentage = filter_on_match_percentage(ordered_best)
    flyeAgplist = read_flye_agp("../../../data/genome/02_deNovoAssembly/scaffolds/makingAGPfile/flyeSspace.agp")

    finallist = process_finallist(flyeAgplist, filteredOnMatchPercentage)
    output_file_path = "../../../data/genome/02_deNovoAssembly/superScaffolds/getOrderdFilteredContigs/ContigsOrderedFilteredToRefseq/ContigsOrderedFilteredToRefseqContigs.agp"
    write_agp_to_file(finallist, output_file_path)


main()
