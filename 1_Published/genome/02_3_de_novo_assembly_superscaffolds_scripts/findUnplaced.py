"""Script for to find unplaced contigs and put gaps in agp

This script finds unplaced contigs and puts them in the agp file and 
puts gaps between contigs.
Author: Evy Kuenen
Date: 8-1-2024

"""
import csv
import os

def read_agp_file(file_path):
    """
    Read AGP file and return a list of placed elements and their lengths.

    Parameters:
    - file_path: Path to the AGP file.

    Returns:
    - Tuple containing a list of placed elements and their total length.
    """
    placed_elements = []
    total_length = 0
    with open(file_path, "r") as agp_file:
        for line in agp_file:
            fields = line.strip().split("\t")
            if fields[4] == "W":
                placed_elements.append(fields[5])
                total_length += int(fields[7])
    return placed_elements, total_length

def read_unplaced_agp(allContigs, allContigsFlye, all_placed, agpverplaatst, mitoplaced, chloropplaced):
    """
    Read AGP file and append unplaced elements to the existing unplaced set.

    Parameters:
    - file_path: Path to the AGP file.
    - unplaced_set: Set containing unplaced elements.

    Returns:
    - List of placed and unplaced AGP data.
    """
    unplaced = []
    placed_and_unplaced_agp = []
    with open(agpverplaatst, "r") as verplaastsingenAgp:
        for line in verplaastsingenAgp:
            fields = line.strip().split("\t")
            placed_and_unplaced_agp.append(fields)
    with open(allContigs, "r") as agp_file:
        for line in agp_file:
            fields = line.strip().split("\t")
            if fields[5] not in all_placed: 
                unplaced.append(fields)
            
        for item in unplaced:
            if item[4] != 'N':
                data1 = ["unplaced_genome", 1, 1000, 1, "U", 1000, "scaffold", "yes", "align_genus"]
                data = ["unplaced_genome", 1, item[7], 1, "W", item[5], 1, item[7], item[8]]  # Use item[7] and item[8]
                placed_and_unplaced_agp.append(data)
                placed_and_unplaced_agp.append(data1)

    with open(mitoplaced, "r") as mitochondrien:
        for line in mitochondrien:
            fields = line.strip().split("\t")
            placed_and_unplaced_agp.append(fields)
    with open(chloropplaced, "r") as chloroplast:
        for line in chloroplast:
            fields = line.strip().split("\t")
            placed_and_unplaced_agp.append(fields)
    return placed_and_unplaced_agp

def write_agp_file(output_file_path, agp_data):
    """
    Write AGP data to an output file.

    Parameters:
    - output_file_path: Path to the output AGP file.
    - agp_data: List of AGP data.
    """
    with open(output_file_path, "w") as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        writer.writerows(agp_data)

def count_unplaced_mpdis(file_path):
    """
    Count the number of unplaced mpdis and their total length from the specified AGP file.

    Parameters:
    - file_path: Path to the AGP file.

    Returns:
    - Tuple containing the number of unplaced mpdis and their total length.
    """
    total_unplaced_mpdis_bp = 0
    unplaced_mpdis = set()
    with open(file_path, "r") as agp_file:
        for _ in range(5):
            next(agp_file)

        for line in agp_file:
            fields = line.strip().split("\t")
            fields_final = [field.split(" | ") for field in fields]
            flat_fields_final = sum(fields_final, [])
            formatted_fields = [item.strip() for sublist in flat_fields_final for item in sublist.split()]
            if formatted_fields[4] == "W":
                unplaced_mpdis.add(formatted_fields[5])
                total_unplaced_mpdis_bp += int(formatted_fields[7])
    return len(unplaced_mpdis), total_unplaced_mpdis_bp

def count_total_mpdis_length(file_path):
    """
    Count the total length of all mpdis from the specified file.

    Parameters:
    - file_path: Path to the file containing mpdi lengths.

    Returns:
    - Total length of all mpdis.
    """
    total_length = 0
    with open(file_path, "r") as mpdis_file:
        for line in mpdis_file:
            if line:
                parts = line.split(': ')
                if len(parts) > 1:
                    total_length += int(parts[2])
    return total_length

def main():
    agpverplaatst = "../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/100Verplaatsingen/100verplaatsingen.agp"
    mitoplaced = "../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/orderedContigsWithOrientationMitochondrienRefSeqNoDoubles.reindexed.agp"
    chloropplaced = "../../../data/chloroplast/mappingAllContigsThatMapTochloroplastToRefSeq/dontMapToRefSeq.agp"
    allContigs = "../../../data/genome/02_deNovoAssembly/scaffolds/makingAGPfile/flyeSspace.agp"
    output_file = "../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenWithUnplaced.agp"
    unplaced_mpdis = "../../../data/genome/02_deNovoAssembly/annotation_releases/annotation_releases/4686/100/GCF_001876935.1_Aspof.V1/GCF_001876935.1_Aspof.V1_assembly_structure/Primary_Assembly/unplaced_scaffolds/AGP/unplaced.scaf.agp"
    total_mpdis = "../../../data/genome/02_deNovoAssembly/annotation_releases/annotation_releases/4686/100/GCF_001876935.1_Aspof.V1/GCF_001876935.1_Aspof.V1_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP/combined_contig_lengths.txt"

    placed_genome, placed_bp_genome = read_agp_file(agpverplaatst)
    placed_mitochondrien, placed_bp_mitochondrien = read_agp_file(mitoplaced)
    placed_chloroplast, placed_bp_chloroplast = read_agp_file(chloropplaced)
    
    all_placed_bp = placed_bp_chloroplast + placed_bp_mitochondrien + placed_bp_genome
    all_placed = placed_chloroplast + placed_mitochondrien + placed_genome
    
    # print("total bp placed new genome: " + str(placed_bp_genome))

    bp_all_contigs, all_contigs = read_agp_file(allContigs)
    # print("total length bp new contigs without gaps: " + str(bp_all_contigs))
    # print("all placed new genome:" + str(len(all_placed)))
    # print("all contigs new genome:" + str(all_contigs))
    # print(len(placed_mitochondrien))
    print(len(placed_genome))
    print('----')
    print(len(all_placed))
   
    allContigsFlye, _ = read_agp_file(allContigs)

    placed_and_unplaced_agp = read_unplaced_agp(allContigs, allContigsFlye, all_placed, agpverplaatst, mitoplaced, chloropplaced)
    write_agp_file(output_file, placed_and_unplaced_agp)

    total_unplaced_mpdis, total_unplaced_mpdis_bp = count_unplaced_mpdis(unplaced_mpdis)
    print("total length unplaced mpdis in bp: " + str(total_unplaced_mpdis_bp))
    print("unplaced mpdis: " + str(total_unplaced_mpdis))
    print("total mpdis: 56182")

    total_lenght_bp_mpdis = count_total_mpdis_length(total_mpdis)
    print("total length mpdis: " + str(total_lenght_bp_mpdis))

main()
