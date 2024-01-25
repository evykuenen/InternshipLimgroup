"""Script for ordering contigs of harkess

This python script groups the current contigs by on which new contig they map and selects alle the current contigs with a coverage above 90% on the new contigs
and a coverage above 1% of the new contig on the current contigs and the overall best scoring current contig on coverage. If the contigs that score above the 90% 
and 1% are maximum 100 contigs removed from the best scoring current contig the contigs get saved as a group. The script reordered the current contigs so that 
if a new contig overlapped multiple current contigs but was interrupted by a current contig the current contig that interrupted the alignment would be parsed 
after the last current contig the new contig mapped to. 

Author: Evy Kuenen
Date: 19-12-2023

"""
from collections import defaultdict
from Bio import SeqIO


def read_coordinate_file(coordinaten_bestand):
    """
    Read coord file and extract contig information.

    Args:
        coordinaten_bestand (str): Path to the coordinate information file.

    Returns:
        tuple: A tuple containing mapping_file, contig_order, mpdi_orientation_mapping, mpdi_contig_count.
    """
    mapping_file = {}
    contig_order = []
    mpdi_orientation_mapping = {}
    mpdi_contig_count = defaultdict(lambda: defaultdict(int))

    skipline = 5
    with open(coordinaten_bestand, "r") as file:
        for _ in range(skipline):
            next(file)

        for line in file:
            data = []
            fields = line.strip().split("\t")
            fields_final = [field.split(" | ") for field in fields]
            flat_fields_final = sum(fields_final, [])
            formatted_fields = [item.strip() for sublist in flat_fields_final for item in sublist.split()]
            contig_name = formatted_fields[12]
            orientation = "-" if int(formatted_fields[2]) > int(formatted_fields[3]) else "+"

            if contig_name not in mapping_file:
                mapping_file[contig_name] = []
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
            data.append(orientation)
            mapping_file[contig_name].append(data)
            data.append(contig_name)
            contig_order.append(data)

            mpdi = formatted_fields[11]
            mpdi_contig_count[mpdi][contig_name] += 1
            if mpdi not in mpdi_orientation_mapping:
                mpdi_orientation_mapping[mpdi] = []

            mpdi_orientation_mapping[mpdi].append((orientation, float(formatted_fields[10]), formatted_fields[12]))

    return mapping_file, contig_order, mpdi_orientation_mapping, mpdi_contig_count


def find_most_mapped(mpdi_contig_count):
    """
    Find the most mapped contigs for each MPDI.

    Args:
        mpdi_contig_count (defaultdict): Contig count per MPDI.

    Returns:
        dict: Most mapped contigs for each MPDI.
    """
    most_mapped = {}
    for mpdi, contig_counts in mpdi_contig_count.items():
        most_mapped_contig = max(contig_counts, key=contig_counts.get)
        most_mapped[mpdi] = []
        most_mapped[mpdi].append(most_mapped_contig)
    return most_mapped


def determine_orientation(mpdi_orientation_mapping, most_mapped):
    """
    Determine the orientation for each MPDI based on most mapped contigs.

    Args:
        mpdi_orientation_mapping (dict): Mapping of MPDI to its orientations.
        most_mapped (dict): Most mapped contigs for each MPDI.

    Returns:
        dict: Updated mapping of MPDI to its orientations.
    """
    for mpdi, orientations in mpdi_orientation_mapping.items():
        if len(orientations) >= 2:
            orientation_counts = defaultdict(int)
            total_coverage = defaultdict(float)

            if mpdi in most_mapped:
                contigs_map = most_mapped[mpdi]

                for orientation, coverage, mpdi_value in orientations:
                    if mpdi_value in contigs_map:
                        orientation_counts[orientation] += 1
                        total_coverage[orientation] += float(coverage)

                total_coverage_minus = total_coverage.get('-', 0.0)
                total_coverage_plus = total_coverage.get('+', 0.0)

                if total_coverage_minus > total_coverage_plus:
                    mpdi_orientation_mapping[mpdi] = '-'
                else:
                    mpdi_orientation_mapping[mpdi] = '+'
            else:
                mpdi_orientation_mapping[mpdi] = orientations[0][0]

    return mpdi_orientation_mapping


def filter_mapping(mapping_file, mpdi_orientation_mapping):
    """
    Filter the mapping file based on coverage criteria.

    Args:
        mapping_file (dict): Mapping of contig names to their data.
        mpdi_orientation_mapping (dict): Mapping of MPDI to its orientations.

    Returns:
        dict: Filtered mapping file.
    """
    filtered_mapping_file = {}
    for contig_name, data_list in mapping_file.items():
        mpdi_groups = {}
        for item in data_list:
            mpdi_name = item[11]
            if mpdi_name not in mpdi_groups:
                mpdi_groups[mpdi_name] = []
            mpdi_groups[mpdi_name].append(item)

        for mpdi_name, mpdi_data in mpdi_orientation_mapping.items():
            highest_coverage_item = max(mpdi_data, key=lambda x: float(x[10]))
            if float(highest_coverage_item[9]) >= 90 and float(highest_coverage_item[10]) >= 1.0:
                if contig_name not in filtered_mapping_file:
                    filtered_mapping_file[contig_name] = []
                filtered_mapping_file[contig_name].append(highest_coverage_item)

    return filtered_mapping_file


def get_close_contigs(mapping_file, filtered_mapping_file, distance_threshold=50):
    """
    Get contigs that are close to each other in the mapping file.

    Args:
        mapping_file (dict): Mapping of contig names to their data.
        filtered_mapping_file (dict): Filtered mapping file.
        distance_threshold (int): Maximum distance threshold for contigs to be considered close.

    Returns:
        dict: New mapping file with close contigs.
    """
    new_mapping_file = {}
    for scaffold_name in filtered_mapping_file:
        values_for_scaffold = filtered_mapping_file[scaffold_name]
        if scaffold_name not in new_mapping_file:
            new_mapping_file[scaffold_name] = []
        for values in values_for_scaffold:
            for item in mapping_file[scaffold_name]:
                if abs(float(item[11][4:]) - float(values[11][4:])) <= distance_threshold:
                    new_mapping_file[scaffold_name].append(item)
    return new_mapping_file


def filter_on_coverage(new_mapping_file):
    """
    Filter contigs based on coverage criteria.

    Args:
        new_mapping_file (dict): New mapping file with close contigs.

    Returns:
        dict: Filtered mapping result.
    """
    filtered_mapping_result = {}
    for scaffold_name, values_for_scaffold in new_mapping_file.items():
        mpdi_groups = {}
        for value in values_for_scaffold:
            mpdi_name = value[11]
            if mpdi_name not in mpdi_groups:
                mpdi_groups[mpdi_name] = []
            mpdi_groups[mpdi_name].append(value)

        for mpdi_name, mpdi_data in mpdi_groups.items():
            highest_coverage_item = max(mpdi_data, key=lambda x: float(x[10]))
            coverage_on_reference = highest_coverage_item[10]
            mpdi_data_with_highest_ref_coverage = [v for v in mpdi_data if v[10] == coverage_on_reference]
            highest_coverage_item_on_ref = max(mpdi_data_with_highest_ref_coverage, key=lambda x: float(x[9]))
            if float(highest_coverage_item_on_ref[10]) >= 1.0:
                if scaffold_name not in filtered_mapping_result:
                    filtered_mapping_result[scaffold_name] = []
                filtered_mapping_result[scaffold_name].append(highest_coverage_item_on_ref)

    for scaffold_name in filtered_mapping_result:
        filtered_mapping_result[scaffold_name] = sorted(filtered_mapping_result[scaffold_name], key=lambda x: x[11])

    return filtered_mapping_result


def get_mpdi_order(filtered_mapping_result):
    """
    Get the MPDI order for the final mapping result.

    Args:
        filtered_mapping_result (dict): Filtered mapping result.

    Returns:
        list: List of MPDIs in the desired order.
    """
    getelde_mpd_set = set()
    all_mpdis = []
    for i in range(1000001, 1042066):
        mpdi = f"MPDI{str(i).zfill(8)}.1"
        all_mpdis.append(mpdi)

    max_gap = 50
    for i in range(1000001, 1042066):
        mpdi = f"MPDI{str(i).zfill(8)}.1"
        for scaffold, mpdi_list in filtered_mapping_result.items():
            if mpdi_list and mpdi_list[0][11] == mpdi:
                scaffold_values = [item[11] for item in mpdi_list]
                getelde_mpd_set.add(mpdi)
                count_minus = sum(1 for item in mpdi_list if item[12] == '-')
                count_plus = sum(1 for item in mpdi_list if item[12] == '+')
                index_item1 = all_mpdis.index(scaffold_values[0])

                for mpdi in scaffold_values:
                    if all_mpdis[index_item1] != mpdi:
                        index_item1 += 1
                        while index_item1 < len(all_mpdis) and int(mpdi[4:12]) - int(all_mpdis[index_item1][4:12]) > max_gap:
                            index_item1 += 1

                if count_minus > count_plus:
                    scaffold_values.reverse()

                all_mpdis[index_item1:index_item1] = scaffold_values
    all_mpdis = list(dict.fromkeys(all_mpdis))
    return all_mpdis


def get_unique_mpdis(all_mpdis):
    """
    Get unique MPDIs from the list.

    Args:
        all_mpdis (list): List of MPDIs.

    Returns:
        list: List of unique MPDIs.
    """
    new_all_mpdis = []
    for mpdi in all_mpdis:
        if mpdi not in new_all_mpdis:
            new_all_mpdis.append(mpdi)
    return new_all_mpdis


def save_fasta_sequences(new_all_mpdis, mpdi_orientation_mapping, input_fasta_bestand, output_fasta_bestand):
    """
    Save FASTA sequences based on MPDI order and orientations.

    Args:
        new_all_mpdis (list): List of MPDIs in the desired order.
        mpdi_orientation_mapping (dict): Mapping of MPDI to its orientations.
        input_fasta_bestand (str): Path to the input FASTA file.
        output_fasta_bestand (str): Path to the output FASTA file.
    """
    fasta_sequences = []
    header_to_numeric = {mpdi: float(mpdi[4:]) for mpdi in new_all_mpdis}

    for record in SeqIO.parse(input_fasta_bestand, "fasta"):
        header = record.id
        if header in header_to_numeric:
            if header in mpdi_orientation_mapping:
                orientation = mpdi_orientation_mapping[header]
            else:
                orientation = '+'
            if orientation == '-':
                record.seq = record.seq.reverse_complement()
            fasta_sequences.append(record)

    fasta_sequences.sort(key=lambda record: new_all_mpdis.index(record.id))
    SeqIO.write(fasta_sequences, output_fasta_bestand, "fasta")

    print(f"FASTA-sequenties zijn opgeslagen in {output_fasta_bestand}")


def main():
    coordinaten_bestand = "../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq.coord"
    mapping_file, contig_order, mpdi_orientation_mapping, mpdi_contig_count = read_coordinate_file(coordinaten_bestand)
    most_mapped = find_most_mapped(mpdi_contig_count)
    mpdi_orientation_mapping = determine_orientation(mpdi_orientation_mapping, most_mapped)
    filtered_mapping_file = filter_mapping(mapping_file, mpdi_orientation_mapping)
    new_mapping_file = get_close_contigs(mapping_file, filtered_mapping_file, distance_threshold=1000)
    filtered_mapping_result = filter_on_coverage(new_mapping_file)
    all_mpdis = get_mpdi_order(filtered_mapping_result)
    new_all_mpdis = get_unique_mpdis(all_mpdis)
    input_fasta_bestand = "../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta"
    output_fasta_bestand = "../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/juiste_volgorde_MPDIs_100verplaatsingen.fasta"
    save_fasta_sequences(new_all_mpdis, mpdi_orientation_mapping, input_fasta_bestand, output_fasta_bestand)


main()
