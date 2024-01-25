"""Script for making agp from a coord file

This script makes from an agp file from a coord file with the right chromosome it matches on
Author: Evy Kuenen
Date: 8-1-2024

"""
from collections import defaultdict
import csv


def get_coord_info():
    """
    Read the coordinate information from a file and data in rows and calculate orientation.

    Returns:
    - A list of lists containing coordinate information.
    """
    coordinfo = []
    skip_lines = 5
    with open("../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/100Verplaatsingen/ContigsOrderedFilteredToRefseqV2WithGroupedCentromeres.coord", "r") as coordfile:
        for _ in range(skip_lines):
            next(coordfile) 
            
        for line in coordfile:
            data = []
            fields = line.strip().split("\t")
            fields_final = [field.split(" | ") for field in fields]
            flat_fields_final = sum(fields_final, [])
            formatted_fields = [item.strip() for sublist in flat_fields_final for item in sublist.split()]
            orientation = "-" if int(formatted_fields[2]) > int(formatted_fields[3]) else "+"  # get orientation
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
            data.append(orientation)
            coordinfo.append(data)
    return coordinfo


def get_best_contig(coordinfo):
    """
    Get the best contig information per mpdi based on COVQ and %IDY.

    Parameters:
    - coordinfo: A list of coordinate information.

    Returns:
    - A list containing the best contig information.
    """
    # get best contig per mpdi
    consolidated_data = defaultdict(list)

    for row in coordinfo:
        tag = row[11]  # Laatste kolom is de TAGS-kolom (mpdi)
        covq = float(100/int(row[8])*int(row[5]))  # COVQ-kolom
        idy = float(row[6])  # %IDY-kolom
        consolidated_data[tag].append((covq, idy, row))

    best_contigs = []

    # Loop door de geconsolideerde gegevens en kies de beste contig per contig hierdoor misschien niet alle MPDI's covered
    for tag, data_list in consolidated_data.items():
        # Sorteer op COVQ in omgekeerde volgorde (hoogste eerst), en dan op %IDY in omgekeerde volgorde
        sorted_data = sorted(data_list, key=lambda x: (x[0], x[1]), reverse=True)
        best_contig = sorted_data[0][2]  # Pak de rij met de hoogste COVQ en %IDY
        best_contigs.append(best_contig)
    newBest_contigs = []
    for item in best_contigs:
        if float(1000001.1) <= float(item[11][5:]) <= float(1005032.1):
            item.append('chr1')
            newBest_contigs.append(item)
        if float(1005033.1) <= float(item[11][5:]) <= float(1008443.1):
            item.append('chr2')
            newBest_contigs.append(item)
        if float(1039209.1) <= float(item[11][5:]) <= float(1042065.1):
            item.append('chr10')
            newBest_contigs.append(item)
        if float(1008444.1) <= float(item[11][5:]) <= float(1012930.1):
            item.append('chr3')
            newBest_contigs.append(item)
        if float(1012931.1) <= float(item[11][5:]) <= float(1018162.1):
            item.append('chr4')
            newBest_contigs.append(item)
        if float(1018163.1) <= float(item[11][5:]) <= float(1023181.1):
            item.append('chr5')
            newBest_contigs.append(item)
        if float(1023182.1) <= float(item[11][5:]) <= float(1025973.1):
            item.append('chr6')
            newBest_contigs.append(item)
        if float(1025974.1) <= float(item[11][5:]) <= float(1031595.1):
            item.append('chr7')
            newBest_contigs.append(item)
        if float(1031596.1) <= float(item[11][5:]) <= float(1036550.1):
            item.append('chr8')
            newBest_contigs.append(item)
        if float(1036551.1) <= float(item[11][5:]) <= float(1039208.1):
            item.append('chr9')
            newBest_contigs.append(item)
    return newBest_contigs


def get_chromosome(newBest_contigs):
    """
    Get unique contigs based on coverage.

    Parameters:
    - newBest_contigs: A list of best contig information.

    Returns:
    - A list containing the chromosome that matches the best with the contig.
    """
    # get contig dict contig {mpdi info}
    best_coverage = {}
    for item in newBest_contigs:
        contig = item[12]
        if contig in best_coverage:
            best_coverage[contig].append(item)
        else:
            best_coverage[contig] = [item]

    # get coverage on chromosome in bp
    for contig, item in best_coverage.items():
        sum_cov = {}
        for data in item:
            chromosome = data[14]
            coverage = int(data[5])
            if chromosome in sum_cov:
                sum_cov[chromosome] += coverage
            else:
                sum_cov[chromosome] = coverage
        
        for data in item:
            chromosome = data[14]
            data[5] = str(sum_cov[chromosome])

    # get best contig with chromosome it has the most match on
    best_contigs = []

    for contig, item in best_coverage.items():
        sorted_data = sorted(item, key=lambda x: (int(x[5]), int(x[4])), reverse=True)
        best_contig = sorted_data[0]
        best_contigs.append(best_contig)
    return best_contigs


def get_mpdi_order():
    """
    Get the order of mpdis from a file.

    Returns:
    - A list of mpdis in order.
    """
    mpdis_in_order = []
    with open('../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/juiste_volgorde_MPDIs_100verplaatsingen.fasta', "r") as mpdi_order:
        for line in mpdi_order:
            if line.startswith('>'):
                line2 = line.replace(">", "").strip()
                mpdis_in_order.append(line2)
    return mpdis_in_order


def get_ordened_filtered_coord(mpdis_in_order, best_contigs):
    """
    Filter and order coordinate data based on the specified order and chromosome information.

    Parameters:
    - mpdis_in_order: List of mpdis in the desired order.
    - best_contigs: List of unique contig information.
    - mpdi_on_chromosome: Dictionary mapping mpdis to chromosomes.

    Returns:
    - List of ordered and filtered coordinate data.
    """
    ordened_filtered_coord = []
    for mpdi in mpdis_in_order:
        for item in best_contigs:
            if mpdi == item[11]:
                ordened_filtered_coord.append(item)
    return ordened_filtered_coord


def make_agp(ordened_filtered_coord):
    """
    Generate AGP (A Golden Path) data based on ordered and filtered coordinate information.

    Parameters:
    - ordened_filtered_coord: List of ordered and filtered coordinate data.

    Returns:
    - List of AGP data.
    """
    agp = []
    for item in ordened_filtered_coord:
        data = []
        data.append(item[14])
        data.append(1)
        data.append(item[8])
        data.append(1)
        data.append("W")
        data.append(item[12])
        data.append(1)
        data.append(item[8])
        data.append(item[13])
        data2 = []
        data2.append(item[14])
        data2.append(1)
        data2.append(item[8])
        data2.append(1)
        data2.append("U")
        data2.append(100)
        data2.append('scaffold')
        data2.append('yes')
        data2.append("align_genus")
        agp.append(data)
        agp.append(data2)
    return agp


def write_file(agp):
    """
    Write AGP data to a file.

    Parameters:
    - agp: List of AGP data.
    """
    with open("../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/100Verplaatsingen/100verplaatsingen.agp", "w") as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        writer.writerows(agp)    


def main():
    coordinfo = get_coord_info()
    newBest_contigs = get_best_contig(coordinfo)
    best_contigs = get_chromosome(newBest_contigs)
    mpdis_in_order = get_mpdi_order()
    ordened_filtered_coord = get_ordened_filtered_coord(mpdis_in_order, best_contigs)
    agp = make_agp(ordened_filtered_coord)
    write_file(agp)


main()
