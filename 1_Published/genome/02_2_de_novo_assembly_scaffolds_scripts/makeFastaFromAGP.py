"""Script for making fasta from agp

This script makes from an agp file a fasta file.
Author: Evy Kuenen
Date: 8-1-2024

"""
def read_agp_file(agp_file_path):
    """
    Read an AGP file and extract the scaffold order.

    Parameters:
    - agp_file_path: The path to the AGP file.

    Returns:
    - A list containing the scaffold order.
    """
    order_fasta = []
    with open(agp_file_path, "r") as agp:
        for line in agp:
            fields = line.strip().split("\t")
            order_fasta.append(fields[5])
    return order_fasta

def read_fasta(fasta_file):
    """
    Read a FASTA file and return a dictionary of sequences.

    Parameters:
    - fasta_file: The path to the FASTA file.

    Returns:
    - A dictionary where headers are keys and sequences are values.
    """
    sequences = {}
    current_header = None
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_header = line[1:]
                sequences[current_header] = ''
            else:
                sequences[current_header] += line
    return sequences

def extract_sequences(fasta_file, headers_list):
    """
    Extract sequences from a FASTA file based on a list of headers.

    Parameters:
    - fasta_file: The path to the FASTA file.
    - headers_list: A list of headers to extract.

    Returns:
    - A dictionary containing extracted sequences.
    """
    fasta_sequences = read_fasta(fasta_file)
    extracted_sequences = {header: fasta_sequences[header] for header in headers_list}
    return extracted_sequences

def write_fasta(output_file, sequences_dict):
    """
    Write sequences to a FASTA file.

    Parameters:
    - output_file: The path to the output FASTA file.
    - sequences_dict: A dictionary of sequences.

    Writes the sequences to the specified output file.
    """
    with open(output_file, 'w') as file:
        for header, sequence in sequences_dict.items():
            file.write(f'>{header}\n{sequence}\n')

def main():
    agp_file_path = "../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/makeAGP/medusa.reindexed.agp"
    fasta_file_path = "../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta"
    output_file_path = "../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/makeAGP/medusaContigs.fa"

    order_fasta = read_agp_file(agp_file_path)
    extracted_sequences = extract_sequences(fasta_file_path, order_fasta)
    write_fasta(output_file_path, extracted_sequences)

main()
