"""Script for extracting contigs of agp file

This script extracts the contigs of an agp file and puts them in a fasta file.
Author: Evy Kuenen
Date: 8-1-2024

"""
from Bio import SeqIO

def extract_contig_sequences(agp_file, fasta_file, output_file):
    """
    Extract contig sequences based on the information provided in the AGP file.

    Parameters:
    - agp_file (str): Path to the AGP file containing contig placement information.
    - fasta_file (str): Path to the FASTA file containing contig sequences.
    - output_file (str): Path to the output file where extracted contig sequences will be written.

    Returns:
    - None
    """
    # Lees AGP-bestand en sla de relevante informatie op
    agp_data = []
    with open(agp_file, 'r') as agp:
        for line in agp:
            if not line.startswith('#'):
                agp_data.append(line.strip().split('\t'))
    # Maak een dictionary om contig-sequenties op te slaan
    contig_sequences = {}
    with open(fasta_file, "r") as fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            contig_sequences[record.id] = record

    # Schrijf de contig-sequenties naar het opgegeven bestand in de volgorde van AGP
    with open(output_file, 'w') as output:
        for agp_entry in agp_data:
            if agp_entry[4] != "U":
                contig_id = agp_entry[5]
                print(agp_entry)
                orientation = agp_entry[8]
                print(contig_id)
                output.write(f'>{contig_id}\n')
                if orientation == '-':
                    output.write(f'{contig_sequences[contig_id].seq.reverse_complement()}\n')
                else:
                    output.write(f'{contig_sequences[contig_id].seq}\n')

def main():
    agp_file_path = '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenWithUnplaced.reindexed.agp'
    fasta_file_path = '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta'
    output_file_path = '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/contigsPythonOrder.fa'

    extract_contig_sequences(agp_file_path, fasta_file_path, output_file_path)

main()