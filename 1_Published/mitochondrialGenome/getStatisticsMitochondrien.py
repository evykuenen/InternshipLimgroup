"""Script for getting statistics of the new mitochondrial genome

This script gets from a fasta file the GC percentage and the base composition.

Author: Evy Kuenen
Date: 18-12-2023

"""
from Bio import SeqIO

def calculate_gc_percentage(genome_filename):
    """
    Calculate the GC percentage of mitochondrial genome.

    Parameters:
        genome_filename (str): The path to the FASTA file containing the genome sequence.

    Returns:
        float: The GC percentage of the genome sequence.
    """
    with open(genome_filename, 'r') as genome_file:
        genome_sequence = SeqIO.read(genome_file, 'fasta').seq

    num_G = genome_sequence.count('G')
    num_C = genome_sequence.count('C')

    seq_without_unknown = ''.join(base for base in genome_sequence if base.upper() != 'N')
    total_bases = len(seq_without_unknown)

    gc_percentage = ((num_G + num_C) / total_bases) * 100

    return gc_percentage

def get_base_compositions(genome_filename):
    """
    Print the percentage composition of each DNA base in the mitochondrial genome.

    Parameters:
        genome_filename (str): The path to the FASTA file containing the genome sequence.
    """
    with open(genome_filename, 'r') as genome_file:
        genome_sequence = SeqIO.read(genome_file, 'fasta').seq

    num_G = genome_sequence.count('G')
    num_C = genome_sequence.count('C')
    num_A = genome_sequence.count('A')
    num_T = genome_sequence.count('T')    
    seq_without_unknown = ''.join(base for base in genome_sequence if base.upper() != 'N')
    total_bases = len(seq_without_unknown)

    g_percentage = (num_G / total_bases) * 100
    c_percentage = (num_C / total_bases) * 100
    a_percentage = (num_A / total_bases) * 100
    t_percentage = (num_T / total_bases) * 100

    print(f"g: {g_percentage}")
    print(f"c: {c_percentage}")
    print(f"a: {a_percentage}")
    print(f"t: {t_percentage}")

def main():
    """
    Main function to execute the GC percentage calculation and base compositions.

    Prints the GC percentage and base compositions of the mitochondrial genome.
    """
    genome_filename = '/home/student/projects/asp-pan-genome-evy/data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/newRefSeqMito.fa'
    result = calculate_gc_percentage(genome_filename)
    get_base_compositions(genome_filename)
    print(f"GC percentage of the reference genome: {result:.2f}%")

main()
