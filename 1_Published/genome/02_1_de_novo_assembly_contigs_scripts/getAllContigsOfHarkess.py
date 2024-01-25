"""Script for obtaining fastas from AGP files

This script gets all the fasta sequences of the current contigs from the AGP file that has been obtained from NCBI

Author: Evy Kuenen
Date: 29-11-2023

"""

from Bio import SeqIO


def get_contigs_from_agp(agp_file, reference_genome_file):
    """Extract contig sequences from an AGP file and a reference genome file.

    Args:
        agp_file (str): Path to the AGP file.
        reference_genome_file (str): Path to the reference genome file.

    Returns:
        dict: A dictionary containing scaffold ID as keys and concatenated contig sequences as values.
    """
    scaffold_sequences = {}
    reference_genome = SeqIO.to_dict(SeqIO.parse(reference_genome_file, "fasta"))

    # Open the AGP file and skip the first 6 lines
    with open(agp_file, "r") as agp:
        for _ in range(6):
            next(agp)

        # Parse the AGP file to extract scaffold sequences
        for line in agp:
            fields = line.strip().split("\t")
            if fields[4] == "W":
                scaffold_id = fields[5]
                contig_id = fields[0]
                contig_start = int(fields[1]) - 1
                contig_end = int(fields[2])
                orientation = fields[8]

                # Extract the contig sequence
                contig_sequence = reference_genome[contig_id].seq

                if scaffold_id not in scaffold_sequences:
                    scaffold_sequences[scaffold_id] = ""

                # Determine the location and orientation of the contig within the scaffold
                if orientation == "+":
                    scaffold_sequences[scaffold_id] += contig_sequence[contig_start:contig_end]
                else:
                    scaffold_sequences[scaffold_id] += contig_sequence[contig_start:contig_end].reverse_complement()
    return scaffold_sequences


def write_to_output(scaffold_sequences, output_file):
    """Write scaffold sequences to an output file.

    Args:
        scaffold_sequences (dict): A dictionary containing scaffold ID as keys and concatenated contig sequences as
        values.
        output_file (str): Path to the output file.

    Returns:
        None
    """
    with open(output_file, "w") as output:
        for scaffold_id, sequence in scaffold_sequences.items():
            output.write(f">{scaffold_id}\n{sequence}\n")


def main():
    agp_file = "../../../data/genome/02_deNovoAssembly/annotation_releases/annotation_releases/4686/100/GCF_001876935.1_Aspof.V1/GCF_001876935.1_Aspof.V1_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP/chr10.comp.agp"
    reference_genome_file = "../../../data/genome/02_deNovoAssembly/annotation_releases/annotation_releases/4686/100/GCF_001876935.1_Aspof.V1/GCF_001876935.1_Aspof.V1_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr10.fna"
    scaffold_sequences = get_contigs_from_agp(agp_file, reference_genome_file)
    output_file = "../../../data/genome/02_deNovoAssembly/superScaffolds/getAllContigsInOldRefseq/chr10AllContigs.fa"
    write_to_output(scaffold_sequences, output_file)


main()
