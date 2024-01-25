"""Script for finding genes in mitochondrien

This script finds all genes that were found with geseq in the mitochodrien from the gff file that was outputted.

Author: Evy Kuenen
Date: 21-12-2024

"""


def parse_gff_file(gff_file):
    """
    Parse the GFF3 file and extract the gene information from the mitochondrien.

    Parameters:
        gff_file (str): The path to the GFF3 file.

    Returns:
        Tuple[int, int, List[str], List[str]]: A tuple containing the total gene count,
            fragment gene count, list of complete gene names, and list of fragment gene names.
    """
    gene_count = 0
    fragment_gene_count = 0
    gene_names = []
    genes_fragment = []
    genes_complete = []

    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith("mitochondrial_asparargusV2	blatN	gene") or line.startswith("mitochondrial_asparargusV2	blatX	gene"):
                gene_count += 1
                gene_name = line.split('\t')[8].split(';')[0].split('=')[1]
                gene_names.append(gene_name)
                if 'fragment' in gene_name:
                    fragment_gene_count += 1
                    genes_fragment.append(gene_name)
                if 'fragment' not in gene_name:
                    genes_complete.append(gene_name)

    return gene_count, fragment_gene_count, genes_complete, genes_fragment

def print_gene_counts(gene_count, fragment_gene_count, genes_complete, genes_fragment):
    """
    Print the counts and names of the genes.

    Parameters:
        gene_count (int): Total count of genes.
        fragment_gene_count (int): Count of fragment genes.
        genes_complete (List[str]): List of complete gene names.
        genes_fragment (List[str]): List of fragment gene names.
    """
    print(f"Aantal gevonden genen: {gene_count}")
    print(f"Aantal gevonden fragment genen: {fragment_gene_count}")
    print("Complete genen:")
    print(genes_complete)
    print("--------")
    print("Fragment genen:")
    print(genes_fragment)

def main():
    gff_file = "/home/student/projects/asp-pan-genome-evy/data/mitochondrien/GeseqResults/job-results-20231024143612/mitoAsparagus_mitochondrial_asparargusV2_GFF3.gff3"
    gene_count, fragment_gene_count, genes_complete, genes_fragment = parse_gff_file(gff_file)
    print_gene_counts(gene_count, fragment_gene_count, genes_complete, genes_fragment)

main()
