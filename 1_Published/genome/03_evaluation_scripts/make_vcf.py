"""Script for making a vcf file from a genetic map

This script makes a vcf file from a genetic map
Author: Evy Kuenen
Date: 11-1-2024

"""
import csv


def read_data(file_path, skip_lines=1):
    """
    Read data from a CSV file.
    Returns a list of lists containing the data.
    """
    data = []

    with open(file_path, 'r') as file_handle:
        for _ in range(skip_lines):
            next(file_handle)

        for line in file_handle:
            fields = line.strip().split("\t")
            data.append(fields)
    return data


def convert_to_vcf(data):
    """
    Convert the input data to VCF format.
    Returns a list of lists representing VCF lines.
    """
    vcf = []

    for item in data:
        lines = []
        lines.append('AsparagusV1_' + str(item[1]).zfill(2))
        lines.append(item[3])
        lines.append(item[0])
        lines.append('N')
        lines.append('N')
        lines.append(1)
        lines.append('.')
        lines.append('.')
        lines.append('GT:PL')
        lines.append('optional')
        vcf.append(lines)
    return vcf


def write_vcf(vcf_data, output_file_path):
    """
    Write VCF data to a specified output file.
    """
    with open(output_file_path, "w") as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        writer.writerows(vcf_data)


def main():
    input_file_path = "../../../data/genome/allmapsGetContigOrderCurrentContig/1_onemap/genetische_kaart_K397.map.csv"
    output_file_path = "../../../data/genome/allmapsGetContigOrderCurrentContig/1_onemap/asparagusV2.vcf"
    data = read_data(input_file_path)
    vcf_data = convert_to_vcf(data)
    write_vcf(vcf_data, output_file_path)


main()
