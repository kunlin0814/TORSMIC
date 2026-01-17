#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script adds gene names to the annovar output from Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt.
"""

import re
import sys

import pandas as pd


def extractAnnovarMutProtein(mut_info):
    try:
        Ensembl_gene = re.findall(r"(ENSCAFG)(\d+):", mut_info)
        Ensembl_trans = re.findall(r"(ENSCAFT)(\d+):", mut_info)

        # Regex for annovar
        total_protein_mut = re.findall(r"p.([A-Z0-9a-z_.*-]*),", mut_info)
        total_ensembl_gene = ["".join(i) for i in Ensembl_gene]
        Total_protein_change = ["".join(i) for i in total_protein_mut]
        total_trans = ["".join(i) for i in Ensembl_trans]
        diff = abs(len(Total_protein_change) != len(total_trans))

        # Add "No_Info_Provided" if annovar ensemble transcripts and consequences have different lengths
        while diff != 0:
            if len(Total_protein_change) > len(total_trans):
                total_trans += ["No_Info_Provided"]
                total_ensembl_gene += ["No_Info_Provided"]
            elif len(Total_protein_change) < len(total_trans):
                Total_protein_change += ["No_Info_Provided"]
            diff = abs(len(Total_protein_change) != len(total_trans))

        # Combine into a single string
        Ensembl_gene = ",".join(total_ensembl_gene)
        Total_protein_change = ",".join(Total_protein_change)
        total_trans = ",".join(total_trans)

        final_return = pd.Series([Ensembl_gene, total_trans, Total_protein_change])

        return final_return
    except:
        print("An error occurred in this sample.")


def append_gene_name(ensembl_gene_string, name_pair_dict):
    ensembl_gene_list = ensembl_gene_string.split(",")
    gene_name_list = []
    for each_ensembl in ensembl_gene_list:
        if each_ensembl in name_pair_dict:
            gene_name_list.append(name_pair_dict[each_ensembl])
        else:
            gene_name_list.append("Unknown")
    ensembl_gene = ensembl_gene_string
    gene_name = ",".join(gene_name_list)
    final_return = pd.Series([ensembl_gene, gene_name])
    return final_return


def get_columns_num(data_file):
    # Delimiter
    data_file_delimiter = "\t"

    # The maximum column count a line in the file could have
    largest_column_count = 0

    # Loop through the data lines
    with open(data_file, "r") as temp_f:
        # Read the lines
        for l in temp_f:
            # Count the column count for the current line
            column_count = len(l.split(data_file_delimiter))

            # Update the largest column count if necessary
            largest_column_count = max(largest_column_count, column_count)

    # Generate column names (0, 1, 2, ..., largest_column_count - 1)
    column_names = [i for i in range(largest_column_count)]
    return column_names


# Get command line arguments
file = sys.argv[1]
name_pair = sys.argv[2]
final_out = sys.argv[3]

# Determine column names
column_names = get_columns_num(file)

# Read annovar info from file
annovar_info = pd.read_csv(file, sep="\t", header=None, names=column_names)

# Read name pair info from file
name_pair = pd.read_csv(name_pair, sep="\t", header=None)
name_pair.columns = ["Ensembl_gene", "Gene_name"]
name_pair_dict = dict(zip(name_pair.Ensembl_gene, name_pair.Gene_name))

# Extract Ensembl gene, Ensembl transcripts, and total protein change
annovar_info.loc[:, ["Ensembl_gene", "Ensembl_transcripts", "Total_protein_change"]] = (
    annovar_info[2].apply(extractAnnovarMutProtein).to_numpy()
)

# Append gene names
annovar_info.loc[:, ["Ensembl_gene", "Gene_name"]] = (
    annovar_info["Ensembl_gene"]
    .apply(append_gene_name, name_pair_dict=name_pair_dict)
    .to_numpy()
)

# Drop unnecessary columns
annovar_info = annovar_info.drop(
    ["Ensembl_gene", "Ensembl_transcripts", "Total_protein_change"], axis=1
)

# Save the final output
annovar_info.to_csv(final_out, sep="\t", index=False, header=None)
