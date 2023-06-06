#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Oct 29 14:36:15 2021

The script will take each annovar output and DBSNP_filtering_gatk as inputs and extract somatic mutation that overlaps with canine pan-cancer 
and c-bioportal data (after translated to human position)
Notice:
1. The final output might have fewer records than the original gatk vcf data because the script will filter out the records that don't have VAF info.
2. We use the genomic mutation (chrom+pos+ref+alt) identified in pan-cancer to identify pan-cancer records.
3. We use transcript mutations (because only one transcript) info to identify mutations that can be found in c-bio and cosmic database.

It will create one output:
1. Final_sample_sum_out (the df contains the annovar info that has the mutation found in pan-cancer (including synonymous mutations), c-bioprotal, cosmic, and remained)
"""

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from natsort import index_natsorted, natsort_keygen, natsorted, order_by_index

# Create the argument parser
parser = argparse.ArgumentParser(description="Script to extract somatic mutations")
# Add the command-line arguments
parser.add_argument("gatk_vcf", type=str, help="Path to the GATK VCF file")
parser.add_argument("annovar_gene_file", type=str, help="Path to the Annovar gene file")
parser.add_argument("sample_name", type=str, help="Name of the sample")
parser.add_argument("final_sample_sum_out", type=str, help="Path to the output file")
parser.add_argument("package_location", type=str, help="Path to the package location")
parser.add_argument(
    "bio_project",
    type=str,
    help="Tumor and bioproject information, separate with '_', ex: MT_PRJNA00001",
)
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit(1)
# Parse the command-line arguments
args = parser.parse_args()

# Input data
gatk_vcf = args.gatk_vcf
annovar_gene_file = args.annovar_gene_file
sample_name = args.sample_name
# Output data
final_sample_sum_out = args.final_sample_sum_out
package_location = args.package_location
bio_project = args.bio_project
# Load the module
module_loc = os.path.join(package_location, "scripts")
sys.path.append(module_loc)
from somatic_germline_module import *

package_location = Path(package_location)
translate_to = "human"

pan_cancer_annovar_file = (
    package_location
    / "data_source"
    / "Ge2_Pass_QC_Pan_Cancer_Final_Mutect_annovar_include_syn_mutation_summary.txt"
)
c_bioportal_file = (
    package_location / "data_source" / "all_studies_c-bio_portal_somatic_mutation.txt"
)
cosmic_file = (
    package_location / "data_source" / "GRCh37_V95_Cosmic_somatic_mutation.txt"
)
c_bio_translate_file = (
    package_location
    / "data_source"
    / "c-bio_Human_GR37_103_canine_3.199_sequenceAlignment.txt"
)
cosmic_translate_file = (
    package_location
    / "data_source"
    / "COSMIC_V95_Human_GR37_93_canine_3.199_sequenceAlignment.txt"
)
retro_gene_file = package_location / "data_source" / "retro_gene_list.txt"
c_biohuman_dog_transcript = (
    package_location
    / "data_source"
    / "c_bioportal_Human_GR37_103_dog_transcript_3.199.txt"
)
cosm_human_dog_transcript = (
    package_location
    / "data_source"
    / "COSMIC_Human_GR37_V95_93_dog_transcript_3.199.txt"
)
## process annovar out and extract all of the annovar information
## the function will auto-extract the target data
target_annovar_info = processAnnovar(annovar_gene_file, retro_gene_file, sample_name)

### process gatk output
target_gatk = process_gatk_output(gatk_vcf)

### merge GATK and annovar information
merge_gatk_annovar = target_annovar_info.merge(
    target_gatk, on=["Line", "Chrom", "Sample_name"], how="left"
)
merge_gatk_annovar = merge_gatk_annovar.loc[
    pd.notna(merge_gatk_annovar["Ref_reads"])
    & pd.notna(merge_gatk_annovar["Alt_reads"])
    & (~merge_gatk_annovar["Ref_reads"].astype(str).str.contains("No info provided"))
    & (~merge_gatk_annovar["Alt_reads"].astype(str).str.contains("No info provided"))
]

merge_gatk_annovar["VAF"] = merge_gatk_annovar["Alt_reads"].astype(float) / (
    merge_gatk_annovar["Ref_reads"].astype(float)
    + merge_gatk_annovar["Alt_reads"].astype(float)
)


# Ref and Alt use GATK format not annovar format, so I extract Ref_y, Alt_y
target_column = [
    "Line",
    "Consequence",
    "Gene_name",
    "Chrom",
    "Start",
    "End",
    "Sample_name",
    "Ensembl_gene",
    "Ensembl_transcripts",
    "Total_protein_change",
    "Gene_mut_info",
    "Transcript_mut_info",
    "Ref_y",
    "Alt_y",
    "VAF_info",
    "Ref_reads",
    "Alt_reads",
    "VAF",
]

target_merge_gatk_annovar = merge_gatk_annovar[target_column]
target_merge_gatk_annovar_columns = target_merge_gatk_annovar.columns.values

target_merge_gatk_annovar = merge_gatk_annovar[target_column]
# change column names: Ref_y,Alt_y to Ref, Alt
target_merge_gatk_annovar.rename(columns={"Ref_y": "Ref", "Alt_y": "Alt"}, inplace=True)
target_merge_gatk_annovar["Chrom_mut_info"] = (
    target_merge_gatk_annovar["Chrom"]
    + "_"
    + target_merge_gatk_annovar["Start"].astype(str)
    + "_"
    + target_merge_gatk_annovar["Ref"]
    + "_"
    + target_merge_gatk_annovar["Alt"]
)

## filtering with canine pan-cancer, change to use genomic_location to idenfiy somatic mutation rather than using gene names mutation or transcripts mutation
## process pan-cancer data
pan_cancer_data = pd.read_csv(pan_cancer_annovar_file, sep="\t")
pan_cancer_source = list(
    set(
        list(
            pan_cancer_data["Chrom"]
            + "_"
            + pan_cancer_data["Pos"].astype(str)
            + "_"
            + pan_cancer_data["Ref"]
            + "_"
            + pan_cancer_data["Alt"]
        )
    )
)
pan_cancer_pass = target_merge_gatk_annovar[
    target_merge_gatk_annovar["Chrom_mut_info"].isin(pan_cancer_source)
]

######### Process human somatic data
## process c-bio files
c_bio_pass = extract_human_somatic(
    c_bioportal_file,
    c_bio_translate_file,
    c_biohuman_dog_transcript,
    target_merge_gatk_annovar,
    translate_to,
)

cosm_pass = extract_human_somatic(
    cosmic_file,
    cosmic_translate_file,
    cosm_human_dog_transcript,
    target_merge_gatk_annovar,
    translate_to,
)


cosm_pass_uniq = cosm_pass.merge(c_bio_pass, how="outer", indicator=True).loc[
    lambda x: x["_merge"] == "left_only"
]
c_bio_pass_uniq = cosm_pass.merge(c_bio_pass, how="outer", indicator=True).loc[
    lambda x: x["_merge"] == "right_only"
]

if not pan_cancer_pass.empty:
    pan_cancer_pass.loc[:, "Source"] = "Pan-cancer"

if not c_bio_pass.empty:
    c_bio_pass.loc[:, "Source"] = "C-bio"

if not cosm_pass_uniq.empty:
    cosm_pass_uniq.loc[:, "Source"] = "Cosmic"

## The majority of Cosmic is overlapped with C-bio, so if I found the same, remove it
cosm_pass_uniq = cosm_pass_uniq.drop(columns="_merge")

if pan_cancer_pass.empty and c_bio_pass.empty and cosm_pass_uniq.empty:
    final_panancer_cbio_cosmic = pd.DataFrame([], columns=["Line", "Chrom_mut_info"])
else:
    final_panancer_cbio_cosmic = pd.DataFrame()
    final_panancer_cbio_cosmic = pd.concat(
        [pan_cancer_pass, c_bio_pass, cosm_pass_uniq], ignore_index=True
    ).drop_duplicates(subset=["Line", "Chrom_mut_info"])

# Keep remaining Annovar files not in pan-cancer and c-bio for future VAF examination
passed_info = final_panancer_cbio_cosmic["Chrom_mut_info"].unique().tolist()

## remained mutations are those not found in c-bio, cosmic, and pan-cancer
remained_df = (
    target_merge_gatk_annovar[
        ~target_merge_gatk_annovar["Chrom_mut_info"].isin(passed_info)
    ]
    .drop(columns=["C_bio_Human_counterpart", "Cosm_Human_counterpart"])
    .assign(Source="Remained")
)

# Combine remaining data with passed data
total_final_out = final_panancer_cbio_cosmic.append(remained_df).drop_duplicates()

# Remove 'VAF_info' column to reduce the file size and sort the DataFrame
total_final_out = (
    total_final_out.drop("VAF_info", axis=1)
    .sort_values(by="Line", key=natsort_keygen())
    .drop(columns="Line")
)

# Assign 'bio_tumor' value to the DataFrame
total_final_out["Bioproject"] = bio_project

total_final_out.to_csv(final_sample_sum_out, sep="\t", index=False)
