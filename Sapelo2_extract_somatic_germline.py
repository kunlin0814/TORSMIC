#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 14:36:15 2021

The script will take each annovar output and DBSNP_filtering_gatk as inputs and extract somatic mutation that overlaps with canine pan-cancer 
and c-bioportal data (after translated to human position)
Notice:
1. The final output might have fewers records than the original gatk vcf data, because the script will filter out the the records that don't have VAF info
2. Now we use genome mutation(chrom+pos+ref+alt) identify in pan-cancer to identify pan-cancer records, 
but still use transcript mutation(because only one transcript, so gene_aa_mutaiotn is the same) info to idenitfy c-bio

It will create one output
1. Final_sample_sum_out ( the df conatins the annovar info that has the mutation found in pan-cancer (include synonymous mutations),c-bioprotal and remained)

### 1/28 we change pan-cancer annovar file as an argument 
### 3/31/2022
## include Cosmic dataset but follow the same logic as c-bioportal
@author: kun-linho
"""

import os
import re
import sys

import pandas as pd
from natsort import index_natsorted, natsort_keygen, natsorted, order_by_index


def extractVAF(gatk_info):
    if len(gatk_info.split(":")) > 1:
        vaf_info = gatk_info.split(":")[1].split(",")
        ref = int(vaf_info[0])
        alt = int(vaf_info[1])
        # vaf_value = float(alt/(ref+alt))
    else:
        ref = "No info provide"
        alt = "No info provide"

    return pd.Series([ref, alt])


### input data
gatk_vcf = sys.argv[1]
annovar_gene_file = sys.argv[2]
sample_name = sys.argv[3]
pan_cancer_annovar_file = sys.argv[4]
## output data
final_sample_sum_out = sys.argv[5]
data_source_folder = sys.argv[6]
module_loc = data_source_folder
sys.path.append(module_loc)
from somatic_germline_module import *

data_source_folder = Path(data_source_folder)

translate_to = "human"


c_bioportal_file = data_source_folder / "all_studies_c-bio_portal_somatic_mutation.txt"
comsmic_file = data_source_folder / "GRCh37_V95_Cosmic_somatic_mutation.txt"
c_bio_translate_file = (
    data_source_folder / "c-bio_Human_GR37_103_canine_3.199_sequenceAlignment.txt"
)
cosmic_translate_file = (
    data_source_folder / "COSMIC_V95_Human_GR37_93_canine_3.199_sequenceAlignment.txt"
)
retro_gene_file = data_source_folder / "retro_gene_list.txt"
c_biohuman_dog_transcript = (
    data_source_folder / "c_bioportal_Human_GR37_103_dog_transcript_3.199.txt"
)
cosm_human_dog_transcript = (
    data_source_folder / "COSMIC_Human_GR37_V95_93_dog_transcript_3.199.txt"
)
## process annovar out and extract all of the annovar information
retro_gene_list = pd.read_csv(retro_gene_file, sep="\n", header=None)
retro_gene_list = retro_gene_list[0].to_list()
## the function will auto-extract the target data
target_annovar_info = processAnnovar(annovar_gene_file, retro_gene_file, sample_name)


### process gatk output
gatk_data = pd.read_csv(gatk_vcf, sep="\t", header=None)
gatk_data.loc[:, "Line"] = ["line" + str(i + 1) for i in range(0, len(gatk_data))]
gatk_data.loc[:, ["Ref_reads", "Alt_reads"]] = (
    gatk_data[9].astype(str).apply(extractVAF).to_numpy()
)
## only extract the needed info from gatk_data
target_gatk = gatk_data.loc[:, [0, 3, 4, 9, 10, "Ref_reads", "Alt_reads", "Line"]]
target_gatk.columns = [
    "Chrom",
    "Ref",
    "Alt",
    "VAF_info",
    "Sample_name",
    "Ref_reads",
    "Alt_reads",
    "Line",
]
merge_gatk_annovar = target_annovar_info.merge(
    target_gatk, on=["Line", "Chrom", "Sample_name"], how="left"
)
merge_gatk_annovar = merge_gatk_annovar.loc[
    (merge_gatk_annovar.Ref_reads != "No info provide")
    & (merge_gatk_annovar.Alt_reads != "No info provide")
]
merge_gatk_annovar.loc[:, "VAF"] = merge_gatk_annovar["Alt_reads"].astype(float) / (
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

# change column names: Ref_y,Alt_y to Ref, Alt
target_merge_gatk_annovar_columns = np.where(
    target_merge_gatk_annovar_columns == "Ref_y",
    "Ref",
    target_merge_gatk_annovar_columns,
)
target_merge_gatk_annovar_columns = np.where(
    target_merge_gatk_annovar_columns == "Alt_y",
    "Alt",
    target_merge_gatk_annovar_columns,
)
target_merge_gatk_annovar.columns = target_merge_gatk_annovar_columns
target_merge_gatk_annovar.loc[:, "Chrom_mut_info"] = (
    target_merge_gatk_annovar["Chrom"]
    + "_"
    + target_merge_gatk_annovar["Start"].astype(str)
    + "_"
    + target_merge_gatk_annovar["Ref"]
    + "_"
    + target_merge_gatk_annovar["Alt"]
)

## process c-bio files
c_bio = pd.read_csv(c_bioportal_file, sep="\t")
## create a total somatic mutation list idenitfied in c-bio
all_c_bio_list = ",".join(c_bio["Mut_type"]).split(",")

## process cosmic files
cosmic = pd.read_csv(comsmic_file, sep="\t")
## create a total somatic mutation list idenitfied in c-bio
all_cosmic_list = ",".join(cosmic["Mut_type"]).split(",")


## filtering with pan-cancer, change to use genomic_location to idenfiy somatic mutation rather than using gene names mutation or transcripts mutation
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

## filter with c-bioprotal but need to translate to human first
## because alignment only align genes with gene_names, we can't remove pos doesn't have gene name
# gene_name_annovar = target_merge_gatk_annovar[target_merge_gatk_annovar['Gene_name']!='-']
#### processing human_dog_translate_allign_file
c_biotranslate_allign_table = pd.read_csv(c_bio_translate_file, sep="\t")
### fill up the gap with '-' , otherwise, it will have many same position even the alignment is a gap
c_biotranslate_allign_table.loc[
    c_biotranslate_allign_table["QueryAA"] == "-", "QueryIdx"
] = "-"
clean_translate_table = c_biotranslate_allign_table[
    c_biotranslate_allign_table["QueryIdx"] != "-"
]
clean_translate_table = clean_translate_table.to_records()
## a function to create four dict that can used to search ini the future
c_biototal_summary_dict = createDictforHumanDogSearch(clean_translate_table)

### creating human_dog_translate file end ###
target_merge_gatk_annovar["C_bio_Human_counterpart"] = target_merge_gatk_annovar[
    "Gene_mut_info"
].apply(
    identify_species_counterparts,
    human_dog_pos_dict=c_biototal_summary_dict["human_dog_pos_dict"],
    dog_human_pos_dict=c_biototal_summary_dict["dog_human_pos_dict"],
    human_aa_dict=c_biototal_summary_dict["human_aa_dict"],
    dog_aa_dict=c_biototal_summary_dict["dog_aa_dict"],
    translate_to=translate_to,
)

c_bio_human_dog_transcript_info = pd.read_csv(
    c_biohuman_dog_transcript, sep="\t", header=None
)
c_bio_human_dog_transcript_info.columns = [
    "Gene_name",
    "Human_transcripts",
    "Dog_transcripts",
]
c_biotranscript_match_target_annovar = target_merge_gatk_annovar.loc[
    target_merge_gatk_annovar["Ensembl_transcripts"].isin(
        c_bio_human_dog_transcript_info["Dog_transcripts"]
    )
]
c_bio_pass = c_biotranscript_match_target_annovar.loc[
    c_biotranscript_match_target_annovar["C_bio_Human_counterpart"].isin(all_c_bio_list)
]

## remove the C_bio_Human_counterpart column to merge with pan-cancer pass
c_bio_pass = c_bio_pass.drop(columns="C_bio_Human_counterpart")

## filter with COSMIC but need to translate to human first
## because alignment only align genes with gene_names, we can't remove pos doesn't have gene name
# gene_name_annovar = target_merge_gatk_annovar[target_merge_gatk_annovar['Gene_name']!='-']
#### processing human_dog_translate_allign_file
cosm_translate_allign_table = pd.read_csv(cosmic_translate_file, sep="\t")
### fill up the gap with '-' , otherwise, it will have many same position even the alignment is a gap
cosm_translate_allign_table.loc[
    cosm_translate_allign_table["QueryAA"] == "-", "QueryIdx"
] = "-"
cosm_clean_translate_table = cosm_translate_allign_table[
    cosm_translate_allign_table["QueryIdx"] != "-"
]
cosm_clean_translate_table = cosm_clean_translate_table.to_records()
## a function to create four dict that can used to search ini the future
cosm_total_summary_dict = createDictforHumanDogSearch(cosm_clean_translate_table)


### creating human_dog_translate file end ###
target_merge_gatk_annovar["Cosm_Human_counterpart"] = target_merge_gatk_annovar[
    "Gene_mut_info"
].apply(
    identify_species_counterparts,
    human_dog_pos_dict=cosm_total_summary_dict["human_dog_pos_dict"],
    dog_human_pos_dict=cosm_total_summary_dict["dog_human_pos_dict"],
    human_aa_dict=cosm_total_summary_dict["human_aa_dict"],
    dog_aa_dict=cosm_total_summary_dict["dog_aa_dict"],
    translate_to=translate_to,
)

cosm_human_dog_transcript_info = pd.read_csv(
    cosm_human_dog_transcript, sep="\t", header=None
)
cosm_human_dog_transcript_info.columns = [
    "Gene_name",
    "Human_transcripts",
    "Dog_transcripts",
]
cosm_transcript_match_target_annovar = target_merge_gatk_annovar.loc[
    target_merge_gatk_annovar["Ensembl_transcripts"].isin(
        cosm_human_dog_transcript_info["Dog_transcripts"]
    )
]
cosm_pass = cosm_transcript_match_target_annovar.loc[
    cosm_transcript_match_target_annovar["Cosm_Human_counterpart"].isin(all_cosmic_list)
]

## remove the Cosm_Human_counterpart column to merge with pan-cancer pass
cosm_pass = cosm_pass.drop(
    columns=["C_bio_Human_counterpart", "Cosm_Human_counterpart"]
)

cosm_pass_uniq = cosm_pass.merge(c_bio_pass, how="outer", indicator=True).loc[
    lambda x: x["_merge"] == "left_only"
]
c_bio_pass_uniq = cosm_pass.merge(c_bio_pass, how="outer", indicator=True).loc[
    lambda x: x["_merge"] == "right_only"
]

if pan_cancer_pass.empty == False:
    pan_cancer_pass.loc[:, "Source"] = "Pan-cancer"

if c_bio_pass.empty == False:
    c_bio_pass.loc[:, "Source"] = "C-bio"

if cosm_pass_uniq.empty == False:
    cosm_pass_uniq.loc[:, "Source"] = "Cosmic"

## The majority of Cosmic is overlapped with C-bio, so if I found the same, remove it
cosm_pass_uniq = cosm_pass_uniq.drop(columns="_merge")

if (
    (pan_cancer_pass.empty == True)
    and (c_bio_pass.empty == True)
    and (cosm_pass_uniq.empty == True)
):
    final_panancer_cbio_cosmic = pd.DataFrame([], columns=["Line", "Chrom_mut_info"])
else:
    final_panancer_cbio_cosmic = pd.DataFrame()
    final_panancer_cbio_cosmic = final_panancer_cbio_cosmic.append(pan_cancer_pass)
    final_panancer_cbio_cosmic = final_panancer_cbio_cosmic.append(c_bio_pass)
    final_panancer_cbio_cosmic = final_panancer_cbio_cosmic.append(
        cosm_pass_uniq
    ).drop_duplicates()

## keep remained annovar files that are not in pan-cancer and c-bio for the future VAF examine
pass_info = final_panancer_cbio_cosmic.Chrom_mut_info.unique().tolist()
remained_df = target_merge_gatk_annovar[
    ~target_merge_gatk_annovar["Chrom_mut_info"].isin(pass_info)
].drop(columns=["C_bio_Human_counterpart", "Cosm_Human_counterpart"])
remained_df.loc[:, "Source"] = "Remained"

total_final_out = final_panancer_cbio_cosmic.append(remained_df).drop_duplicates()
total_final_out = total_final_out.drop(["VAF_info"], axis=1)
total_final_out = total_final_out.sort_values(by="Line", key=natsort_keygen()).drop(
    columns="Line"
)
total_final_out.to_csv(final_sample_sum_out, sep="\t", index=False)
