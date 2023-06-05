# Tumor-only RAN-seq Somatic Mutation Identification in Canine (TORSMIC)

## Overview

TORSMIC is a pipeline designed for identifying somatic mutations in canine tumor samples using RNA-seq data. It utilizes a combination of tools and scripts to perform variant calling, annotation, and classification. The pipeline aims to provide accurate and reliable somatic mutation identification in a tumor-only setting without matching normal samples.

## Requirements

To run the TORSMIC pipeline, ensure that the following requirements are met:

- Python 3.7
- natsort >= 8.3
- Pandas >= 1.3
- Numpy >= 1.24
- Scikit-learn >= 1.0
- Java >= SE8
- Annovar (version after 2017 Jul16)
- GATK/3.8-1
- picard/2.21.6
- R Data.table 1.14

## General Usage

The TORSMIC pipeline utilizes a Unix/Linux shell script and involves several steps, as outlined below:

1. Run `CMT-002_somatic_mutation_pipeline.sh` for each sample to perform variant calling and somatic mutation identification.
2. Merge the results obtained from individual samples.
3. Utilize `pipeline_ml_mutation_filtering.sh` to apply machine learning-based mutation filtering and classification to obtain the final mutation calls.

### Steps to Run the Package

1. Download the package (Linux/Unix platform is required).
2. Modify the shell script `CMT-002_somatic_mutation_pipeline.sh` for each sample and specify the actual path for the following directories.

The overall directory structure should resemble the following:

```
base_folder/
base_folder/star_align_bam_dir/each_sample/each_sample.bam
base_folder/somatic_output_folder/each_sample/each_sample_final_sample_somatic_sum.txt
```

```bash
package_location=''                   # Path to the tumor-only somatic mutation identification package
bsample=''                            # Sample name you want to use
base_folder=''                        # Parent directory of the final results for each sample
bam_file_folder=''                    # Path to the folder containing the bam files aligned with STAR
bio_tumor=''                          # Specify the tumor type first and the project name separated by '_' (e.g., OSA_PRJNA000001), see also Note:
reference=''                          # Path to the directory of canfam3 reference sequence
annovar_index=''                      # Path to the Annovar package for gene annotation
```

Note:
To achieve the best results, we need to set the tumor type the same as our previous naming conventions. Our previously analyzed results include the following tumor types:
a. Bladder tumor - BLA
b. Glioma - GLM
c. Hemangiosarcoma - HSA
d. Mammary Tumor - MT
e. Oral melanoma - OM
f. Osteosarcoma - OSA
g. Prostate cancer - PRO

If your sample is derived from the same tumor type as listed above, please use the same acronym. However, if your dataset doesn't match the tumor types mentioned above, you can name your tumor type as desired (e.g., OSA_PRJNA000001).

3. Concatenate all the final results derived from `CMT-002_somatic_mutation_pipeline.sh` (results file ending with `*final_sample_somatic_sum.txt`) into a single table.
4. Run `pipeline_ml_mutation_filtering.sh` by specifying the actual path for the following directories.

```bash
package_location=''                   # Path to the tumor-only somatic mutation identification package
new_add_data=''                       # Path to the table where you concatenated the final results for each sample
pipeline_out_file_name=''             # Output file from the pipeline without machine learning predictions
ml_out_file_name=''                   # Output file from the pipeline with machine learning predictions
```

5. Required package/software loading: The shell script contains lines for loading required packages/software that are unique to the UGA (Sapelo2) platform. If your platform uses a different approach to load packages/software, make the corresponding changes in the script for the following lines.

```bash
module load GATK/3.8-1-Java-1.8.0_144
module load picard/2.21.6-Java-11
module load Java
module load R/4.0.0-foss-2019b            # (Requires data.table and tidyverse library installation)
module load Miniconda3/4.9.2              # (Custom Conda environment with Python 3)
module load Perl/5.26.1-GCCcore-6.4.0
```

6. Once the shell script execution is finished, the result file "\*ml_filtering.txt" will appear in the specified directory. The file includes various columns representing sample information, variant details, and model predictions. The columns present in the file are as follows:

```
[1] Sample_name
[2] Bioproject
[3] VAF   (Variant allele frequency)
[4] Source    (Mutations that can be matched to the human somatic database or not)
[5] Consequence
[6] Gene_name
[7] Chrom
[8] Pos   (Mutation start position)
[9] End   (Mutation end position)
[10] Ensembl_gene
[11] Ensembl_transcripts
[12] Total_protein_change
[13] Ref    (Reference base)
[14] Alt    (Alternative base)
[15] Ref_reads    (Number of reference reads detected with GATK)
[16] Alt_reads    (Number of alternative reads detected with GATK)
[17] Model_prediction (Final mutation results, e.g., germline, somatic, or WT)
```

These columns provide valuable information about each mutation detected in the samples, including genomic coordinates, variant allele frequency, gene annotation, consequence, and the final mutation classification based on the machine learning model.

Please refer to the generated "\*ml_filtering.txt" file for detailed information about the mutations identified in the tumor samples.

## Additional Files

In addition to the TORSMIC pipeline, you will need the following files to run the pipeline successfully:

```
[1] Two human-dog protein sequence alignment files
[2] Annovar package for obtaining dog gene annotation
[3] A db_snp file to filter known germline variants
```

Due to their large file size, these files cannot be included in the GitHub repository. Please contact me at the provided email address if you require these files.

Contact Information:
Kun-Lin Ho
Email: abc730814@gmail.com

We will be happy to help and provide any necessary support.
