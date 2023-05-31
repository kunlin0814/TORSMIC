# Tumor only RAN-seq somatic mutation identification in canine (TORSMIC)

## **Requirements**

The following package/software is required to run TORSMIC pipeline

```
  Python 3.7
  natsort >= 8.3
  Pandas >= 1.3
  Numpy >= 1.24
  Scikit-learn >= 1.0
  Java >= SE8
  Annovar (version after 2017 Jul16)
  GATK/3.8-1
  picard/2.21.6
  R Data.table 1.14
```

## **General Usage**

The package use a Unix/Linux shell script.

1. Use CMT-002_somatic_mutation_pipeline.sh to run each sample.
2. Merge the results from all of the samples
3. Use pipeline_ml_mutation_filtering.sh to get the final mutation classification.

### Steps to run the package

1. Downloading the package (Linux/Unix platform is required)
2. Modify the shell scripts CMT-002_somatic_mutation_pipeline.sh for each sample and specify the actual path for the following directories.

The overall structure looks like below:

base_foldery/
base_foldery/star_align_bam_dir/each_sample
base_foldery/somatic_output_folder/each_sample

```
package_location=''                   # path to the tumor-only somatic mutation identification package
bsample=''                            # sample name you want to put
base_folder=''                        # the parent directory of the final results for each sample
bam_file_folder=''                    # path to the folder of the bam file aligned with STAR
bio_tumor=''                          # put the tumor type first and the project name sep with '_' ex: OSA_PRJNA000001
reference=''                          # path to the directory of canfam3
annovar_index=''                      # path to the annovar package for gene annotation
```

Notice:
We will add new data to our previous results to identify somatic mutations. To get the best results, we need to set the tumor type the same as our previous naming conventions. Our previous analyzed results include the following tumor types:
a. Bladder tumor -BLA
b. Glioma - GLM
c. Hemangiosarcoma - HSA
d. Mammary Tumor - MT
e. Oral melanoma - OM
f. Osteosarcoma - OSA
f. Prostate cancer - PRO
If your sample is derived from the same tumor type shown above, please use the same acronym, but if your dataset doesn't have the tumor type shown above, then you can name your tumor type as you like. Ex: OSA_PRJNA000001

3. Concat all of the final results derived from CMT-002_somatic_mutation_pipeline.sh (file names ending with \*final_sample_somatic_sum.txt) into a big table
4. Run pipeline_ml_mutation_filtering.sh by specifyingthe actual path for the following directories.

```
package_location=''                   # path to the tumor-only somatic mutation identification package
new_add_data=''                       # path to the table where you concat the final results for each sample
pipeline_out_file_name=''             # output file from the pipeline without machine learning predictions
ml_out_file_name=''                   # output file from the pipeline with machine learning predictions
```

5.  Required package/software loading: the shell script contain lines for loading required package/software that are unique to the UGA (Sapelo2) platform. If your platform uses a different approach of loading package/software, you will need to make corresponding changes in the script for lines below.

```
module load GATK/3.8-1-Java-1.8.0_144
module load picard/2.21.6-Java-11
module load Java
module load R/4.0.0-foss-2019b            (need to install data.table and tidyverse library)
module load Miniconda3/4.9.2              (a custom conda envirnoment python3)
module load Perl/5.26.1-GCCcore-6.4.0
```

6.  Once the shell script execution is finished, the result file, "\*ml_filtering.txt", will appear in the directory you specified, which contains the following columns, see also example_final_somatic_mutation_identification.txt in the demo directory.

```
[1] Sample_name
[2] Bioproject
[3] VAF   (Variants allele frequency)
[4] Source    (mutations that can be matched to human somatic database or not )
[5] Consequence
[6] Gene_name
[7] Chrom
[8] Pos   (mutations start position)
[9] End   (mutations end position)
[10] Ensembl_gene
[11] Ensembl_transcripts
[12] Total_protein_change
[13] Ref    (Reference base)
[14] Alt    (Alternative base)
[15] Ref_reads    (No. of reference reads detected with GATK)
[16] Alt_reads    (No. of alt reads detected with GATK)
[17] Model_prediction (final mutations results, ex: germline or somatic or WT)

```

## **Other required file**

To run the pipeline, you will also need four extra files
[1] two human-dog protein sequence alignment files  
[2] annovar package to get dog gene annotation,
[3] a db_snp files to filter known germline variants

Because the files size above are too large to put in GitHub, contact me if you need those files.

Kun-Lin Ho <abc730814@gmail.com>
