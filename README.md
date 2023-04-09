# Tumor only RAN-seq somatic mutation identification in canine

## **Requirements**

The following package/software is required to run the Tumor only RNA-seq somatic mutation identification pipeline

```
  Python 3.7
  Pandas >= 1.3
  Numpy >= 1.24
  Scikit-learn >= 1.0
  SAMtools 1.9
  Java >= SE8
  Annovar (version 2017 Jul16)
  GATK/3.8-1
  picard/2.21.6
  R Data.table 1.14
```

## **Usage**

The package use a Unix/Linux shell script (CMT-002_somatic_mutation_pipeline.sh) for each sample.
After merging all of the sample results, then use pipeline_ml_mutation_filtering.sh to get the final results.

1. Downloading the package (Linux/Unix platform is required)
2. Modify the shell script (CMT-002_somatic_mutation_pipeline.sh) by specifying the actual path to the whole package directory, the reference directory, the input BAM file mapped with STAR 2-pass alignments, the sample name, and the final pipeline results output

### Steps to run the package

1. Downloading the package (Linux/Unix platform is required)
2. Modify the shell scripts CMT-002_somatic_mutation_pipeline.sh for each sample by specifying the actual path to the whole package directory, the reference directory, the input BAM file mapped with STAR 2-pass alignments, the sample name, and the final pipeline results output, as below.

```
package_location=''                   # path to the tumor-only somatic mutation identification package
bsample=''                            # sample name you want to put
base_folder=''                        # the parent folder of the final results for each sample
bam_file_folder=''                    # path to the folder of the bam file aligned with STAR
ref=''                                # path to the folder of canfam3
annovar_index=''                      # path to the annovar package for gene annotation
```

3. Concat all of the final results (file names ending with \*final_sample_somatic_sum.txt) into a big table
4. Run pipeline_ml_mutation_filtering.sh by specifying the actual path to the whole package directory, the table where you merge all the results from the step3, the bio_tumor_label, the final pipeline results without machine learning prediction and the final pipeline results with machine learning prediction , as below.

```
package_location=''                   # path to the tumor-only somatic mutation identification package
new_add_data=''                       # path to the table where you concat the final results for each sample
bio_tumor_label=''                    # put the tumor type first and the project name sep with '_' ex: OSA_PRJNA000001
pipeline_out_file_name=''             # output file from the pipeline without machine learning predictions
ml_out_file_name=''                   # output file from the pipeline with machine learning predictions
```

4.  Required package/software loading: the shell script contain lines for loading required package/software that are unique to the UGA (Sapelo2) platform. If your platform uses a different approach of loading package/software, you will need to make corresponding changes in the script for lines below.

```
module load GATK/3.8-1-Java-1.8.0_144
module load picard/2.21.6-Java-11
module load Java
module load SAMtools/1.9-GCC-8.3.0
module load R/4.0.0-foss-2019b (need to install data.table and tidyverse library)
module load Miniconda3/4.9.2 (for python3)
module load Perl/5.26.1-GCCcore-6.4.0
```

5.  Once the shell script execution is finished, the result file, "Sample_name.final_sample_somatic_sum.txt", will appear in the directory you specified, containing the following columns, see also example_final_somatic_mutation_identification.txt.

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

To run the pipeline, you will need two extra files (human-dog protein sequence alignment files), but the file size is too large to put in GitHub.
Contact us if you need those files.

Kun-Lin Ho <kh31516@uga.edu> Shaying Zhao <szhao@uga.edu>
