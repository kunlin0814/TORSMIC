#!/bin/bash
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=4           # Number of CPU cores per task (4)
#SBATCH --mem=70G                   # Job memory limit (10 GB)
#SBATCH --time=100:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=CMT-002_extract_somatic.%j.out    # Standard output log
#SBATCH --error=CMT-002_extract_somatic.%j.err     # Standard error log

package_location='/home/kh31516/kh31516_Lab_Share_script/IdentifiyNeoAntigene/TumorOnlySomatic'
base_folder=/scratch/kh31516/UGA
bsample=CMT-002                     # sample name
bio_tumor='MT_PRJNA000001'          # put the tumor type first and then project name,  separate with '_' 


somatic_output_folder=${base_folder}/somatic_results/${bsample} ## the directory where you want to put your result for each sample. It can be separated from the bam file directory  
scripts_location=${package_location}/scripts                    ## the directory where the scripts of the package locates
data_source_location=${package_location}/data_source            ## the directory where the datasource of the package locates


### For GATK calling on STAR 2-pass alignment 
bam_file_folder=${base_folder}/results/${bsample}/STAR ## the directory that contains bam align with STAR 2-pass  
reference='/work/szlab/dog_resouces/source'  ## the directory that contains canFam3 reference sequence 

### required file 
cds_file=${data_source_location}'/UniqueCainineCdsInterval.interval_list'
annovar_index='/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf' ## the directory that contains canine annovar annotation files  
db_snp='/work/szlab/Lab_shared_PanCancer/source/DbSNP_canFam3_version151-DogSD_Broad_March2022.vcf' 


mkdir -p $somatic_output_folder

## a custom conda environment for running Python 
ml Miniconda3/4.9.2
source activate /home/kh31516/myenv/py38

##  A function that performs annovar annotation and then append the gene name
annovar_gene_annotation(){
vcf_file=$1
final_gatk_out=$2
final_output_name=$3
final_output_folder=$4

module load Perl/5.26.1-GCCcore-6.4.0

awk '$7 == "PASS" {print $0}' ${vcf_file} > ${vcf_file}-PASS

perl $annovar_index/convert2annovar.pl -format vcf4old ${vcf_file}-PASS > ${vcf_file}-PASS-avinput

perl $annovar_index/annotate_variation.pl --buildver canFam3 ${vcf_file}-PASS-avinput $annovar_index

# Use ensemble ID to append gene names
python ${package_location}/Update_Add_GeneName.py \
${vcf_file}-PASS-avinput.exonic_variant_function \
${data_source_location}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt \
${vcf_file}-PASS-avinput.exonic_variant_function_WithGeneName

## append the sample name to the gatk mutation calling result and annovar annotation result
cat ${vcf_file}-PASS | awk -v awkvar="${bsample}" -F "\t" 'BEGIN {FS=OFS ="\t"} {print $0 OFS awkvar}' \
> ${final_output_folder}/${final_gatk_out}

cat ${vcf_file}-PASS-avinput.exonic_variant_function_WithGeneName | awk -v awkvar="${bsample}" -F "\t" 'BEGIN {FS=OFS ="\t"} {print $0 OFS awkvar}' \
> ${final_output_folder}/${final_output_name}

}
####### Picard #######
# Sort the bam file
picard(){
ml picard/2.21.6-Java-11
bam_file_folder=$1

java -Xms4g -Xmx16g -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${bam_file_folder}/${bsample}.bam \
O=${bam_file_folder}/${bsample}-rg_added_sorted.bam SO=coordinate RGID=id RGLB=library \
RGPL=platform RGPU=machine RGSM=${bsample}
# Mark duplicates
java -Xms4g -Xmx16g -jar $EBROOTPICARD/picard.jar  MarkDuplicates I=${bam_file_folder}/${bsample}-rg_added_sorted.bam \
O=${bam_file_folder}/${bsample}_rg_added_sorted_dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
M=${bam_file_folder}/${bsample}-output.metrics

}

####### GATK mutation calling and  Annovar annotation #######
gatk_annovar(){
ml GATK/3.8-1-Java-1.8.0_144
bam_file_folder=$1
#Split'N'Trim and reassign mapping qualities
java -Xms16g -Xmx48g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R $reference/canFam3.fa \
--maxReadsInMemory 1000000 \
-I ${bam_file_folder}/${bsample}_rg_added_sorted_dedupped.bam \
-o ${bam_file_folder}/${bsample}-rg_added_sorted_dedupped_split.bam \
-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
# Generating interval file for sort.bam
java -Xms16g -Xmx64g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa \
-I ${bam_file_folder}/${bsample}-rg_added_sorted_dedupped_split.bam \
-o ${bam_file_folder}/${bsample}-rg_added_sorted_dedupped_split.bam.intervals \
--allow_potentially_misencoded_quality_scores
# Realign
java -Xms16g -Xmx48g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa \
--maxReadsInMemory 1000000 \
-I ${bam_file_folder}/${bsample}-rg_added_sorted_dedupped_split.bam \
-targetIntervals ${bam_file_folder}/${bsample}-rg_added_sorted_dedupped_split.bam.intervals \
-o ${bam_file_folder}/${bsample}-rg_added_sorted_dedupped_split.realigned.bam \
--allow_potentially_misencoded_quality_scores
# Variant calling
java -Xms16g -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference/canFam3.fa \
--maxTotalReadsInMemory 10000000 \
-I ${bam_file_folder}/${bsample}-rg_added_sorted_dedupped_split.realigned.bam \
-dontUseSoftClippedBases -stand_call_conf 20.0 \
-o ${bam_file_folder}/${bsample}-rg_added_sorted_dedupped_split.realigned.bam.vcf
# Variant filtering
java -Xms16g -Xmx48g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R $reference/canFam3.fa \
-V ${bam_file_folder}/${bsample}-rg_added_sorted_dedupped_split.realigned.bam.vcf \
-window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${bam_file_folder}/${bsample}-rg_added_sorted_dedupped_split.realigned.bam.filter.vcf

### Annovar annotation to the GATK mutation alling 
annovar_gene_annotation ${bam_file_folder}/${bsample}-rg_added_sorted_dedupped_split.realigned.bam.filter.vcf \
${bsample}-gatk_file_withSamplename \
${bsample}-annovar_WithSampleName \
${somatic_output_folder}

}


####################################################
## Steps for picard sorting and mark duplicates 
picard ${bam_file_folder}
## Steps for GATK mutation calling and then annotate gene with anovar
gatk_annovar ${bam_file_folder}


####################################################
## Steps for Somatic mutation identification

## required input
# gatk_vcf = sys.argv[1]
# annovar_gene_file = sys.argv[2]
# sample_name = sys.argv[3]
# final_sample_sum_out = sys.argv[4]
# data_source_folder = sys.argv[5]
# bio_tumor = sys.argv[6]

## get all of the information before the pipeline
python ${scripts_location}/Sapelo2_extract_somatic_germline.py \
${somatic_output_folder}/${bsample}-gatk_file_withSamplename \
${somatic_output_folder}/${bsample}-annovar_WithSampleName \
${bsample} \
${somatic_output_folder}/Before_pipeline_${bsample}_sample_somatic_sum.txt \
${package_location} \
${bio_tumor} 

### Load Java module
ml Java
# remove germline mutations found in the databse
java -Xmx32g -cp  ${scripts_location}/ DbSNP_filtering \
$db_snp \
${bam_file_folder}/${bsample}-rg_added_sorted_dedupped_split.realigned.bam.filter.vcf-PASS \
${bam_file_folder}/DB_SNP_filtering_${bsample}-rg_added_sorted_dedupped_split.realigned.bam.filter.vcf \
${bam_file_folder}/Germline_filtered_${bsample}-rg_added_sorted_dedupped_split.realigned.bam.filter.vcf


### limit to CDS region
python ${scripts_location}/Limit_vcf_to_CDS.py \
${bam_file_folder}/DB_SNP_filtering_${bsample}-rg_added_sorted_dedupped_split.realigned.bam.filter.vcf \
$cds_file \
${bam_file_folder}/CDS_DB_SNP_filtering_${bsample}-rg_added_sorted_dedupped_split.realigned.bam.filter.vcf

## annovar annotation to CDS variants
annovar_gene_annotation ${bam_file_folder}/CDS_DB_SNP_filtering_${bsample}-rg_added_sorted_dedupped_split.realigned.bam.filter.vcf \
CDS_DB_SNP_filtering_${bsample}-gatk_file_withSamplename \
CDS_DB_SNP_filtering_${bsample}-annovar_WithSampleName \
${somatic_output_folder}

### input data
# vcf_file = sys.argv[1]
# annovar_gene_file = sys.argv[2]
# SRR_run = sys.argv[3]
## output data
# final_sample_sum_out = sys.argv[4]

### classify the variants c-bio, cosmic, and remained
python ${scripts_location}/Sapelo2_extract_somatic_germline.py \
${somatic_output_folder}/CDS_DB_SNP_filtering_${bsample}-gatk_file_withSamplename \
${somatic_output_folder}/CDS_DB_SNP_filtering_${bsample}-annovar_WithSampleName \
${bsample} \
${somatic_output_folder}/${bsample}_final_sample_somatic_sum.txt \
${package_location} \
${bio_tumor} 
