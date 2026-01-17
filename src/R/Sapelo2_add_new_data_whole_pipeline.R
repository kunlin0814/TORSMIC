library(data.table)
library(tidyverse)

args    <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

package_location <- as.character(args[1])
new_data_file <- as.character(args[2])
out_file_name <- as.character(args[3])


module_location <- paste(package_location,'scripts','somatic_germline_module.R',sep ="/")
data_source <- paste(package_location,'data_source',sep="/")
setwd(data_source)
source(module_location)

sample_cut_off <- 10
second_sample_cut_off <- 5
propor_cut_off <- 0.5 # if 0.4 <= VAF <=0.6 in >0.5 among variants, then germline 
second_propor_cut_off <- 0.4 # if 0.4 <= VAF <=0.6 in >0.4 , then germline
midd_range <- c(0.4, 0.6)
upper_range <- c(0.9,1)
vaf_cut <- 0.90
## if the VAF of the mutations ratio >0.2, then the mut is germline 
vaf_ratio_cut <- 0.2
## 5/24/22 add new cut off to see if the variant exist in each bioproject and total samples
each_tumor_type_min_sample <- 10
total_sample_cut <- 0.15 # PIK3CA 0.13
gt1_tumor_type_cut <- 0.3 ## BRAF 0.31, PIK3CA is 0.27
one_tumor_type_cut <- 0.45

filtering_5steps <- T
sample_recurrent <- T
breed_filtering <- T
wes_rna_filtering <- T
new_data <- fread(new_data_file)
file <- fread("PassQC_Total_biosample_somatic_germline_sum_03_02_23.txt.gz")
file[Source=='C-bio'| Source=="Cosmic",Source:="Human",]
old_sample <- unique(file$Sample_name)

### GATK sometimes call the mutation with , ex: C >T,G, so we need to separate by "," for further analysis
file <- data.table(separate_rows(file,Alt,sep =","))
file[,Chrom_mut_info:=paste(Chrom,Start,Ref,Alt, sep ="_")]

### check if tumor type is enter
col_miss <- setdiff(colnames(file),colnames(new_data))
# ### add tumor type for those newly added data
# new_data[,(col_miss):=bio_tumor_label]
new_data[Source=='C-bio'| Source=="Cosmic",Source:="Human",]
if (nrow(new_data)!=0){
  file <- rbind(file,new_data)
}else{
  print("No mutations found")
}

## here is you decided whether you want to merge different tumor types
file[,c("Cancer","Project") := tstrsplit(Bioproject,"_",fixed=T)]
file[Cancer =="UB", Cancer:="BLA"]
## merge different tumor type
merge_bio <- file[,unique(Project),keyby = .(Cancer)]
total_disease <- unique(merge_bio$Cancer)
for (i in total_disease){
  total_bioproject <- unique(merge_bio[Cancer==i]$V1)
  if (length(total_bioproject)>1){
    final_label <- paste(c(i,total_bioproject),collapse="_")
    file[Project %in% total_bioproject,Bioproject:=final_label]
  }
}
file[ ,`:=`(Cancer = NULL, Project = NULL)]


## this is the breed-enriched (specific) variants identified in WES
breed_enrich_variants_file <- 'BreedEnrichedMut10Breeds.txt'


### Use MT RNA-WES information to filter
wes_rna_overlap <- fread("Total_WES_overlap_GATK_info_02_09_22.txt")
wes_rna_overlap <- filtering_5_step(wes_rna_overlap)
wes_rna_overlap_somatic <- unique(wes_rna_overlap$Chrom_mut_info)

wes_rna_seq <- 'Total_WES_uniqe_to_pipeline_GATK_info_11_07_22_training_summarize.txt'
wes_rna_seq_training_data <- fread(wes_rna_seq)
germline_info <- unique(wes_rna_seq_training_data[Results=="Germline",]$Chrom_mut_info)
somatic_info <- unique(wes_rna_seq_training_data[Results=="Somatic",]$Chrom_mut_info)
wt_info <- unique(wes_rna_seq_training_data[Results=="WT",]$Chrom_mut_info)


source <- c('Pan-cancer','Human','Remained')

if (filtering_5steps){
  file <- filtering_5_step(file)
}
## calculate the sample_number 
file[,Sample_number := length(unique(Sample_name)),keyby = .(Chrom_mut_info)]
target_proportation_file <- unique(file[,.(Chrom_mut_info,VAF,Sample_name)])
## find the ratio of samples with VAF in 0.4-0.6 and >=0.9 and the ratio of samples with VAF in >= 0.9
## separate into two steps to make the program to faster
file_samples_sum <-  target_proportation_file[,check_proportion(.SD,vaf_cut=vaf_cut,
                                                                midd_range=midd_range,upper_range=upper_range),
                                              by =.(Chrom_mut_info),.SDcols=c('VAF', 'Sample_name')]

target_proportation_file <- NULL
### Merge all the info
file <- merge(file,file_samples_sum,by =  'Chrom_mut_info', all.x = T)

## assign all somatic first and then assign the germline in the future
file[,Status:="Somatic",]
file_samples_sum <- NULL
##### process pan-cancer,
pan_unique_data <- unique(file[Source %in% source[1],
                               .(Chrom_mut_info,Sample_name,VAF,Source,Sample_number,Homo_ratio,Homo_Hetero_Ratio, Status)])

##### condition 1:  >= 10 samples
## condition 1-1: if >= propor_cut_off and samples >= 10, then germline
pan_unique_data[Sample_number>=sample_cut_off & Homo_Hetero_Ratio >=propor_cut_off, Status:="Germline"]

## 12/17 add extra category 5-9
## condition 2: 
## 12/7 if the mutation in at least 5 sample, and the herto_homo ratio >= second_propor_cut_off, then they are germline
pan_unique_data[Sample_number<sample_cut_off &Sample_number>=second_sample_cut_off 
                & Homo_Hetero_Ratio >= second_propor_cut_off,
                Status:="Germline"]

## condition 3:<= 4 samples, 12/17, add extra category <=4 samples
## 08/09/22
## if the mutation has le 4 samples, then if one vaf >0.9, then germline 
pan_unique_data[Sample_number<second_sample_cut_off & 
                  (VAF %between% midd_range),Status:="Germline", ]

###### Combine C-bio and Cosmic as Human and process Human
## Human only needs to included that cut-off for the ratio of VAF > 0.9 > 0.18, no need to use sample distribution

human_unique_data <- unique(file[Source %in% source[2],
                                 .(Chrom_mut_info,Sample_name,VAF,Source,Sample_number,Homo_ratio,Homo_Hetero_Ratio,Status)])
# ## assign all c-bio mutation as Somatic mutation and then later filter out with VAF>0.9 if the ratio > 0.18 for sample >=5
human_unique_data[Sample_number >= second_sample_cut_off & Homo_ratio>=vaf_ratio_cut, Status:='Germline']

## 2-2. if sample is less than 5, then if one vaf > 0.9 , then all of variants are germline
human_unique_data[Sample_number < second_sample_cut_off & Homo_ratio>=vaf_ratio_cut, Status:='Germline']

##### process remained 
remained_unique_data <- unique(file[Source %in% source[3],.(Chrom_mut_info,Sample_name,VAF,Source,Sample_number,Homo_ratio,Homo_Hetero_Ratio,Status)])
unique_remained_training_data <- unique(file[Source %in% source[3],.(Chrom_mut_info,Sample_name,Ref_reads,Alt_reads,Gene_mut_info,VAF,Source,Sample_number,Homo_ratio,Homo_Hetero_Ratio,Status)])

## condition 1:  >= 10 samples ,
## condition 1-1: if Homo_Hetero_Ratio >= propor_cut_off and samples >= 10, then germline
remained_unique_data[Sample_number>=sample_cut_off & Homo_Hetero_Ratio >=propor_cut_off,Status:="Germline"]
## 12/17 add extra category 5-9
## condition 2-1: 
## 12/7 if the mutation has at least 5 smaples and <10 samples and the herto_homo ratio >= second_propor_cut_off, then they are germline
remained_unique_data[Sample_number<sample_cut_off & Sample_number>=second_sample_cut_off
                     & Homo_Hetero_Ratio >= second_propor_cut_off, Status:="Germline",]


## condition 3:<= 4 samples, 12/17, add extra category <=4 samples
## if the mutation has le 4 samples, then check individual VAF in the 
remained_unique_data[Sample_number<second_sample_cut_off &  (VAF %between% midd_range),Status:="Germline", ]

## combine all
### combine Pan-cancer, c-bio and remained.
final_sum <- rbind(pan_unique_data,human_unique_data)
final_sum <- data.table(unique(rbind(final_sum,remained_unique_data)))

## append back to original file
file[,Status:=NULL,]
overlap_col <- intersect(colnames(file), colnames(final_sum))
file <- merge(file, final_sum, by = overlap_col, all.x = T)
## add new Homo_ratio cut-off and C-bio hoping to remove VAF>0.90
file[Homo_ratio>=vaf_ratio_cut , Status:='Germline']
### summarized that mutation recurrent in all sample level and each tumor level
### all sample level 
total_sample_num <- length(unique(file$Sample_name))
file[,Variant_in_all_sample_ratio := Sample_number/total_sample_num]
### each tumor type level
file[,Variants_in_number_tumor_type:=length(unique(Bioproject)), keyby =.(Chrom_mut_info)]

each_tumor_type_samples <- file[,.(Each_tumor_type_total_sample=length(unique(Sample_name))),keyby = .(Bioproject)]
target_tumor_type <- unique(file[,.(Chrom_mut_info,Sample_name,Bioproject)])

target_tumor_type[,Variants_in_each_tumor_type := length(unique(Sample_name)),
                  keyby = .(Bioproject,Chrom_mut_info )]

target_tumor_type <- merge(target_tumor_type, each_tumor_type_samples,
                           by = "Bioproject", all.x = T)

target_tumor_type$Each_tumor_type_ratio <- target_tumor_type$Variants_in_each_tumor_type/target_tumor_type$Each_tumor_type_total_sample
overlap_col <- intersect(colnames(file),colnames(target_tumor_type ))
file <- merge(file,target_tumor_type, by = overlap_col, all.x = T )
before_sample_recurrent_filtering <- copy(file[Status=='Somatic',])

###### sample_recurrent filtering
## add new cut off to see if the variant exist in each tumor type and total Samples
## filter out the variants that exist in many samples
## consider each tumor type (tumor type must contains at least 10 samples
## if not, then we don't use each tumor type cut to filter variants
## in case after QC, some tumor only have 1 or 2 samples)

if (sample_recurrent){
  ## extract variants presented in all samples gt than cut-off (0.15) 
  gt_all <- unique(before_sample_recurrent_filtering[Variant_in_all_sample_ratio>=total_sample_cut]$Chrom_mut_info)
  
  ## look at how many samples have specific variants within each tumor type
  ## if each tumor type contains less than 10 samples or the variants can't reach gt2 tumor type cutoff (0.3), then we don't examine that variants)
  target_tumor_type_chrom_mut <- unique(before_sample_recurrent_filtering[Each_tumor_type_total_sample>=each_tumor_type_min_sample & Each_tumor_type_ratio>=one_tumor_type_cut,]$Chrom_mut_info)
  target_tumor_type <-  before_sample_recurrent_filtering[Chrom_mut_info %in% target_tumor_type_chrom_mut, ]   
  ## summarize the variants existing in how many tumor and if those variants gt than one tumor type cut or two tumor type cut 
  each_tumor_target_mut_sum <-  target_tumor_type[,count_chrom_mut_info_project(.SD), keyby = .(Chrom_mut_info)]
  ### if variants in gt two tumor types, if variants have ratios > gt2 tumor types, the variants are germline
  germline_in_gt2_tumor <- each_tumor_target_mut_sum[In_number_tumor_type>=2 & num_tumor_type_with_gt2_cut_off >=2 ,]$Chrom_mut_info
  ### if variants in one tumor type, then variants can't have ratio gt than 1 tumor type cut off (set 0.4 for each)
  germline_in_1_tumor <- each_tumor_target_mut_sum[num_tumor_type_with_gt1_cut_off >=1,]$Chrom_mut_info
  ### merge variants from total sample cut or each bioproject cut
  overlap <- unique(c(germline_in_gt2_tumor, germline_in_1_tumor,gt_all))
  file[Chrom_mut_info %in% overlap, Status:="Germline",]
}

if (breed_filtering){
  breed_enrich_variants <- fread(breed_enrich_variants_file)
  breed_enrich_variants[,Locus:=gsub(":","_",Locus)][,Mutation:=gsub(">","_",Mutation)][,Chrom_mut_info:=paste(Locus, Mutation,sep ="_")]
  breed_enrich_variants <- unique(breed_enrich_variants$Chrom_mut_info)
  file[Chrom_mut_info %in% breed_enrich_variants, Status:="Germline",]
}

if (wes_rna_filtering){
  file[Chrom_mut_info %in% wt_info, Status:="WT"]
  file[Chrom_mut_info %in% germline_info, Status:="Germline"]
  file[Chrom_mut_info %in% somatic_info, Status:="Somatic"]
  file[Chrom_mut_info %in% wes_rna_overlap_somatic , Status:="Somatic"]
}

file <- file[grepl("Somatic",Status,ignore.case = T),]
file <- unique(file[!Sample_name %in% old_sample,])

fwrite(file, file = out_file_name, 
       sep = "\t",eol = "\n",col.names = T, quote = F)
