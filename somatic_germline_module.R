library(data.table)
library(tidyverse)
library(gridExtra)
library(HDInterval)

append_element <- function(each_table){
  
  long_format <- paste(unique(each_table$Gene_mut_info), collapse = ",")
  sum <- list('Gene_mut_source'= long_format)
  return (sum)
}

count_chrom_mut_info_project <- function(each_chrom_mut_info_project, two_project_cut = each_bioproject_cut, one_project_cut=0.4){
  bioproject_num <- length(unique(each_chrom_mut_info_project$Bioproject))
  each_chrom_mut_info_project_ratio <- unique(each_chrom_mut_info_project[,.(Each_bioproject_ratio,Bioproject)])$Each_bioproject_ratio
  gt2_cut_off <- each_chrom_mut_info_project_ratio[each_chrom_mut_info_project_ratio >= two_project_cut]
  gt1_cut_off <- each_chrom_mut_info_project_ratio[each_chrom_mut_info_project_ratio >= one_project_cut]
  num_project_gt2_cut_off <- length(gt2_cut_off)
  num_project_gt1_cut_off <- length(gt1_cut_off)
  gt2_ratio <- paste(gt2_cut_off,collapse = ";")
  
  final_sum <- list(In_number_bioproject = bioproject_num,
                    num_project_with_gt2_cut_off = num_project_gt2_cut_off,
                    num_project_with_gt1_cut_off = num_project_gt1_cut_off,
                    gt2_ratio =gt2_ratio)
  return (final_sum)
}




filtering_5_step <- function(final_results){
  final_results[,Total_tumor_reads:=Ref_reads+Alt_reads]
  final_results <- final_results[VAF>=0.05 ,]
  final_results <- final_results[Total_tumor_reads >= 10,]
  final_results <- final_results[!(Alt_reads <=5 & VAF <0.15), ]
  final_results <- final_results[!(Total_tumor_reads < 20 & VAF< 0.2), ]
  final_results[,Total_tumor_reads:=NULL]
  
  return (final_results)
}

### divided VAF plot into different group
classify_sample_num <- function(sample_num){
  final_result <- " "
  if (sample_num < 5){
    final_result <- "<5"
  }else if (5<=sample_num & sample_num<10){
    final_result <- "5-9"
  }else if (sample_num>=10){
    final_result <- ">=10"
  }
  return (final_result)
}

check_proportion <- function(group_table,vaf_cut, midd_range,upper_range){
  #group_table <- setDT(group_table)
  each_group_target <-  group_table
  VAF_values <- group_table$VAF
  total_sample <- length(VAF_values)
  hetero_ratio <- (sum((VAF_values >= midd_range[1] & VAF_values <= midd_range[2]))/total_sample)
  homo_ratio <- (sum((VAF_values >= upper_range[1] & VAF_values <= upper_range[2]))/total_sample)
  meet_criteria <- sum((VAF_values >= midd_range[1] & VAF_values <= midd_range[2]) | (VAF_values >= upper_range[1] & VAF_values <= upper_range[2]))
  
  ratio <- meet_criteria/total_sample
  return (list('Homo_Hetero_Ratio'= ratio,
               'Homo_ratio'= homo_ratio,
               'Hetero_ratio' = hetero_ratio))
}
