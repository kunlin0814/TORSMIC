import pandas as pd

base_folder=config['base']
package_location=config['in']['package_location']
bsample=config['in']['bsample']                    # sample name
bio_tumor=config['in']['bio_tumor']          # put the tumor type first and then project name,  separate with '_' 
annovar_index = config['in']['annovar_index'] ## the directory that contains canine annovar annotation files  
db_snp=config['in']['dbsnp']
# '/work/szlab/Lab_shared_PanCancer/source/DbSNP_canFam3_version151-DogSD_Broad_March2022.vcf' 
scripts_location = config['in']['scripts_location']                    ## the directory where the scripts of the package locates
data_source_location = config['in']['data_source_location']          ## the directory where the datasource of the package locates
gatk_vcf_dir = config['in']['gatk_vcf_dir'] ## the directory that contains bam align with STAR 2-pass  


somatic_output_dir = config['out']['somatic_output_dir']  ## the directory where you want to put your result for each sample. It can be separated from the bam file directory  


### required file 

add_gene_script = f"{package_location}/scripts/Update_Add_GeneName.py"
gene_ensembl_pair = f"{package_location}/data_source/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt"
cds_file = f"{package_location}/data_source/UniqueCainineCdsInterval.interval_list"


rule all:
    input:
        torsmic_classify_results = f"{somatic_output_dir}/{bsample}_cds_dbsnp_mut_summary.txt"


##  A function that performs annovar annotation and then append the gene name
rule filter_with_dbsnp:
    input:
        dbsnp_file = db_snp,
        gatk_filtered_vcf = f"{gatk_vcf_dir}/{bsample}_rg_added_sorted_dedupped_split.realigned.bam.filter.vcf"  
    output:
        dbsnp_filtered = f"{somatic_output_dir}/DB_SNP_filtering_{bsample}.filter.vcf",
        germline_filtered_info = f"{somatic_output_dir}/Germline_filtered_{bsample}.filter.vcf"
    conda:
        config['envs']
    shell: 
        '''
        java -Xmx{resources.mem} -cp {scripts_location} DbSNP_filtering
        {input.dbsnp_file} \
        {input.gatk_filtered_vcf} \
        {output.dbsnp_filtered} \
        {output.germline_filtered_info}
        '''

rule limit_cds:
    input:
        dbsnp_filtered = rules.filter_with_dbsnp.output.dbsnp_filtered
    output:
        cds_dbsnp_filtered = f"{somatic_output_dir}/{bsample}_cds_dbsnp.filter.vcf"
    conda:
        config['envs']
    params:
        scripts_location = scripts_location
    shell:
        '''
        python {params.scripts_location}/Limit_vcf_to_CDS.py \
        {input.dbsnp_filtered} \
        {input.cds_file} \
        {output.cds_dbsnp_filtered}
        '''

rule annovar_annotation:
    input:
        vcf_file = somatic_output_dir+ '/{file_name}.filter.vcf'
        #rules.limit_cds.cds_dbsnp_filtered
    output:
        vcf_pass = temp(somatic_output_dir+"/{file_name}.filter.vcf-PASS"),
        vcf_pass_avinput = temp(somatic_output_dir+"/{file_name}.filter.vcf-PASS-avinput"),
        annovar_out = temp(somatic_output_dir+"/{file_name}.filter.vcf-PASS-avinput.exonic_variant_function"),
        add_genename = temp(somatic_output_dir+"/{file_name}.filter.vcf-PASS-avinput.exonic_variant_function_WithGeneName"),
        vcf_pass_with_sample = somatic_output_dir+"/{file_name}_filter_vcf_withsample.vcf",
        annovar_with_sample = somatic_output_dir+"/{file_name}_annovar_gene_withsample"
    params:
        annovar_index=config['in']['annovar_index'],
        add_gene_script = add_gene_script,
        gene_ensembl_pair = gene_ensembl_pair,
        bsample_name = bsample
    conda:
        config['envs']
    shell:
        '''
        awk '$7 == "PASS" {{print $0}}' {input.vcf_file} > {output.vcf_pass}
        perl {params.annovar_index}/convert2annovar.pl -format vcf4old {output.vcf_pass} > {output.vcf_pass_avinput}

        # annovar annotat
        perl {params.annovar_index}/annotate_variation.pl --buildver canFam3 {output.vcf_pass_avinput} {params.annovar_index}

        python {params.add_gene_script} \
        {output.annovar_out} \
        {params.gene_ensembl_pair} \
        {output.add_genename}

        ## append the sample names to the gatk and annovar gene_names
        cat {output.vcf_pass} | awk -v awkvar="{params.bsample_name}" -F "\t" 'BEGIN {{FS=OFS ="\t"}} {{print $0 OFS awkvar}}' \
        > {output.vcf_pass_with_sample}

        cat {output.add_genename} | \
        awk -v awkvar="{params.bsample_name}" -F "\t" 'BEGIN {{FS=OFS ="\t"}} {{print $0 OFS awkvar}}' \
        > {output.annovar_with_sample}

        '''


rule classify_mutation:
    input:
        vcf_pass_with_sample = somatic_output_dir+"/{file_name}_filtered_vcf_withsample.vcf",
        annovar_with_sample = somatic_output_dir+"/{file_name}_annovar_gene_withsample"
    output:
       mut_classify_results = somatic_output_dir+"/{file_name}_mut_summary.txt"
       #f"{somatic_output_dir}/{bsample}_cds_dbsnp_mutation_summary.txt"
    conda:
        config['envs']
    params:
        package_location = package_location,
        bio_tumor_symbol = bio_tumor
    shell:
        '''
        python "${scripts_location}/Sapelo2_extract_somatic_germline.py" \
        {somatic_output_dir}/{input.vcf_pass_with_sample} \
        {somatic_output_dir}/{input.annovar_with_sample} \
        {bsample} \
        {somatic_output_dir}/{output.mut_classify_results} \
        {params.package_location} \
        {params.bio_tumor_symbol}
        '''

