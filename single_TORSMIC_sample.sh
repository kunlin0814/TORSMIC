#!/bin/bash
#SBATCH --partition=iob_p
#SBATCH --job-name=test_RNA
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --tasks-per-node=1
#SBATCH --mem=60G
#SBATCH --time=30:00:00
#SBATCH --mail-user=kh31516@uga.edu
#SBATCH --mail-type=END
#SBATCH --output=/scratch/kh31516/snakeTest/scripts/test_RNA.out

base="/scratch/${USER}/snakeTest"
bio_tumor=OSA_PRJNA78827
bsample='SAMN00013513' 

package_location='/work/szlab/Kun-Lin_pipeline/TORSMIC'
gatk_vcf_dir='/scratch/kh31516/snakeTest/results/OSA_PRJNA78827/SAMN00013513/STAR'
dbsnp='/work/szlab/Lab_shared_PanCancer/source/DbSNP_canFam3_version151-DogSD_Broad_March2022.vcf'

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate TORSMIC

mkdir -p ${base}
cd ${base}
mkdir -p logs/ config/ data/ out/ results/

#### test WES germline calling pipeline per case
config=$base/config/torsmic_config.json

python ${package_location}/scripts/TORSMIC_make_snakemake_config.py \
    --base ${base} \
    --out ${config} \
    --bio_tumor ${bio_tumor} \
    --package_location ${package_location} \
    --bsample ${bsample} \
    --gatk_vcf_dir ${gatk_vcf_dir} \
    --dbsnp ${dbsnp} \
    --threads 8 \
    --memory "60G"

# snakemake --dag -s "${base}/scripts/Snakefile_STAR_VCF_per_case" --configfile ${config}| dot -Tpng > workflow_dag.png

snakemake \
    --cores ${SLURM_NTASKS} \
    --use-conda \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --configfile ${config} \
    --snakefile "${base}/scripts/snake_TORSMIC_each_sample"


### test QC
# config=$base/config/test_qc.json

# python ${project_dir}/scripts/qc/make_snakemake_config.py \
#     --project_dir ${project_dir} \
#     --out ${config} \
#     --outdir ${base} \
#     --Bioproject "PRJEB53653" \
#     --Normal_Run "ERR9923018" \
#     --Tumor_Run "ERR9923074" \
#     --CaseName "24" \
#     --metadata ${project_dir}/metadata/data_collection_new.csv \
#     --readlength ${project_dir}/metadata/data_new_readlength.csv \
#     --threads 8 \
#     --memory "60G"

# snakemake \
#     --cores ${SLURM_NTASKS} \
#     --use-conda \
#     --rerun-incomplete \
#     --rerun-triggers mtime \
#     --configfile ${config} \
#     --snakefile ${project_dir}"/scripts/qc/Snakefile"