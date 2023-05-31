#!/bin/bash
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=1           # Number of CPU cores per task (4)
#SBATCH --mem=50G                   # Job memory limit (10 GB)
#SBATCH --time=10:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=Step1_extract_somatic.%j.out    # Standard output log
#SBATCH --error=Step1_extract_somatic.%j.err     # Standard error log

#### The pathes that need to be changed
package_location='/home/kh31516/kh31516_Lab_Share_script/IdentifiyNeoAntigene/TORSMIC' ## the location where you put the package
new_add_data='/scratch/kh31516/UGA/UGA_OSA_somatic_germline.txt.gz'                  ## the location of the table which you concat all of the results from all the samples
pipeline_out_file_name='/scratch/kh31516/UGA/UGA_OSA_pipeline_filtering.txt' ## the output file from the pipeline without maching learning prediction
ml_out_file_name='/scratch/kh31516/UGA/UGA_OSA_ml_filtering.txt' ## the output file from the pipeline with maching learning prediction
####

scripts_location=${package_location}/scripts
data_source_location=${package_location}/data_source

ml R/4.0.0-foss-2019b
# package_location <- as.character(args[1])
# new_data_file <- as.character(args[2])
# out_file_name <- as.character(args[3])

Rscript --vanilla ${scripts_location}/Sapelo2_add_new_data_whole_pipeline.R \
${package_location} \
${new_add_data} \
${pipeline_out_file_name} 

ml Miniconda3/4.9.2
source activate /home/kh31516/myenv/py38

cd $scripts_location

# model_folder = sys.argv[1]
# data_test_file = sys.argv[2]
# final_output = sys.argv[3]
python ${scripts_location}/Sapelo2_ML_test_pipeline_data.py \
${data_source_location} \
${pipeline_out_file_name}  \
${ml_out_file_name}