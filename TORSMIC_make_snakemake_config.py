import argparse
import json
import math
import os
import re
import sys

# this script makes a config file for snakemake in JSON
# snakefile for germline/somantic variant calling pipeline

def main():

    args = parse_args()

    config = {}
    config["base"] = args.base
    ## envs
    config["envs"] = "TORSMIC"
    ## parameters for inputs
    config["in"] = {
        "package_location":args.package_location,
        "bio_tumor": args.bio_tumor,
        "bsample": args.bsample,
        "annovar_index": args.annovar_dir,
        "dbsnp" : args.dbsnp,
        "scripts_location" :args.package_location+"/"+"scripts",
        "data_source_location":args.package_location+"/"+"data_source",
        "gatk_vcf_dir" : args.gatk_vcf_dir
    }

    ## paths for outputs
    config["out"] = {
        "somatic_output_dir": args.base+"/data/"+args.bio_tumor+"/"+"somatic_out"
    }

    ## parameters for cluster job resources/envs
    mem_int = int(re.sub("G", "", args.memory))
    mem_2parallel = math.floor(mem_int * 0.5)
    mem_4parallel = math.floor(mem_int * 0.25)
    config["resources"] = {
        "threads": args.threads,
        "mem": args.memory,
        "threads_2parallel": math.floor(args.threads * 0.5),
        "mem_2parallel": str(mem_2parallel)+"G",
        "threads_4parallel": math.floor(args.threads * 0.25),
        "mem_4parallel": str(mem_4parallel)+"G"
    }


    
    with open(args.out, "w") as conf:
        json.dump(config, conf, indent=4)



def parse_args():
    parser = argparse.ArgumentParser(prog='TORSMIC_make_snakemake_config.py', description="Create config file in json format for snakemake pipeline.")
    ## required ##
    essential_args = parser.add_argument_group('Essential')
    essential_args.add_argument("--base", type=str, help="Path to the outmost directory", required=True)
    essential_args.add_argument("--package_location", type=str, help="Path to TORSMIC directory", required=True)
    essential_args.add_argument("--gatk_vcf_dir", type=str, help="VCF file after GATK filtering.", required=True)
    essential_args.add_argument("--bio_tumor", type=str, help="bio_tumor source of the data.", required=True)
    essential_args.add_argument("--bsample", type=str, help="Biosample name", required=True)
    essential_args.add_argument("--dbsnp", type=str, help="A dbsnp file", required=True)
    essential_args.add_argument("--out", type=str, help="File name of the output config file in json. Required.", required=True)
    ## optional ##
    optional_args = parser.add_argument_group('Optional')
    optional_args.add_argument("--annovar_dir", type=str, help="Annovar script folder.", required=False)
    ## resources options ##
    resources_args = parser.add_argument_group('Optional resources')
    resources_args.add_argument("--threads", type=int, help="The number of processors to use for individual cluster jobs. [default = 4]", default=4, required=False)
    resources_args.add_argument("--memory", type=str, help="The number of memory in 'G' to use for individual cluster jobs. [default = '64G']", default="64G", required=False)

    args = parser.parse_args()

    ## required ##
    # parse project repo path
    args.base = os.path.abspath(args.base)
    args.package_location = os.path.abspath(args.package_location)
    args.gatk_vcf_dir = os.path.abspath(args.gatk_vcf_dir)
    args.dbsnp = os.path.abspath(args.dbsnp)
    # parse output path for the json config
    args.out = os.path.abspath(args.out)


    ## optional ##
    # parse reference
    # parse annovar_dir
    if args.annovar_dir is not None:
        args.annovar_dir = os.path.abspath(args.annovar_dir)
    else:
        args.annovar_dir = os.path.abspath("/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf")

    # parse threads
    args.threads = int(args.threads)
    
    # parse memory
    if re.search("G$", args.memory) is None:
        sys.stderr.write("Only memory in 'G' is allowed for --memory option, i.e, '20G'")
        sys.exit(1)


    return args


if __name__ == "__main__":
    main()