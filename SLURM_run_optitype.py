#!/usr/bin/env python3

"""
For a given list of samples, generates individual sbatch scripts to run optitype.

/groups/wyattgrp/users/amunzur/toolkit/SLURM_run_optitype.py \
--dir_fastqs /groups/wyattgrp/users/amunzur/hla_pipeline/results/data/fq/raw \
--sample_name GU-18-296_WBC_09Jul2018 \
--dir_output_main /groups/wyattgrp/users/amunzur/hla_pipeline/results/optitype \
--dir_batch_scripts /groups/wyattgrp/users/amunzur/hla_pipeline/workflow/batch_scripts/optitype \
--dir_logs /groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/optitype \
--threads 16 \
--memory 64G \
--timelimit 23:00:00
"""

import os
import argparse
import sys

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Run Optitype.")
    parser.add_argument("--dir_fastqs", required=True, help="Absolute path where all FastQ files are stored.")
    parser.add_argument("--sample_name", required=True, help="Sample name.")
    parser.add_argument("--dir_output_main", required=True, help="Main output directory, will make a new sample directory within here.")
    parser.add_argument("--dir_batch_scripts", required=True, help="Main directory where all batch scripts will be saved.")    
    parser.add_argument("--dir_logs", required=True, help="Main directory where all batch scripts will be saved.")    
    parser.add_argument("--threads", required=True, help="Number of threads.")    
    parser.add_argument("--memory", required=True, help="Memory. Provide with suffix, like 8G")    
    parser.add_argument("--timelimit", required=True, help="Timelimit.")    
    
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    sample_name = os.path.basename(args.sample_name)
    dir_output_sample=os.path.join(args.dir_output_main, sample_name)
    dir_logs=args.dir_logs
    path_log=os.path.join(args.dir_logs, sample_name+".log")
    path_batch_file=os.path.join(args.dir_batch_scripts, sample_name+".batch")
    
    dir_fastqs=args.dir_fastqs
    path_fastq1=os.path.join(dir_fastqs, f"{sample_name}_1.fq.gz")
    path_fastq2=os.path.join(dir_fastqs, f"{sample_name}_1.fq.gz")
        
    if not os.path.exists(dir_output_sample):
        os.makedirs(dir_output_sample)
    
    if not os.path.exists(dir_logs):
        os.makedirs(dir_logs)
    
    if not os.path.exists(args.dir_batch_scripts):
        os.makedirs(args.dir_batch_scripts)
    
    command=f"python ~/anaconda3/envs/mhc2_genotyping_optitype/bin/OptiTypePipeline.py \
        -i {path_fastq1} {path_fastq2} \
        --dna \
        --enumerate 2 \
        --verbose \
        --outdir {dir_output_sample}"
    
    if os.path.exists(path_batch_file):
        os.remove(path_batch_file)
    
    with open(path_batch_file, 'a') as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --job-name=WGS_OPI_{sample_name}\n')
        file.write(f'#SBATCH --cpus-per-task={args.threads}\n')
        file.write(f'#SBATCH --mem={args.memory}\n')
        file.write(f'#SBATCH --time={args.timelimit}\n')
        file.write(f'#SBATCH --error {path_log}\n')
        file.write(f'#SBATCH --output {path_log}\n')
        file.write('\n')
        file.write('source /home/amunzur/anaconda3/etc/profile.d/conda.sh\n')
        file.write('conda activate mhc2_genotyping_optitype\n')
        file.write('\n')
        file.write('#RUN OPTITYPE\n')
        file.write('\n')
        file.write(command)
    
    print(f"BATCH: {path_batch_file}")

if __name__ == "__main__":
    main()

# Generate batch scripts for multiple samples at once
# while IFS= read -r line || [[ -n "$line" ]]; do
#     echo "Processing: $line"
#     /groups/wyattgrp/users/amunzur/toolkit/SLURM_run_optitype.py \
#     --dir_fastqs /groups/wyattgrp/users/amunzur/hla_pipeline/results/data/fq/raw \
#     --sample_name ${line} \
#     --dir_output_main /groups/wyattgrp/users/amunzur/hla_pipeline/results/optitype \
#     --dir_batch_scripts /groups/wyattgrp/users/amunzur/hla_pipeline/workflow/batch_scripts/optitype \
#     --dir_logs /groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/optitype \
#     --threads 16 \
#     --memory 64G \
#     --timelimit 23:00:00
# done < /groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test_WBC_samples.tsv


# for file in /groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/bams_without_alt_contigs/alignments/*WBC*WES*.bam; do
# 	/groups/wyattgrp/users/amunzur/toolkit/SLURM_run_optitype.py \
# 	--dir_fastqs /groups/wyattgrp/users/amunzur/hla_pipeline/results/data/fq/raw \
# 	--sample_name $(basename "${file}" .bam) \
# 	--dir_output_main /groups/wyattgrp/users/amunzur/hla_pipeline/results/optitype \
# 	--dir_batch_scripts /groups/wyattgrp/users/amunzur/hla_pipeline/workflow/batch_scripts/optitype \
# 	--dir_logs /groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/optitype \
# 	--threads 16 \
# 	--memory 64G \
# 	--timelimit 23:00:00
# done
