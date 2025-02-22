#!/usr/bin/env python3

"""
For a given list of samples, generate sbatch scripts to run lilac on them.

/groups/wyattgrp/users/amunzur/toolkit/SLURM_run_lilac.py \
--path_wbc_bam  \
--dir_logs /groups/wyattgrp/users/amunzur/ecdna_project/logs/lilac \
--dir_output_main /groups/wyattgrp/users/amunzur/hla_pipeline/results/lilac \
--dir_batch_scripts /groups/wyattgrp/users/amunzur/hla_pipeline/workflow/batch_scripts/lilac \
--threads 4 \
--memory 16G \
--timelimit 29:00
"""

import os
import argparse
import sys

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Run gridss.")
    parser.add_argument("--path_wbc_bam", required=True, help="Abs path to the WBC bam")
    parser.add_argument("--dir_logs", required=True, help="Directory where the sample log files will be saved.")
    parser.add_argument("--dir_output_main", required=True, help="Output directory, will make a new sample directory within here.")
    parser.add_argument("--dir_batch_scripts", required=True, help="Abs path to the batch file to be run with Slurm.")    
    parser.add_argument("--threads", required=True, help="Number of threads.")    
    parser.add_argument("--memory", required=True, help="Memory. Provide with suffix, like 8G")    
    parser.add_argument("--timelimit", required=True, help="Timelimit.")    
    
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    sample_name = os.path.basename(args.path_wbc_bam).replace(".bam", "")
    dir_output_sample=os.path.join(args.dir_output_main, sample_name)
    path_log=os.path.join(args.dir_logs, sample_name+".log")
    path_batch_file=os.path.join(args.dir_batch_scripts, sample_name+".batch")
    # threads=args.threads
    wbc_bam=args.path_wbc_bam
    
    if not os.path.exists(dir_output_sample):
        os.makedirs(dir_output_sample)
    
    path_reference_no_alt_contigs="/groups/wyattgrp/users/jbacon/reference/hg38_emseq/GRCh38_no_alt_plus_hs38d1.fa"
    lilac_resource_dir="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/lilac"
    
    command=f"java -jar /home/amunzur/lilac_v1.5.2.jar \
        -sample {sample_name} \
        -ref_genome {path_reference_no_alt_contigs} \
        -ref_genome_version V38 \
        -resource_dir {lilac_resource_dir} \
        -reference_bam {wbc_bam} \
        -output_dir {dir_output_sample} \
        -threads {args.threads}"
    
    # WRITE TO BATCH FILE
    if not os.path.exists(args.dir_batch_scripts):
        os.makedirs(args.dir_batch_scripts)
    
    if os.path.exists(path_batch_file):
        os.remove(path_batch_file)
    
    with open(path_batch_file, 'a') as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --job-name={sample_name}\n')
        file.write(f'#SBATCH --cpus-per-task={args.threads}\n')
        file.write(f'#SBATCH --mem={args.memory}\n')
        file.write(f'#SBATCH --time={args.timelimit}\n')
        file.write(f'#SBATCH --error {path_log}\n')
        file.write(f'#SBATCH --output {path_log}\n')
        file.write('\n')
        file.write('source /home/amunzur/anaconda3/etc/profile.d/conda.sh\n')
        file.write('conda activate lilac\n')
        file.write('\n')
        file.write('#RUN LILAC\n')
        file.write('\n')
        file.write(command)
    
    print(f"BATCH: {path_batch_file}")

if __name__ == "__main__":
    main()





