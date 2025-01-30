#!/usr/bin/env python3

"""
For a given list of samples, generate sbatch scripts to run gridss on them.

/groups/wyattgrp/users/amunzur/toolkit/SLURM_run_gridss.py \
--dir_working /groups/wyattgrp/users/amunzur/ecdna_project/sv \
--dir_logs /groups/wyattgrp/users/amunzur/ecdna_project/logs/gridss \
--path_wbc_bam /groups/wyattgrp/users/amunzur/ecdna_project/alignments/TheraP-284_WBC-2019Aug19.bam \
--path_cfdna_bam /groups/wyattgrp/users/amunzur/ecdna_project/alignments/TheraP-284_cfDNA-2020Jul21.bam \
--dir_batch_scripts /groups/wyattgrp/users/amunzur/ecdna_project/scripts/batch_scripts_gridss
"""

import pandas as pd
import os
import argparse
import sys

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Run gridss.")
    parser.add_argument("--dir_working", required=True, help="Some intermediate gridss files will be saved here.")
    parser.add_argument("--dir_logs", required=True, help="Where the Slurm error messages and outputs will be saved.")
    parser.add_argument("--path_wbc_bam", required=True, help="Abs path to the WBC bam")
    parser.add_argument("--path_cfdna_bam", required=True, help="Abs path to the cfDNA bam")
    parser.add_argument("--dir_batch_scripts", required=True, help="Abs path to the batch file to be run with Slurm.")
    
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    sample_name = os.path.basename(args.path_cfdna_bam).replace(".bam", "")
    dir_working_sample=os.path.join(args.dir_working, sample_name)
    
    path_log=os.path.join(args.dir_logs, sample_name+".log")
    path_assembly=os.path.join(dir_working_sample, sample_name+"_assembly.bam")
    path_output=os.path.join(dir_working_sample, sample_name+".vcf")
    path_batch_file=os.path.join(dir_working_sample, sample_name+".batch")
    
    if not os.path.exists(dir_working_sample):
        os.makedirs(dir_working_sample)
    
    path_reference="/groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa"
    
    gridss_command=f"gridss --jar ~/anaconda3/envs/gridss/share/gridss-2.13.2-3/gridss.jar \
    --workingdir {dir_working_sample} \
    --threads 16 \
    --jvmheap 32G \
    --otherjvmheap 32G \
    --reference {path_reference} \
    --output {path_output} \
    --assembly {path_assembly} \
    {args.path_wbc_bam} \
    {args.path_cfdna_bam}"
    
    # WRITE TO BATCH FILE
    if not os.path.exists(args.dir_batch_scripts):
        os.makedirs(args.dir_batch_scripts)
    
    if os.path.exists(path_batch_file):
        os.remove(path_batch_file)
    
    with open(path_batch_file, 'a') as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --job-name=GRIDSS\n')
        file.write('#SBATCH --cpus-per-task=16\n')
        file.write('#SBATCH --mem=256\n')
        file.write('#SBATCH --time=UNLIMITED\n')
        file.write(f'#SBATCH --error {path_log}\n')
        file.write(f'#SBATCH --output {path_log}\n')
        file.write('\n')
        file.write('source /home/amunzur/anaconda3/etc/profile.d/conda.sh\n')
        file.write('conda activate gridss\n')
        file.write('\n')
        file.write('#RUN GRIDSS\n')
        file.write('\n')
        file.write(gridss_command)
    
    print(" ")
    print(f"BATCH: {path_batch_file}")

if __name__ == "__main__":
    main()