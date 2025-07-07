#!/usr/bin/env python3

"""
For a given list of samples, generate sbatch scripts to run gridss on them.

/groups/wyattgrp/users/amunzur/toolkit/SLURM_run_gridss.py \
--dir_working /groups/wyattgrp/users/amunzur/ecdna_project/sv \
--dir_logs /groups/wyattgrp/users/amunzur/ecdna_project/logs/gridss \
--path_wbc_bam /groups/wyattgrp/users/amunzur/ecdna_project/alignments/TheraP-284_WBC-2019Aug19.bam \
--path_cfdna_bam /groups/wyattgrp/users/amunzur/ecdna_project/alignments/TheraP-284_cfDNA-2020Jul21.bam \
--dir_batch_scripts /groups/wyattgrp/users/amunzur/ecdna_project/scripts/batch_scripts/gridss
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
    dir_working_gripss=os.path.join(args.dir_working, sample_name, "gripss")
    
    path_log=os.path.join(args.dir_logs, sample_name+".log")
    path_assembly=os.path.join(dir_working_sample, sample_name+"_assembly.bam")
    path_output=os.path.join(dir_working_sample, sample_name+".vcf")
    path_batch_file=os.path.join(args.dir_batch_scripts, sample_name+".batch")
    
    if not os.path.exists(dir_working_sample):
        os.makedirs(dir_working_sample)
    
    if not os.path.exists(dir_working_gripss):
        os.makedirs(dir_working_gripss)
    
    path_reference="/groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa"
    pon_sgl_file="/groups/wyattgrp/users/amunzur/useful_data/gripss_pon/v5_34/ref/38/sv/sgl_pon.38.bed.gz"
    pon_sv_file="/groups/wyattgrp/users/amunzur/useful_data/gripss_pon/v5_34/ref/38/sv/sv_pon.38.bedpe.gz"
    known_hotspot_file="/groups/wyattgrp/users/amunzur/useful_data/gripss_pon/v5_34/ref/38/sv/known_fusions.38.bedpe"
    repeat_mask_file="/groups/wyattgrp/users/amunzur/useful_data/gripss_pon/v5_34/ref/38/sv/repeat_mask_data.38.fa.gz"
    
    gridss_command=f"gridss --jar ~/anaconda3/envs/gridss/share/gridss-2.11.1-1//gridss.jar\\\n\
    --workingdir {dir_working_sample}\\\n\
    --threads 16\\\n\
    --jvmheap 32G\\\n\
    --otherjvmheap 32G\\\n\
    --reference {path_reference}\\\n\
    --output {path_output}\\\n\
    --assembly {path_assembly}\\\n\
    {args.path_wbc_bam}\\\n\
    {args.path_cfdna_bam}\n"

    wbc_name=os.path.basename(args.path_wbc_bam).replace(".bam", "")
    cfDNA_name=os.path.basename(args.path_cfdna_bam).replace(".bam", "")

    gripss_command=f'gripss -Xmx64g\\\n\
    -sample {cfDNA_name}\\\n\
    -reference {wbc_name}\\\n\
    -ref_genome_version 38\\\n\
    -ref_genome {path_reference}\\\n\
    -filter_sgls\\\n\
    -pon_sgl_file {pon_sgl_file}\\\n\
    -pon_sv_file {pon_sv_file}\\\n\
    -known_hotspot_file {known_hotspot_file}\\\n\
    -repeat_mask_file {repeat_mask_file}\\\n\
    -vcf {path_output}\\\n\
    -output_dir {dir_working_gripss}\n'

    # WRITE TO BATCH FILE
    if not os.path.exists(args.dir_batch_scripts):
        os.makedirs(args.dir_batch_scripts)
    
    if os.path.exists(path_batch_file):
        os.remove(path_batch_file)
    
    with open(path_batch_file, 'a') as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --job-name=GRIDSS\n')
        file.write('#SBATCH --cpus-per-task=8\n')
        file.write('#SBATCH --mem=64G\n')
        file.write('#SBATCH --time=UNLIMITED\n')
        file.write(f'#SBATCH --error {path_log}\n')
        file.write(f'#SBATCH --output {path_log}\n')
        file.write('\n')
        # file.write("start_time=$(date +%s)\n")
        file.write('source /home/amunzur/anaconda3/etc/profile.d/conda.sh\n')
        file.write('\n')
        file.write('#RUN GRIDSS\n')
        file.write('conda activate gridss\n')
        file.write('\n')
        # file.write(gridss_command)
        file.write('\n')
        file.write('#RUN GRIPSS\n')
        file.write('conda activate gripss\n')
        file.write('\n')
        file.write(gripss_command)
        file.write("\n")
        # file.write("end_time=$(date +%s)")
        # file.write("elapsed_time=$((end_time - start_time))\n")
        # file.write(f"echo 'Job finished in $elapsed_time seconds' >> {path_log}\n")

    
    print(f"BATCH: {path_batch_file}")

if __name__ == "__main__":
    main()

