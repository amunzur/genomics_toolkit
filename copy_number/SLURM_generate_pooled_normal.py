#!/home/amunzur/anaconda3/envs/snakemake/bin/python
# ca snakemake
"""
Generates the cnvkit command to generate a pooled normal. You need to have a conda env named cnv with cnvkit installed.
"""
    
import os
import sys
import argparse
import re
import pandas as pd

def main():
    default_conda_path = os.path.expanduser("~/anaconda3/etc/profile.d/conda.sh")
    
    parser = argparse.ArgumentParser(description="Generate cnvkit pooled normal.")
    parser.add_argument("--path_conda", required=False, default=default_conda_path, help=f"Absolute path to your conda.sh file. Defaults to {default_conda_path}")
    parser.add_argument("--path_cnns", required=True, help="Path to a text file with absolute paths to cnn files to use when making the pooled normal, or a string of abs paths separated by a string.")
    parser.add_argument("--path_hg38", required=True, help="The ref genome the cnn files came from.")
    parser.add_argument("--dir_logs", required=True, help="Where the Slurm error messages and outputs will be saved.")
    parser.add_argument("--path_output", required=True, help="Output files will be saved here for each allele type.")
    parser.add_argument("--path_sbatch", required=True, help="Generated sbatch file will be saved here.")
    parser.add_argument("--sbatch_time_string", required=False, default="29:00", help="Timelimit for sbatch file. Will default to 30 mins.")
    parser.add_argument("--sbatch_partition", required=False, default="long", help="Partition. Choose from: long, big-mem, normal, express, debug. Defaults to long.")
        
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    path_conda=args.path_conda
    path_cnns=args.path_cnns
    path_hg38=args.path_hg38
    dir_logs=args.dir_logs
    path_output=args.path_output
    path_sbatch=args.path_sbatch
    sbatch_time_string=args.sbatch_time_string
    sbatch_partition=args.sbatch_partition
    
    # Check if a string of paths or a text file provided.
    if os.path.isfile(path_cnns):
        cnn_file=pd.read_csv(path_cnns, sep="\t")
        path_cnns_list=[line.strip() for line in open(path_cnns)]
        cnn_paths_string=" ".join(path_cnns_list)
    else:
        cnn_paths_string=path_cnns
    
    path_logs=os.path.join(dir_logs, "make_pooled_norm.log")
    
    cnvkit_command=f"cnvkit.py reference {cnn_paths_string} -f {path_hg38} -o {path_output} --male-reference --sample-sex male"
    os.makedirs(os.path.dirname(path_sbatch), exist_ok=True)
        
    if os.path.exists(path_sbatch): os.remove(path_sbatch)
    jobname = f"Make_pool_norm"
    
    with open(path_sbatch, 'a') as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --job-name={jobname}\n')
        file.write('#SBATCH --cpus-per-task=4\n')
        file.write('#SBATCH --mem=8G\n')
        file.write(f'#SBATCH --time={sbatch_time_string}\n')
        file.write(f'#SBATCH --error {path_logs}\n')
        file.write(f'#SBATCH --output {path_logs}\n')
        file.write(f'#SBATCH --partition={sbatch_partition}\n')
        file.write('\n')
        file.write(f'source {path_conda}\n')
        file.write('conda activate cnv\n')
        file.write('\n')
        file.write('# Run CNVKIT reference\n')
        file.write(cnvkit_command)
        
        print(f"BATCH: {path_sbatch}")

if __name__ == "__main__":
    main()

# /groups/wyattgrp/users/amunzur/toolkit/SLURM_generate_pooled_normal.py \
# --path_cnns /groups/wyattgrp/users/amunzur/lu_chip/results/copynum/path_for_pooled_norm.txt \
# --path_hg38 /groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa \
# --dir_logs /groups/wyattgrp/users/amunzur/lu_chip/results/logs_slurm/cnvkit_pooled_normal \
# --path_output /groups/wyattgrp/users/amunzur/lu_chip/results/copynum/pooled_normal_all_wbc_except_154.cnn \
# --path_sbatch /groups/wyattgrp/users/amunzur/lu_chip/workflow/scripts/batch_scripts/cnvkit_generate_pooled_normal/script.sbatch