#!/usr/bin/env python3

"""
Generate batch script to LOHHLA.

How to use:
/groups/wyattgrp/users/amunzur/toolkit/SLURM_run_lohhla.py \
    --tumor_name 
    --outputDir 
    --BAMDir 
    --normalBAMfile 
    --hlaPath 
    --hla_fasta 
    --CopyNumLoc 
    --path_log 
    --path_batch_file 
"""

import pandas as pd
import os
import argparse
import sys

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Run gridss.")
    parser.add_argument("--tumor_name", required=True, help="Name of the tumor sample, withouth the .bam extension.")
    parser.add_argument("--outputDir", required=True, help="LOHHLA outputs will be saved here.")
    parser.add_argument("--BAMDir", required=True, help="Where LOHHLA checks for tumor and normal bams, with bai files.")
    parser.add_argument("--normalBAMfile", required=True, help="Abs path to the normal bam, must be in the BAMDir.")
    parser.add_argument("--hlaPath", required=True, help="Absolute path to the winners file from Polysolver.")
    parser.add_argument("--hla_fasta", required=True, help="Absolute path to patient specific HLA fasta file.")
    parser.add_argument("--CopyNumLoc", required=True, help="Absolute path to ploidy and purity estimates, in the format LOHHLA requires.")
    parser.add_argument("--path_log", required=True, help="Absolute path to the SLURM log file.")
    parser.add_argument("--path_batch_file", required=True, help="Absolute path to the SLURM batch file.")
    
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    gatkDir="/home/amunzur/anaconda3/envs/polysolver/jar"
    novoDir="/home/amunzur/anaconda3/envs/lohhla/bin"
    
    lohhla_command = f"Rscript /groups/wyattgrp/users/amunzur/gillian_proj/lohhla/LOHHLAscript.R \
      --patientId {args.tumor_name} \
      --outputDir {args.outputDir} \
      --BAMDir {args.BAMDir} \
      --normalBAMfile {args.normalBAMfile} \
      --hlaPath {args.hlaPath} \
      --CopyNumLoc {args.CopyNumLoc} \
      --HLAfastaLoc {args.hla_fasta} \
      --mappingStep TRUE \
      --minCoverageFilter 10 \
      --fishingStep TRUE \
      --cleanUp FALSE \
      --gatkDir {gatkDir} \
      --novoDir {novoDir}"
    
    # WRITE TO BATCH FILE
    if not os.path.exists(os.path.dirname(args.path_batch_file)):
        os.makedirs(os.path.dirname(args.path_batch_file))
    
    if os.path.exists(args.path_batch_file):
        os.remove(args.path_batch_file)
    
    with open(args.path_batch_file, 'a') as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --job-name=LOHHLA\n')
        file.write('#SBATCH --cpus-per-task=8\n')
        file.write('#SBATCH --mem=8G\n')
        file.write('#SBATCH --time=UNLIMITED\n')
        file.write(f'#SBATCH --error {args.path_log}\n')
        file.write(f'#SBATCH --output {args.path_log}\n')
        file.write('\n')
        file.write('source /home/amunzur/anaconda3/etc/profile.d/conda.sh\n')
        file.write('conda activate lohhla\n')
        file.write('\n')
        file.write('#RUN LOHHLA\n')
        file.write('\n')
        file.write(lohhla_command)
    
    print(" ")
    print(f"BATCH: {args.path_batch_file}")


if __name__ == "__main__":
    main()
