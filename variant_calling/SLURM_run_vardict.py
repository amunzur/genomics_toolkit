#!/usr/bin/python
"""
For a given list of samples, generate per sample batch scripts to run Vardict. This is single sample mode.
Need to have VarDictJava in home dir or provide the path.
"""
# /groups/wyattgrp/users/amunzur/toolkit/SLURM_run_vardict.py \
    # --dir_batch_scripts /groups/wyattgrp/users/amunzur/hla_pipeline/workflow/batch_scripts/vardict \
    # --path_bam /groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/sorted/GU-19-534_WBC_03Sep2019.bam \
    # --dir_logs /groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/run_vardict \
    # --path_hg38 /groups/wyattgrp/users/jbacon/reference/hg38_emseq/GRCh38_no_alt_plus_hs38d1.fa \
    # --threshold_min_vaf 0.35 \
    # --min_alt_reads 5 \
    # --path_bed /groups/wyattgrp/users/amunzur/hla_pipeline/resources/hla_panel/probe_locations.bed \
    # --dir_output /groups/wyattgrp/users/amunzur/hla_pipeline/results/vardict/wbc_snps \

import os
import sys
import argparse

def main():
    # Argument parsing
    home_dir=os.path.expanduser("~")
    
    parser = argparse.ArgumentParser(description="Run VarDictJava on a single sample.")
    parser.add_argument("--path_hg38", required=True, help="Path to the reference fasta.")
    parser.add_argument("--threshold_min_vaf", required=True, help="Minimum VAF required to detect mutations.")
    parser.add_argument("--min_alt_reads", required=True, help="Minimum number of alt reads needed to detect a mutation.")
    parser.add_argument("--path_bam", required=True, help="Path to the bam file of interest for Vardict.")
    parser.add_argument("--path_bed", required=True, help="Bed file to regions of interest.")
    parser.add_argument("--dir_output", required=True, help="Outputted VCF will be saved here.")
    parser.add_argument("--dir_batch_scripts", required=True, help="Batch script will be written here.")
    parser.add_argument("--dir_logs", required=True, help="Vardict logs will be here.")
    parser.add_argument("--dir_vardictjava", required=False, default=home_dir, help="Path to vardict. Defaults to home directory.")
    parser.add_argument("--sbatch_time_string", required=False, default="29:00", help="Timelimit for sbatch file. Will default to 30 mins.")
    parser.add_argument("--sbatch_partition", required=False, default="", help="Partition. Choose from: long, big-mem, normal, express, debug. Defaults to long.")
    
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    path_hg38=args.path_hg38
    threshold_min_vaf=args.threshold_min_vaf
    min_alt_reads=args.min_alt_reads
    path_bam=args.path_bam
    path_bed=args.path_bed
    sample_name=os.path.basename(path_bam).replace(".bam", "")
    path_output=os.path.join(args.dir_output, sample_name+".vcf")
    dir_batch_scripts=args.dir_batch_scripts
    dir_logs=args.dir_logs
    dir_vardictjava=args.dir_vardictjava
    sbatch_time_string=args.sbatch_time_string
    sbatch_partition=args.sbatch_partition
    
    #
    path_sbatch=os.path.join(dir_batch_scripts, sample_name+".batch")
    path_logs=os.path.join(dir_logs, sample_name+".log")
    if os.path.exists(path_sbatch): os.remove(path_sbatch)
    jobname = "Vardict"

    with open(path_sbatch, 'a') as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --job-name={jobname}\n')
        file.write('#SBATCH --cpus-per-task=4\n')
        file.write('#SBATCH --mem=4G\n')
        file.write(f'#SBATCH --time={sbatch_time_string}\n')
        file.write(f'#SBATCH --error {path_logs}\n')
        file.write(f'#SBATCH --output {path_logs}\n')
        if sbatch_partition!="":
            file.write(f'#SBATCH --partition={sbatch_partition}\n')
        
        file.write('\n')
        file.write('export PATH="/usr/local/bin:$PATH"') # This ensure Rscript is correctly located
        file.write('\n')
        file.write('# Run Vardict\n')
        file.write(f'/home/amunzur/VarDictJava/VarDict/vardict \
                   -G {path_hg38} \
                    -f {threshold_min_vaf} \
                    -N {sample_name} \
                    -r {min_alt_reads} \
                    -b {path_bam} \
                    -k 0 -c 1 -S 2 -E 3 -g 4 {path_bed} | \
                    {dir_vardictjava}/VarDictJava/build/install/VarDict/bin/teststrandbias.R | \
                    {dir_vardictjava}/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl > {path_output}')
        print(f"sbatch {path_sbatch}")

if __name__ == "__main__":
    main()

# for file in /groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/sorted/*WBC*WGS*.bam; do
#     if [[ "$file" != *WES* && "$file" != *WGS* ]]; then
#         /groups/wyattgrp/users/amunzur/toolkit/SLURM_run_vardict.py \
#             --dir_batch_scripts /groups/wyattgrp/users/amunzur/hla_pipeline/workflow/batch_scripts/vardict \
#             --path_bam "$file" \
#             --dir_logs /groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/run_vardict \
#             --path_hg38 /groups/wyattgrp/users/jbacon/reference/hg38_emseq/GRCh38_no_alt_plus_hs38d1.fa \
#             --threshold_min_vaf 0.35 \
#             --min_alt_reads 5 \
#             --path_bed /groups/wyattgrp/users/amunzur/hla_pipeline/resources/hla_panel/probe_locations.bed \
#             --dir_output /groups/wyattgrp/users/amunzur/hla_pipeline/results/vardict/wbc_snps
#     fi
# done
