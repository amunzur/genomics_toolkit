#!/home/amunzur/anaconda3/envs/snakemake/bin/python
"""
For a given list of samples, generate per sample batch scripts to run Vardict. This is single sample mode.
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
    parser = argparse.ArgumentParser(description="Run AA pipeline.")
    parser.add_argument("--path_hg38", required=True, help="Path to the reference fasta.")
    parser.add_argument("--threshold_min_vaf", required=True, help="Minimum VAF required to detect mutations.")
    parser.add_argument("--path_bam", required=True, help="Path to the bam file of interest for Vardict.")
    parser.add_argument("--path_bed", required=True, help="Bed file to regions of interest.")
    parser.add_argument("--dir_output", required=True, help="Outputted VCF will be saved here.")
    parser.add_argument("--dir_batch_scripts", required=True, help="Batch script will be written here.")
    parser.add_argument("--dir_logs", required=True, help="Vardict logs will be here.")

    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    path_hg38=args.path_hg38
    threshold_min_vaf=args.threshold_min_vaf
    path_bam=args.path_bam
    path_bed=args.path_bed
    sample_name=os.path.basename(path_bam).replace(".bam", "")
    path_output=os.path.join(args.dir_output, sample_name+".vcf")
    dir_batch_scripts=args.dir_batch_scripts
    dir_logs=args.dir_logs
    #
    path_sbatch=os.path.join(dir_batch_scripts, sample_name+".batch")
    path_logs=os.path.join(dir_logs, sample_name+".log")
    if os.path.exists(path_sbatch): os.remove(path_sbatch)
    jobname = "freebayes"
    
    os.makedirs(args.dir_batch_scripts, exist_ok=True)
    os.makedirs(args.dir_output, exist_ok=True)
    os.makedirs(dir_logs, exist_ok=True)

    with open(path_sbatch, 'a') as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --job-name={jobname}\n')
        file.write('#SBATCH --cpus-per-task=4\n')
        file.write('#SBATCH --mem=4G\n')
        file.write('#SBATCH --time=30:00\n')
        file.write(f'#SBATCH --error {path_logs}\n')
        file.write(f'#SBATCH --output {path_logs}\n')
        file.write('\n')
        file.write('# Run Vardict\n')
        file.write('source /home/amunzur/anaconda3/etc/profile.d/conda.sh\n')
        file.write("conda activate chip_variantcalling\n")
        file.write(f'freebayes {path_bam} \
                    -f {path_hg38} \
                    -t {path_bed} \
                    --pooled-continuous \
                    --min-alternate-fraction {threshold_min_vaf} \
                    --min-alternate-count 1 | bcftools norm -m-both -f {path_hg38} -o {path_output}')
        print(f"sbatch {path_sbatch}")

if __name__ == "__main__":
    main()

# for file in /groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/sorted/*WBC*WGS*.bam; do
#     if [[ "$file" != *WES* && "$file" != *WGS* ]]; then
#         /groups/wyattgrp/users/amunzur/toolkit/SLURM_run_mutect2.py \
#             --dir_batch_scripts /groups/wyattgrp/users/amunzur/hla_pipeline/workflow/batch_scripts/mutect2 \
#             --path_bam "$file" \
#             --dir_logs /groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/run_mutect2 \
#             --path_hg38 /groups/wyattgrp/users/jbacon/reference/hg38_emseq/GRCh38_no_alt_plus_hs38d1.fa \
#             --threshold_min_vaf 0.35 \
#             --path_bed /groups/wyattgrp/users/amunzur/hla_pipeline/resources/hla_panel/probe_locations.bed \
#             --dir_output /groups/wyattgrp/users/amunzur/hla_pipeline/results/mutect2/snps
#     fi
# done