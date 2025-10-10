#!/home/amunzur/anaconda3/envs/snakemake/bin/python
# ca snakemake

"""
Generates a SLURM batch script to run CNVkit's coverage command on a BAM file.

Inputs:
- BAM file for a single sample
- BED file with probe locations (target regions)
- Output directory for `.cnn` file
- Directories for batch script and SLURM logs
"""
    
import os
import sys
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Run cnvkit coverage on single sample.")
    parser.add_argument("--path_bam", required=True, help="Abs path to bam.")
    parser.add_argument("--dir_output", required=True, help="Outputted .cnn file will be saved here. Abs path please.")
    parser.add_argument("--path_panel", required=True, help="Path to bed file with locations to the probes.")
    parser.add_argument("--dir_batch_scripts", required=True, help="Generated sbatch file will be saved here.")
    parser.add_argument("--dir_logs", required=True, help="Where the Slurm error messages and outputs will be saved.")
    parser.add_argument("--sbatch_time_string", required=False, default="29:00", help="Timelimit for sbatch file. Will default to 30 mins.")
    
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    path_bam=args.path_bam
    dir_output=args.dir_output
    path_panel=args.path_panel
    dir_batch_scripts=args.dir_batch_scripts
    dir_logs=args.dir_logs
        
    sample_name=os.path.basename(path_bam).replace(".bam", "")
    path_sbatch=os.path.join(dir_batch_scripts, f"{sample_name}.batch")
    path_logs=os.path.join(dir_logs, f"{sample_name}.log")
    path_output=os.path.join(dir_output, f"{sample_name}.cnn")

    cnvkit_command=f"cnvkit.py coverage {path_bam} {path_panel} -o {path_output} -p 1 2> {path_logs}"
    
    os.makedirs(dir_batch_scripts, exist_ok=True)
        
    if os.path.exists(path_sbatch): os.remove(path_sbatch)
    jobname = f"cnvkit_cov_{sample_name}"
    
    with open(path_sbatch, 'a') as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --job-name={jobname}\n')
        file.write('#SBATCH --cpus-per-task=2\n')
        file.write('#SBATCH --mem=4G\n')
        file.write('#SBATCH --time=35:00\n')
        file.write(f'#SBATCH --error {path_logs}\n')
        file.write(f'#SBATCH --output {path_logs}\n')
        file.write('\n')
        file.write('source /home/amunzur/anaconda3/etc/profile.d/conda.sh\n')
        file.write('conda activate cnv\n')
        file.write('\n')
        file.write('# Run CNVKIT coverage\n')
        file.write(cnvkit_command)
        
        print(f"sbatch {path_sbatch}")

if __name__ == "__main__":
    main()

# Single sample
# /path/to/amunzur/toolkit/SLURM_run_cnvkit_coverage.py \
# --path_bam /path/to/amunzur/lu_chip/results/bam/SSCS2_final/sample_name.bam \
# --dir_output /path/to/amunzur/lu_chip/results/copynum/coverage/snps \
# --path_panel /path/to/jbacon/reference/bed_files/all/SNP_spikein_probes.bed \
# --dir_batch_scripts /path/to/amunzur/lu_chip/workflow/scripts/batch_scripts/cnvkit_coverage_snps \
# --dir_logs /path/to/amunzur/lu_chip/results/logs_slurm/cnvkit_coverage_snps

# Multiple samples
# while IFS=$'\t' read -r sample_name; do
#     # Skip blank lines or header
#     if [[ -n "$sample_name" && "$sample_name" != "sample_names" ]]; then
#         /path/to/amunzur/toolkit/SLURM_run_cnvkit_coverage.py \
#         --path_bam /path/to/amunzur/lu_chip/results/bam/crpc_panel/${sample_name}.bam \
#         --dir_output /path/to/amunzur/lu_chip/results/copynum/coverage/snps \
#         --path_panel /path/to/jbacon/reference/bed_files/all/SNP_spikein_probes.bed \
#         --dir_batch_scripts /path/to/amunzur/lu_chip/workflow/scripts/batch_scripts/cnvkit_coverage_snps \
#         --dir_logs /path/to/amunzur/lu_chip/results/logs_slurm/cnvkit_coverage_snps
#     fi
# done < path_to_samples_list.tsv