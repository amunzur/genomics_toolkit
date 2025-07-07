#!/home/amunzur/anaconda3/envs/snakemake/bin/python

"""
Generates a SLURM batch script to run AmpliconArchitect (AA) on a cfDNA sample using matched WBC as normal.

Inputs:
- cfDNA BAM file
- WBC BAM file
- GRIPSS VCF (structural variants)
- Segmented copy number file (with integer copy states)
- Downsampling depth (e.g. 20 for 20X; use -1 to skip)
- Output directories for AA results, logs, and batch scripts

Outputs:
- A SLURM-compatible bash script to run AA and AmpliconClassifier
- A log file path for SLURM
- A cfDNA-only subset of the GRIPSS VCF

Usage:
python SLURM_run_AA.py \
  --dir_batch_scripts /path/to/scripts \
  --path_wbc_bam /path/to/WBC.bam \
  --path_cfdna_bam /path/to/cfDNA.bam \
  --dir_logs /path/to/logs \
  --dir_AA_output /path/to/output \
  --downsampling_depth 20 \
  --path_gripss_vcf /path/to/gripss.filtered.vcf.gz \
  --path_segmented_cn /path/to/segmented_copy_number.tsv
"""
    
import os
import sys
import argparse
import re

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Run AA pipeline.")
    parser.add_argument("--dir_batch_scripts", required=True, help="X")
    parser.add_argument("--path_wbc_bam", required=True, help="Abs path to WBC bam.")
    parser.add_argument("--path_cfdna_bam", required=True, help="Abs path to cfDNA bam.")
    parser.add_argument("--dir_logs", required=True, help="Where the Slurm error messages and outputs will be saved.")
    parser.add_argument("--dir_AA_output", required=True, help="Working dir for AA. Outputs will be saved here.")
    parser.add_argument("--downsampling_depth", required=True, help="-1 for no downsampling. Give integer. '20' for 20x.")
    parser.add_argument("--path_gripss_vcf", required=True, help="Path to the GRIPSS filtered vcf for SVs.")
    parser.add_argument("--path_segmented_cn", required=True, help="Path to the segmented CN file. Last col need to be integer copy states.")   

    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    dir_batch_scripts=args.dir_batch_scripts
    path_wbc_bam=args.path_wbc_bam
    path_cfdna_bam=args.path_cfdna_bam
    dir_AA_output=args.dir_AA_output
    downsampling_depth=args.downsampling_depth
    path_gripss_vcf=args.path_gripss_vcf
    path_segmented_cn=args.path_segmented_cn
    dir_logs=args.dir_logs

    PATH_gripss_vcf_cfDNA_only=path_gripss_vcf.replace("gripss.filtered.vcf.gz", "cfDNA_only.gripss.filtered.vcf")

    cfdna_sample=os.path.basename(path_cfdna_bam).replace(".bam", "")
    path_sbatch=os.path.join(dir_batch_scripts, cfdna_sample+".batch")
    path_logs=os.path.join(dir_logs, cfdna_sample+".log")
    os.makedirs(dir_AA_output, exist_ok=True)
    os.makedirs(os.path.dirname(dir_batch_scripts), exist_ok=True)
    
    if os.path.exists(path_sbatch): os.remove(path_sbatch)
    
    jobname = re.sub(r"_cfDNA.*", "", cfdna_sample)
    # downsampling_depth = int(sys.argv[1])

    with open(path_sbatch, 'a') as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --job-name={downsampling_depth}X-{jobname}\n')
        file.write('#SBATCH --cpus-per-task=24\n')
        file.write('#SBATCH --mem=100G\n')
        file.write('#SBATCH --time=UNLIMITED\n')
        file.write(f'#SBATCH --error {path_logs}\n')
        file.write(f'#SBATCH --output {path_logs}\n')
        file.write('\n')
        file.write('source /home/amunzur/anaconda3/etc/profile.d/conda.sh\n')
        file.write('conda activate bcftools\n')
        file.write('\n')
        file.write('# Subset the VCF file to the cfDNA sample only. AA requires single sample VCFs.\n')
        file.write(f'bcftools view -s {cfdna_sample} {path_gripss_vcf} > {PATH_gripss_vcf_cfDNA_only}\n')
        file.write('\n')
        file.write('# Run AA\n')
        file.write('conda activate ampliconArchitect\n')
        file.write(f'/home/amunzur/AmpliconSuite-pipeline/AmpliconSuite-pipeline.py \
    -s {cfdna_sample} \
    -t 24 \
    --sorted_bam {path_cfdna_bam} \
    --normal_bam {path_wbc_bam} \
    --ref "GRCh38" \
    --output_directory {dir_AA_output} \
    --run_AA \
    --run_AC \
    --downsample {downsampling_depth} \
    --cnv_bed {path_segmented_cn} \
    --sv_vcf {PATH_gripss_vcf_cfDNA_only}')
        
    print(f"BATCH: {path_sbatch}")

if __name__ == "__main__":
    main()