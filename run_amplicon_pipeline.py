#!/usr/bin/env python3

# Writes batch script to run Amplicon architect.

# /groups/wyattgrp/users/amunzur/toolkit/calculate_CN_segment.py \
#     --path_segmentation /groups/wyattgrp/users/amunzur/ecdna_project/igv_tracks/OZM-054-018-003-cfDNA-2016Jan13_segments.tsv \
#     --path_output /groups/wyattgrp/users/amunzur/ecdna_project/igv_tracks/OZM-054-018-003-cfDNA-2016Jan13_segments_CN.tsv \
#     --diploid_level -0.707 \
#     --ctdna_fraction 0.691 

import sys
import os
import pandas as pd
import argparse

def calculate_CN(path_segmentation, diploid_level, ctdna_fraction):
    df = pd.read_csv(path_segmentation, sep = "\t")
    df = df[~df["CHR"].isin(["chrX", "chrY", "chrM"])]
    df["CN"] = df.apply(lambda row:  (2 * (ctdna_fraction + 2**(row["LOGRATIO"] - diploid_level) - 1)) / ctdna_fraction, axis=1)
    df = df[['CHR', 'START', 'END', 'LOGRATIO', 'CN']]
    return(df)

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Calculate raw copy numbers from segmentation file.")
    parser.add_argument("--path_segmentation", required=True, help="Absolute path to the segmentation file")
    parser.add_argument("--path_output", help="Optional, if given saves here.")
    parser.add_argument("--diploid_level", required=True, help="Supply the sample's diploid level")
    parser.add_argument("--ctdna_fraction", required=True, help="Supply the sample's ctdna fraction as an integer, not as a fraction.")
    
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    if args.path_output is None:
        path_output = args.path_segmentation.replace(".tsv", "_CN.tsv")
    else: 
        path_output = args.path_segmentation.replace(".tsv", "_CN.tsv")
    
    # print(f"path_segmentation='{args.path_segmentation}'")
    # print(f"path_output='{args.path_output}'")
    # print(f"diploid_level='{args.diploid_level}'")
    # print(f"ctdna_fraction='{args.ctdna_fraction}'")
    
    if float(args.ctdna_fraction) > 1:
        raise ValueError("Please provide the ctDNA fraction as an integer, not fraction.")
    elif float(args.ctdna_fraction) == 0:
        raise ValueError("Provided ctDNA fraction is 0.")
    
    df = calculate_CN(args.path_segmentation, float(args.diploid_level), float(args.ctdna_fraction))
    df.to_csv(path_output, index = None, header = None, sep = "\t")
    print(f"Output saved to {path_output}")

if __name__ == "__main__":
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
$ cat generate_AA_batch_scripts.py
"""
For a given list of samples, generate per sample batch scripts to run AA.
"""
import pandas as pd
import os

DIR_working = "/groups/wyattgrp/users/amunzur/bladder_ecdna"
PATH_paired_samples = "/groups/wyattgrp/users/amunzur/bladder_ecdna/sample_lists/paired_samples.tsv"
DIR_bams = "/groups/wyattgrp/users/amunzur/bladder_ecdna/alignments_sorted"
slurm_dir = "/groups/wyattgrp/users/amunzur/bladder_ecdna/logs/AA_no_cnvkit"

samps = pd.read_csv(PATH_paired_samples, sep = "\t")

for i, row in samps.iterrows():
    test = row["TEST"]
    ref = row["REF"]
    jobname = test.split("_")[0].replace("GU-", "")
    PATH_batch = os.path.join(DIR_working, "batch_scripts/run_AA", test + ".bash")
    PATH_slurm_error = os.path.join(slurm_dir, test + ".err")
    PATH_slurm_out = os.path.join(slurm_dir, test + ".out")
    PATH_cfDNA_bam = os.path.join(DIR_working, "alignments_sorted", test + ".bam")
    PATH_wbc_bam = os.path.join(DIR_working, "alignments_sorted", ref + ".bam")
    AA_output = os.path.join(DIR_working, "results/AA_no_cnvkit", test)
    #AA_output = os.path.join("/groups/wyattgrp/users/asli_test/AA_no_cnvkit", test)
    PATH_gripss_vcf = os.path.join(DIR_working, f"results/gripss/{test}/{test}.gripss.filtered.vcf.gz")
    PATH_gripss_vcf_cfDNA_only = os.path.join(DIR_working, f"results/gripss/{test}/{test}.gripss.filtered_cfDNA_only.vcf")
    PATH_cn = os.path.join(DIR_working, f"igv_tracks_1500bp/{test}_segments.tsv")
    if os.path.exists(PATH_batch):
        os.remove(PATH_batch)
    with open(PATH_batch, 'a') as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --job-name=W-{jobname}\n')
        file.write('#SBATCH --cpus-per-task=64\n')
        file.write('#SBATCH --mem=512G\n')
        file.write('#SBATCH --time=30-00:00:00\n')
        file.write(f'#SBATCH --error {PATH_slurm_error}\n')
        file.write(f'#SBATCH --output {PATH_slurm_out}\n')
        file.write('\n')
        file.write('source /home/amunzur/anaconda3/etc/profile.d/conda.sh\n')
        file.write('conda activate bcftools\n')
        file.write('\n')
        file.write('# Subset the VCF file to the cfDNA sample only. AA requires single sample VCFs.\n')
        file.write(f'bcftools view -s {test} {PATH_gripss_vcf} > {PATH_gripss_vcf_cfDNA_only}\n')
        file.write('\n')
        file.write('# Run AA\n')
        file.write('conda activate ampliconArchitect\n')
        file.write(f'/home/amunzur/AmpliconSuite-pipeline/AmpliconSuite-pipeline.py') \






-s {test} \
-t 32 \
--sorted_bam {PATH_cfDNA_bam} \
--normal_bam {PATH_wbc_bam} \
--ref "GRCh38" \
--output_directory {AA_output} \
--run_AA \
--run_AC \
--downsample -1 \
--cnv_bed {PATH_cn} \
--sv_vcf {PATH_gripss_vcf_cfDNA_only}')
