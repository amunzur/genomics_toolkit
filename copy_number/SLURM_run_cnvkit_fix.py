#!~/anaconda3/envs/snakemake/bin/python

"""
Given a sample name and path to a pooled normal, runs cnvkit fix on the sample.
"""
    
import os
import sys
import argparse
import re

def main():
    parser = argparse.ArgumentParser(description="Run cnvkit fix.")
    parser.add_argument("--path_cnn", required=True, help="This is the path to the sample's cnn file.")
    parser.add_argument("--path_pooled_reference_normal", required=True, help="Pooled normal")
    parser.add_argument("--dir_wbc_vcf", default=None, help="Will look for a VCF file here for snps. If argument not passed, vcf SNPs won't be used during segmentation.")
    parser.add_argument("--dir_logs", required=True, help="Where the Slurm error messages and outputs will be saved.")
    parser.add_argument("--dir_batch_scripts", required=True, help="Generated sbatch file will be saved here.")
    parser.add_argument("--dir_output", required=True, help="An output file will be generated here.")

    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    path_cnn=args.path_cnn
    path_pooled_reference_normal=args.path_pooled_reference_normal
    dir_wbc_vcf=args.dir_wbc_vcf
    dir_logs=args.dir_logs
    dir_batch_scripts=args.dir_batch_scripts
    dir_output=args.dir_output
    
    sample_name=os.path.basename(path_cnn).replace(".cnn", "")
    patient_id=re.sub(r'(_cfDNA|_WBC|_FiT).*', '', sample_name)
    
    if dir_wbc_vcf is not None:
        path_vcf_list=[os.path.join(dir_wbc_vcf, f) for f in os.listdir(dir_wbc_vcf) if patient_id in f and "WBC" in f]
        if len(path_vcf_list)==1:
            path_vcf=path_vcf_list[0]
        else:
            path_vcf=None
    else:
        path_vcf=None
    
    path_output=os.path.join(dir_output, sample_name+".cnn.fix")
    path_output_segment=os.path.join(dir_output, sample_name+".cnn.fix.segment")
    path_sbatch=os.path.join(dir_batch_scripts, sample_name+".batch")
    path_logs=os.path.join(dir_logs, sample_name+".log")
    
    command_fix=f"cnvkit.py fix {path_cnn} /groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/blank {path_pooled_reference_normal} -o {path_output}"
    
    if path_vcf is not None:
        command_segment=f"cnvkit.py segment {path_output} -o {path_output_segment} --vcf {path_vcf}"
    else:
        command_segment=f"cnvkit.py segment {path_output} -o {path_output_segment}"
    
    command_plotting_scatter=f"cnvkit.py scatter {path_output} -s {path_output_segment} -o /groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/figures/{sample_name}_scatter.pdf"
    command_plotting_diagram=f"cnvkit.py diagram {path_output} -s {path_output_segment} -o /groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/figures/{sample_name}_diagram.pdf"
    
    os.makedirs(dir_batch_scripts, exist_ok=True)
    os.makedirs(dir_logs, exist_ok=True)
    os.makedirs(dir_output, exist_ok=True)
    
    if os.path.exists(path_sbatch): os.remove(path_sbatch)
    
    jobname = f"cnvfix_{sample_name}"
    # downsampling_depth = int(sys.argv[1])

    with open(path_sbatch, 'a') as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --job-name={jobname}\n')
        file.write('#SBATCH --cpus-per-task=4\n')
        file.write('#SBATCH --mem=8G\n')
        file.write('#SBATCH --time=30:00\n')
        file.write(f'#SBATCH --error {path_logs}\n')
        file.write(f'#SBATCH --output {path_logs}\n')
        file.write('\n')
        file.write('source /home/amunzur/anaconda3/etc/profile.d/conda.sh\n')
        file.write('conda activate cnv\n')
        file.write('\n')
        file.write('# Run CNVKIT fix\n')
        file.write(f"{command_fix}\n")
        file.write(f"{command_segment}\n")
        file.write(f"{command_plotting_scatter}\n")
        file.write(f"{command_plotting_diagram}\n")
        
    print(f"BATCH: {path_sbatch}")

if __name__ == "__main__":
    main()

# /groups/wyattgrp/users/amunzur/toolkit/SLURM_run_cnvkit_fix.py \
    # --path_cnn /groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/coverage/genes/GU-21-633_FiT_22Oct2021_1.cnn \
    # --path_pooled_reference_normal /groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/pooled_normals/all_wbc_samples_reference.cnn \
    # --dir_logs dir_logs /groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/cnvkit_fix \
    # --dir_batch_scripts /groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/cnvkit_fix \
    # --dir_output /groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/fix
    # --dir_batch_scripts /groups/wyattgrp/users/amunzur/hla_pipeline/workflow/batch_scripts/cnvkit_fix

# while IFS=$'\t' read -r Patient_ID WBC_name Tumor_name; do
#     /groups/wyattgrp/users/amunzur/toolkit/SLURM_run_cnvkit_fix.py \
#         --path_cnn "/groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/coverage/genes/${Tumor_name}.cnn" \
#         --path_pooled_reference_normal /groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/pooled_normals/all_wbc_samples_reference.cnn \
#         --dir_wbc_vcf /groups/wyattgrp/users/amunzur/hla_pipeline/results/vardict/wbc_snps \
#         --dir_logs /groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/cnvkit_fix \
#         --dir_batch_scripts /groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/cnvkit_fix \
#         --dir_output /groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/fix \
#         --dir_batch_scripts /groups/wyattgrp/users/amunzur/hla_pipeline/workflow/batch_scripts/cnvkit_fix
# done < /groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv
