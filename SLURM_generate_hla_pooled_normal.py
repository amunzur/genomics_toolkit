#!/home/amunzur/anaconda3/envs/snakemake/bin/python

"""
Based on HLA types, generates the cnvkit command to generate a pooled normal.
"""
    
import os
import sys
import argparse
import re
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Run AA pipeline.")
    parser.add_argument("--path_compiled_hla_types", required=True, help="Path to the file that has all HLA types from all samples compiled into one file.")
    parser.add_argument("--dir_cnn", required=True, help="Abs path to dir where all bams are.")
    parser.add_argument("--path_hg38", required=True, help="The ref genome the cnn files came from.")
    parser.add_argument("--dir_logs", required=True, help="Where the Slurm error messages and outputs will be saved.")
    parser.add_argument("--dir_cnvkit_pooled_normals", required=True, help="Output files will be saved here for each allele type.")
    parser.add_argument("--dir_batch_scripts", required=True, help="Generated sbatch file will be saved here.")

    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    path_compiled_hla_types=args.path_compiled_hla_types
    dir_cnn=args.dir_cnn
    path_hg38=args.path_hg38
    dir_logs=args.dir_logs
    dir_cnvkit_pooled_normals=args.dir_cnvkit_pooled_normals
    dir_batch_scripts=args.dir_batch_scripts
    
    hla_types=pd.read_csv(path_compiled_hla_types)[["Patient", "LILAC_targeted"]].groupby("LILAC_targeted")
    for i, group in hla_types:
        pts=group["Patient"].unique()
        hla_type = group["LILAC_targeted"].str.replace("*", "_", regex=False).values[0].replace(":", "_")
        cnn_paths = [os.path.join(dir_cnn, f) for pt in pts for f in os.listdir(dir_cnn) if pt in f and "WBC" in f]
        cnn_paths_string=" ".join(cnn_paths)
        path_output=os.path.join(dir_cnvkit_pooled_normals, f"{hla_type}_reference.cnn")
        cnvkit_command=f"cnvkit.py reference {cnn_paths_string} -f {path_hg38} -o {path_output} --male-reference --sample-sex male"
        #
        path_sbatch=os.path.join(dir_batch_scripts, hla_type+".batch")
        path_logs=os.path.join(dir_logs, hla_type+".log")
        os.makedirs(dir_cnvkit_pooled_normals, exist_ok=True)
        os.makedirs(os.path.dirname(dir_batch_scripts), exist_ok=True)
        
        if os.path.exists(path_sbatch): os.remove(path_sbatch)
        
        jobname = f"Pooled_norm_{hla_type}"
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
            file.write('# Run CNVKIT reference\n')
            file.write(cnvkit_command)
            
        print(f"BATCH: {path_sbatch}")

if __name__ == "__main__":
    main()

# /groups/wyattgrp/users/amunzur/toolkit/SLURM_generate_hla_pooled_normal.py \
    # --path_compiled_hla_types /groups/wyattgrp/users/amunzur/hla_pipeline/results/compiled_hla_types.csv \
    # --path_hg38 /groups/wyattgrp/users/jbacon/reference/hg38_emseq/GRCh38_no_alt_plus_hs38d1.fa
    # --dir_cnn /groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/coverage/genes \
    # --dir_logs /groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/cnvkit_pooled_normal \
    # --dir_cnvkit_pooled_normals /groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/pooled_normals
    # --dir_batch_scripts /groups/wyattgrp/users/amunzur/hla_pipeline/workflow/batch_scripts/cnvkit_generate_pooled_normal
