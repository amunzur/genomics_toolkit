import os
import subprocess
import pandas as pd
import argparse
import sys

"""
SLURM_run_DASH.py

This script automates the generation of SLURM batch scripts to run the DASH pipeline. 
My modifications to the original DASH algorith are here: https://github.com/amunzur/DASH
It prepares and formats necessary inputs such as read counts, purity/ploidy estimates, and HLA types, and then constructs 
a command to execute the DASH algorithm. A sample-specific SLURM script is created and printed for submission.

Inputs required:
- Normal (WBC) and tumor FastQ files (combined R1+R2).
- Output directory (non-sample specific).
- Directory for writing generated batch scripts.

Major components:
- Read count fetching for tumor and normal samples.
- Ploidy and purity lookup from a TSV.
- Reformatting of LILAC HLA types to Polysolver-compatible format.
- Creation of a SLURM batch script to run DASH for each sample.

Usage example:
python SLURM_run_DASH.py \
    --normal_fastq /path/to/normal.fq.gz \
    --tumor_fastq /path/to/tumor.fq.gz \
    --output_dir /path/to/output/dir \
    --DIR_batch_scripts /path/to/batch_scripts/
"""

def get_read_counts(sample_name):
    """
    Returns number of reads present in the given FastQ file.
    """
    dir_read_counts="/groups/wyattgrp/users/amunzur/hla_project/snakemake_hla/results/metrics/trimmed_combined_read_counts"
    path_nreads=os.path.join(dir_read_counts, sample_name+".txt")
    
    if not os.path.exists(path_nreads):
        raise FileNotFoundError(f"Missing read count file: {path_nreads}")
    else:
        with open(path_nreads) as file:
            lines = [line.rstrip() for line in file]
    
    nreads=lines[0]
    return(nreads)

def return_ploidy_purity(tumor_name, path_sample_purity_ploidy):
    df=pd.read_csv(path_sample_purity_ploidy, sep="\t")
    sample_subset=df[df["sample"]==tumor_name]
    if not df.empty:
        ploidy=sample_subset["ploidy"].values[0]
        purity=sample_subset["ctdna"].values[0]
        return(ploidy, purity)
    else:
        raise ValueError(f"Sample not found in the {path_sample_purity_ploidy}.")

def reformat_lilac_hla_types(n_name):
    """
    Reformats LILAC's HLA calls to match Polysolver's format. Dash expects HLA calls to be provided in Polysolver's format.
    n_name is WBC name. No file extensions.
    """
    
    dir_lilac_hla_calls="/groups/wyattgrp/users/amunzur/hla_project/data/lilac"
    path_lilac_calls=os.path.join(dir_lilac_hla_calls, n_name, n_name+".lilac.tsv")
    
    lilac_calls=pd.read_csv(path_lilac_calls, sep="\t")
    alleles=lilac_calls["Allele"]
    
    # Group by locus
    grouped = alleles.groupby(alleles.str[0])
    
    # Format output
    output_lines = []
    for gene, group in grouped:
        hla_gene = f"HLA-{gene}"
        converted = group.str.lower().str.replace("*", "_").str.replace(":", "_").apply(lambda x: f"hla_{x}")
        
        if len(converted) < 2:
            raise ValueError(f"Less than 2 alleles found for {hla_gene} in {path_lilac_calls}")
        
        output_lines.append(f"{hla_gene}\t{converted.iloc[0]}\t{converted.iloc[1]}\n")
    
    path_reformatted_output=path_lilac_calls.replace(".lilac.tsv", "_reformatted.lilac.tsv")
    
    f = open(path_reformatted_output, "w")
    f.writelines(output_lines)
    f.close()
    
    return(path_reformatted_output)

def main():
    parser = argparse.ArgumentParser(description="Run the DASH pipeline for HLA LOH detection.")
    parser.add_argument("--normal_fastq", required=True, help="Path to processed WBC FastQ. Read1 and 2 must be combined into one file.")
    parser.add_argument("--tumor_fastq", required=True, help="Path to processed tumor FastQ. Read1 and 2 must be combined into one file.")
    parser.add_argument("--output_dir", required=True, help="A sample directory will be made inside this dir.")
    parser.add_argument("--DIR_batch_scripts", required=True, help="A sample directory will be made inside this dir.")
    
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    # Generate some of the inputs based on user provided arguments
    t_name=os.path.basename(args.tumor_fastq).replace(".fq.gz", "")
    n_name=os.path.basename(args.normal_fastq).replace(".fq.gz", "")
    
    normal_read_count = get_read_counts(n_name)
    tumor_read_count = get_read_counts(t_name)
    
    all_allele_reference="/groups/wyattgrp/users/amunzur/hla_project/references/abc_complete.fasta"
    model_filename="/groups/wyattgrp/users/amunzur/hla_project/resources/training.xgboost_model.2021_05_10.p"
    sample_output_dir=os.path.join(args.output_dir, t_name)
    path_lilac_hla_calls=reformat_lilac_hla_types(n_name)
    
    dir_hla_mutations=os.path.join(f"/groups/wyattgrp/users/amunzur/hla_pipeline/results/polysolver/hla_mutations/{t_name}/hla_mutations")
    
    path_sample_purity_ploidy="/groups/wyattgrp/users/amunzur/hla_project/sample_ploidy.tsv"
    ploidy, purity=return_ploidy_purity(t_name, path_sample_purity_ploidy)
    
    DIR_log="/groups/wyattgrp/users/amunzur/hla_project/logs/run_dash"
    
    dash_cmd = "python /groups/wyattgrp/users/amunzur/hla_pipeline/software/DASH/Algorithm/DASH.manuscript.py " + \
           "--path_polysolver_winners {} ".format(path_lilac_hla_calls) + \
           "--ploidy {} ".format(ploidy) + \
           "--purity {} ".format(purity) + \
           "--normal_fastq {} ".format(args.normal_fastq) + \
           "--tumor_fastq {} ".format(args.tumor_fastq) + \
           "--hla_somatic_mutations {} ".format(dir_hla_mutations) + \
           "--normal_read_count {} ".format(normal_read_count) + \
           "--tumor_read_count {} ".format(tumor_read_count) + \
           "--all_allele_reference {} ".format(all_allele_reference) + \
           "--model_filename {} ".format(model_filename) + \
           "--output_dir {} ".format(sample_output_dir)
    
    job_name = "DASH_{}".format(t_name)
    
    # Generate a batch script for each sample
    batch_script_path = os.path.join(args.DIR_batch_scripts, t_name + "_run_DASH.bash")
    # subprocess.call("rm " + batch_script_path)
    with open(batch_script_path, 'w') as batch_file:
        batch_file.write('#!/bin/bash\n')
        batch_file.write('#SBATCH --job-name={}\n'.format(job_name))
        batch_file.write('#SBATCH --output={}/{}.log\n'.format(DIR_log, t_name))
        batch_file.write('#SBATCH --error={}/{}.log\n'.format(DIR_log, t_name))
        batch_file.write('#SBATCH --time=5:00:00\n')
        batch_file.write('#SBATCH --mem=16G\n')
        batch_file.write('#SBATCH --cpus-per-task=4\n')
        batch_file.write('export MPLCONFIGDIR=/groups/wyattgrp/users/amunzur/matplotlib_cache\n')
        batch_file.write('source /home/amunzur/anaconda3/etc/profile.d/conda.sh\n')
        batch_file.write('conda activate dash\n')
        batch_file.write(dash_cmd)
    
    slurm_cmd = f"sbatch {batch_script_path}"
    print(slurm_cmd)


if __name__ == "__main__":
    main()

# Run example
# /groups/wyattgrp/users/amunzur/toolkit/hla/SLURM_run_DASH.py \
# --normal_fastq path/to/fastq \
# --tumor_fastq path/to/fastq \
# --output_dir path/to/dir_output \ # Should not be sample specific
# --DIR_batch_scripts path/to/batch \ 

