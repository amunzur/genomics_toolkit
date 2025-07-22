import os
import subprocess
import pandas as pd

DIR_results="/groups/wyattgrp/users/amunzur/hla_pipeline/results"
DIR_resources="/groups/wyattgrp/users/amunzur/hla_pipeline/resources"
PATH_samples_list="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv"
DIR_log="/groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/run_dash"
DIR_batch_scripts="/groups/wyattgrp/users/amunzur/hla_pipeline/workflow/batch_scripts/dash"
path_sample_purity_ploidy="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_ploidy.tsv"

samples_df = pd.read_csv(PATH_samples_list, sep="\t")
tumors_list = samples_df["Tumor_name"].tolist()

def get_wbc_name(tumor_name): 
    samples_df = pd.read_csv(PATH_samples_list, sep="\t")
    mask = samples_df["Tumor_name"] == tumor_name
    wbc_row = samples_df[mask]
    if wbc_row.empty:
        raise ValueError("Could not find WBC sample for tumor '{}'".format(tumor_name))
    elif wbc_row.shape[0] != 1:
        raise ValueError("Tumor '{}' appears multiple times in the samples df".format(tumor_name))
    wbcname = wbc_row["WBC_name"].iloc[0]
    return(wbcname)

def get_WBC_hla_results(tumor_name):
    wbcname = get_wbc_name(tumor_name)
    DIR_results = "/groups/wyattgrp/users/amunzur/hla_pipeline/results"
    PATH_hla_types = os.path.join(DIR_results, "polysolver/hla_types", wbcname, "winners.hla.txt")
    return(PATH_hla_types)

# Given the name of a cfDNA sample, return the path to the combine fq file for the WBC sample.
def wbc_trimmed_combined_fq(tumor_name): 
    wbcname = get_wbc_name(tumor_name)
    DIR_results = "/groups/wyattgrp/users/amunzur/hla_pipeline/results"    
    PATH_wbc_fq = os.path.join(DIR_results, "data/fq/trimmed_combined", wbcname + ".fq.gz")
    return(PATH_wbc_fq)

def get_read_counts(tumor_name, sample_type):
    wbcname = get_wbc_name(tumor_name)
    DIR_results = "/groups/wyattgrp/users/amunzur/hla_pipeline/results"
    DIR_readcounts = os.path.join(DIR_results, "metrics", "trimmed_combined_read_counts")
    if sample_type == "cfDNA":
        path = os.path.join(DIR_readcounts, tumor_name + ".txt")
    else: 
        path = os.path.join(DIR_readcounts, wbcname + ".txt")
    with open(path) as file: 
        lines = [line.strip() for line in file.readlines()] # Modify this line
    readcount = lines[0]
    return(readcount)

def return_ploidy_purity(tumor_name, path_sample_purity_ploidy):
    df=pd.read_csv(path_sample_purity_ploidy, sep="\t")
    sample_subset=df[df["sample"]==tumor_name]
    ploidy=sample_subset["ploidy"].values[0]
    purity=sample_subset["ctdna"].values[0]
    
    return(ploidy, purity)

# tumor or cfDNA
for tumor_name in tumors_list: 
    dir_sequenza = os.path.join(DIR_results + "/sequenza_sln/{}".format(tumor_name, tumor_name))
    path_polysolver_winners = get_WBC_hla_results(tumor_name)
    normal_fastq = wbc_trimmed_combined_fq(tumor_name)
    tumor_fastq = os.path.join(DIR_results + "/data/fq/trimmed_combined/{}.fq.gz".format(tumor_name))
    hla_somatic_mutations = os.path.join(DIR_results + "/polysolver/hla_mutations/{}/hla_mutations".format(tumor_name))
    normal_read_count = get_read_counts(tumor_name, "WBC")
    tumor_read_count = get_read_counts(tumor_name, "cfDNA")
    all_allele_reference = os.path.join(DIR_resources + "/abc_complete.fasta")
    model_filename = os.path.join(DIR_resources + "/training.xgboost_model.2021_05_10.p")
    output_dir = os.path.join(DIR_results + "/dash/{}".format(tumor_name))
    ploidy, purity=return_ploidy_purity(tumor_name, path_sample_purity_ploidy)
    
    dash_cmd = "python /groups/wyattgrp/users/amunzur/hla_pipeline/software/DASH/Algorithm/DASH.manuscript.py " + \
           "--path_polysolver_winners {} ".format(path_polysolver_winners) + \
           "--ploidy {} ".format(ploidy) + \
           "--purity {} ".format(purity) + \
           "--normal_fastq {} ".format(normal_fastq) + \
           "--tumor_fastq {} ".format(tumor_fastq) + \
           "--hla_somatic_mutations {} ".format(hla_somatic_mutations) + \
           "--normal_read_count {} ".format(normal_read_count) + \
           "--tumor_read_count {} ".format(tumor_read_count) + \
           "--all_allele_reference {} ".format(all_allele_reference) + \
           "--model_filename {} ".format(model_filename) + \
           "--output_dir {} ".format(output_dir)
    job_name = "DASH_{}".format(tumor_name)
    # Generate a batch script for each sample
    batch_script_path = os.path.join(DIR_batch_scripts, tumor_name + "_run_DASH.bash")
    # subprocess.call("rm " + batch_script_path)
    with open(batch_script_path, 'w') as batch_file:
        batch_file.write('#!/bin/bash\n')
        batch_file.write('#SBATCH --job-name={}\n'.format(job_name))
        batch_file.write('#SBATCH --output={}/{}.log\n'.format(DIR_log, tumor_name))
        batch_file.write('#SBATCH --error={}/{}.log\n'.format(DIR_log, tumor_name))
        batch_file.write('#SBATCH --time=5:00:00\n')
        batch_file.write('#SBATCH --mem=16G\n')
        batch_file.write('#SBATCH --cpus-per-task=4\n')
        batch_file.write('export MPLCONFIGDIR=/groups/wyattgrp/users/amunzur/matplotlib_cache\n')
        batch_file.write('source /home/amunzur/anaconda3/etc/profile.d/conda.sh\n')
        batch_file.write('conda activate dash\n')
        batch_file.write(dash_cmd)
    slurm_cmd = "printf sbatch {}".format(batch_script_path)
    # print(slurm_cmd)
    # subprocess.call(slurm_cmd, shell=True)





    # slurm_cmd = "printf sbatch --job-name={} --output={} --error={}/{}.log {}".format(job_name, output_dir, DIR_log, tumor_name, dash_cmd)
    # print(slurm_cmd)
    # subprocess.call(slurm_cmd, shell=True)

