#!/home/amunzur/anaconda3/envs/snakemake/bin/python
"""
For a given list of samples, generate per sample batch scripts to run DASH. This is single sample mode.
"""
DIR_hla_pipeline="/groups/wyattgrp/users/amunzur/hla_pipeline"
/groups/wyattgrp/users/amunzur/toolkit/SLURM_run_DASH.py \
--dir_sequenza ${DIR_hla_pipeline}/results/sequenza_sln/AE-092_cfDNA_Baseline_2015Sep28 \
--path_polysolver_winners ${DIR_hla_pipeline}/results/polysolver/hla_types/AE-092_WBC_Baseline_2015Sep28 \
--normal_fastq ${DIR_hla_pipeline}/results/data/fq/trimmed/AE-092_cfDNA_Baseline_2015Sep28.fq.gz \
--tumor_fastq ${DIR_hla_pipeline}/results/data/fq/trimmed_combined/AE-092_WBC_Baseline_2015Sep28.fq.gz \
--hla_somatic_mutations ${DIR_hla_pipeline}/results/polysolver/hla_mutations/AE-092_cfDNA_Baseline_2015Sep28/hla_mutations \
--normal_read_count ${DIR_hla_pipeline}/results/metrics/trimmed_combined_read_counts/AE-092_WBC_Baseline_2015Sep28.txt \
--tumor_read_count ${DIR_hla_pipeline}/results/metrics/trimmed_combined_read_counts/AE-092_cfDNA_Baseline_2015Sep28.txt \
--all_allele_reference /groups/wyattgrp/users/amunzur/gillian_proj/hla-polysolver/data/abc_complete.fasta \
--model_filename ${DIR_hla_pipeline}/resources/training.xgboost_model.2021_05_10.p \
--output_dir ${DIR_hla_pipeline}/results/dash/AE-092_cfDNA_Baseline_2015Sep28


import os
import sys
import argparse




def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Run AA pipeline.")
    parser.add_argument("--dir_sequenza", required=True, help="DIR to the sample specific sequenza_sln dir.")
    parser.add_argument("--path_polysolver_winners", required=True, help="Path to the Polysolver winnes.txt file.")
    parser.add_argument("--normal_fastq", required=True, help="Path to the combined FastQ for the WBC sample.")
    parser.add_argument("--tumor_fastq", required=True, help="Path to the combined FastQ for the WBC sample")
    parser.add_argument("--hla_somatic_mutations", required=True, help="DIR to HLA mutations fr")
    parser.add_argument("--normal_read_count", required=True, help="Outputted VCF will be saved here.")
    parser.add_argument("--tumor_read_count", required=True, help="Batch script will be written here.")
    parser.add_argument("--all_allele_reference", required=True, help="Vardict logs will be here.")
    parser.add_argument("--model_filename", required=True, help="Batch script will be written here.")
    parser.add_argument("--output_dir", required=True, help="Vardict logs will be here.")

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
        file.write('#SBATCH --time=30:00\n')
        file.write(f'#SBATCH --error {path_logs}\n')
        file.write(f'#SBATCH --output {path_logs}\n')
        file.write('\n')
        file.write('# Run Vardict\n')
        file.write(f'/home/jbacon/mambaforge/envs/pipeline/bin/vardict-java \
                   -G {path_hg38} \
                    -f {threshold_min_vaf} \
                    -N {sample_name} \
                    -r {min_alt_reads} \
                    -b {path_bam} \
                    -k 0 -c 1 -S 2 -E 3 -g 4 {path_bed} | \
                    /home/amunzur/VarDictJava/build/install/VarDict/bin/teststrandbias.R | \
                    /home/amunzur/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl > {path_output}')
        print(f"sbatch {path_sbatch}")

if __name__ == "__main__":
    main()

