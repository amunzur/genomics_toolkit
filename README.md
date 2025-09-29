# Genomic Analysis Toolkit

This repository contains a collection of scripts and tools I developed over time to assist with various aspects of genomic data analysis, including copy number analysis, ecDNA characterization, HLA typing, variant processing, and custom visualization tools.

## Directory Structure

- **copy_number/**  
  Scripts for copy number variation analysis and processing, including CNVkit workflows and pooled normal generation.

- **ecDNA/**  
  Pipeline scripts for analyzing extrachromosomal DNA (ecDNA) events, including GridSS and AmpliconArchitect integration.

- **hla/**  
  Tools and batch scripts for HLA typing and related analyses (e.g., DASH, OptiType, LOHHLA), as well as GC correction and depth metrics.

- **metrics/**  
  Utility scripts for computing coverage metrics such as average sequencing depth.

- **process_variants_steps/**  
  Stepwise variant processing and curation scripts, including filtering, combining variant callers, mutation curation, and IGV snapshot generation. Includes helpful utilities for variant data management.
  
 `STEP1_filter_variants.R`: This script will combine all per sample variant tables into one unfiltered csv file. Here is an example on how to run it:
  Rscript STEP1_filter_variants.R \
  --consensus SSCS2 \
  --caller Freebayes \
  --type chip \
  --combine_tumor_and_wbc TRUE \
  --min_alt_reads 5 \
  --forced_rerun FALSE \
  --working_dir /groups/wyattgrp/users/amunzur/rumble/Clonal-hematopoiesis-pipeline/results \
  --PATH_sample_list path/to/paired_samples.tsv \
  --PATH_bg /path/to/background_error_rate_file_for_CH.tsv \
  --PATH_bg_ctDNA /path/to/background_error_rate_file_for_ctDNA.tsv
  
  `consensus` should be set to `SSCS2`, because the entirety of the pipeline functions based on SSCS2 reads. For `caller` Mutect2, Vardict or Freebayes can be provided, capitilization is unimportant. Set `type` to `chip` or `ctdna` based on the type of analysis you are performing. If `combine_tumor_and_wbc` is set to `TRUE`, only variants present in both ctDNA (or another source of tumor material) and WBC samples will be retained, when set to `FALSE` all variants will be retained. `min_alt_reads` is the desired minimum number of altered reads carrying the mutation. Mutations with read support less than this number will be filtered out. `forced_rerun` set to `TRUE` will run the entire script again from scratch to recreate the `CHIP_before_filtering.csv` file. If set to false, it will use the existing `CHIP_before_filtering.csv` to continue filtering. I added this option to save time when experimenting with various ways of filtering. `working_dir` should be set to the `results` directory of the pipeline. `PATH_sample_list` is a two column tab seperated file with WBC and cfDNA as column names. It is important to have ONE pair of samples per row. In case there are multiple ctDNA samples with a single WBC available, that's OK, you can repeat the same WBC name in different rows for each ctDNA sample. When `combine_tumor_and_wbc` is set to true, these pairings will be used. `PATH_bg` and `PATH_bg_ctDNA` are background error rate files which can be generated using the `bg_error_scripts` in this repository. When running germline variants analysis, if the `exclude_nonpathogenic_germline` flag is passed benign germline variants will not be reported. When the `do_not_impose_vaf_upper_limit` flag is passed upper limit for calling CH variants will be removed. This is useful to detect progression CH variants that may have any VAF, including in the germline thresholds. Using the `keyword_for_output_files` you can add keywords to the output files, this is useful when applying the script to different scenarious.

- **variant_calling/**  
  Batch scripts and wrappers for variant calling tools like FreeBayes, Mutect2, and VarDict.

- **visualization/**  
  Scripts for generating various genomic data visualizations, including oncoprints, lollipop plots, mutation heatmaps, and patient-specific mutation plots.

## Usage

Each directory contains scripts primarily intended to be run via batch systems (e.g., SLURM). Most scripts are designed to be modular and reusable for different datasets and projects.

Example usage snippets can be found as comments inside individual scripts.

## Requirements

- Python (>=3.6), with packages: pandas, numpy, matplotlib, seaborn, argparse  
- R (for variant filtering scripts)  
- External tools as required: CNVkit, GridSS, AmpliconArchitect, FreeBayes, Mutect2, VarDict, OptiType, LOHHLA, DASH, etc.

## Repository file tree

```bash
├── copy_number
│   ├── calculate_CN_segment.py
│   ├── SLURM_generate_pooled_normal.py
│   ├── SLURM_run_cnvkit_coverage.py
│   └── SLURM_run_cnvkit_fix.py
├── define_classes.py
├── ecDNA
│   ├── run_amplicon_pipeline.py
│   ├── SLURM_run_AA.py
│   └── SLURM_run_gridss.py
├── hla
│   ├── compute_background_profile_GC_correction.bash
│   ├── SLURM_HLA_generate_pooled_normal.py
│   ├── SLURM_run_DASH.py
│   ├── SLURM_run_lilac.py
│   ├── SLURM_run_lohhla.py
│   └── SLURM_run_optitype.py
├── metrics
│   └── compute_average_depth.bash
├── process_variants_steps
│   ├── combine_variant_callers_UTILITIES.R
│   ├── generate_sample_matrix.py
│   ├── STEP1_filter_variants.R
│   ├── STEP2_combine_variant_callers.R
│   ├── STEP3_make_IGV_snapshots_tumor_wbc.py
│   ├── STEP4_curate_mutations.py
│   └── UTILITIES.R
├── README.md
├── rename_fastq.py
├── variant_calling
│   ├── SLURM_run_freebayes.py
│   ├── SLURM_run_mutect2.py
│   └── SLURM_run_vardict.py
└── visualization
    ├── Make_basic_plots_one_group.py
    ├── make_lolliop_plots.py
    ├── make_oncoprint_functions.py
    ├── make_OP_collapsed.py
    ├── make_OP_show_all_muts.py
    ├── per_patient_plots_JUST_CHIP.py
    ├── plotting_UTILITIES.py
    ├── SUPP_master_OP.py
    ├── UTILITIES_make_OP.py
    ├── utilities_plotting_functions.py
    ├── UTILITIES_therap.py
    └── zoomed_in_gene_cn_plot_cnvkit.py
```
