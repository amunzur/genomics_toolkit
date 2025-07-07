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

