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

## Information on process_variants_steps

### `STEP1_filter_variants.R`
This script will combine all per sample variant tables into one unfiltered csv file. Here is an example on how to run it:
 
```bash
 Rscript STEP1_filter_variants.R \
 --consensus SSCS2 \
 --caller Freebayes \
 --type chip \
 --combine_tumor_and_wbc TRUE \
 --min_alt_reads 5 \
 --forced_rerun FALSE \
 --working_dir /path/to/main/pipeline/folder/results \
 --PATH_sample_list path/to/paired_samples.tsv \
 --PATH_bg /path/to/background_error_rate_file_for_CH.tsv \
 --PATH_bg_ctDNA /path/to/background_error_rate_file_for_ctDNA.tsv
 ```
 
 - `consensus` should be set to `SSCS2`, because the entirety of the pipeline functions based on SSCS2 reads. 
 - For `caller` Mutect2, Vardict or Freebayes can be provided, capitilization is unimportant. Set `type` to `chip` or `ctdna` based on the type of analysis you are performing.
 - If `combine_tumor_and_wbc` is set to `TRUE`, only variants present in both ctDNA (or another source of tumor material) and WBC samples will be retained, when set to `FALSE` all variants will be retained.
 - `min_alt_reads` is the desired minimum number of altered reads carrying the mutation. Mutations with read support less than this number will be filtered out.
 - `forced_rerun` set to `TRUE` will run the entire script again from scratch to recreate the `CHIP_before_filtering.csv` file. If set to false, it will use the existing `CHIP_before_filtering.csv` to continue filtering. I added this option to save time when experimenting with various ways of filtering.
 - `working_dir` should be set to the `results` directory of the pipeline.
 - `PATH_sample_list` is a two column tab seperated file with WBC and cfDNA as column names. It is important to have ONE pair of samples per row. In case there are multiple ctDNA samples with a single WBC available, that's OK, you can repeat the same WBC name in different rows for each ctDNA sample. When `combine_tumor_and_wbc` is set to true, these pairings will be used. 
 - `PATH_bg` and `PATH_bg_ctDNA` are background error rate files which can be generated using the `bg_error_scripts` in this repository.
 - When running germline variants analysis, if the `exclude_nonpathogenic_germline` flag is passed benign germline variants will not be reported.
 - When the `do_not_impose_vaf_upper_limit` flag is passed upper limit for calling CH variants will be removed. This is useful to detect progression CH variants that may have any VAF, including in the germline thresholds.
 - Using `keyword_for_output_files` you can add keywords to the output files, this is useful when applying the script to different scenarious.

This script will create the following three output files. Variants in `min_alt_reads_5_CHIP_after_filtering.csv` will be filtered to retain only those present in both ctDNA and WBC samples, saved to `min_alt_reads_5_CHIP_final.csv`. 

```bash
CHIP_before_filtering.csv
min_alt_reads_5_CHIP_after_filtering.csv
min_alt_reads_5_CHIP_final.csv

```

### `STEP2_combine_variant_callers.R`
Combines variants from the three callers into one file, keeping variants detected by 2/3 callers.

```bash
Rscript STEP2_combine_variant_callers.R \
--dir_working /path/to/main/pipeline/folder \
--mutation_type chip \
--DIR_bams /path/to/main/pipeline/folder/results/data/bam/sorted \
--path_output /path/ch_mutations.csv \
--WBC_only \
--input_file_keyword WBConly
```
- For `dir_working` provide the main working directory for the pipeline. Don't provide the results folder.
- `chip` or `somatic` can ba provided for `mutation_type`.
- `DIR_bams`: Absolute path to where all the aligned bams are located. cfDNA and WBC bams should be both here.
- `path_output` is the output file where the mutations will be saved.
- `WBC_only` pass this flag if your data is WBC only with no cfDNA. If both cfDNA and WBC are available omit the flag.
- `input_file_keyword` should be the same as the keyword used for `STEP1_filter_variants.R`.

### `STEP3_make_IGV_snapshots_tumor_wbc.py`
Generates IGV screenshots of all variants in the file generated by `STEP2_combine_variant_callers.R`. User then manually needs to manually inspect the PNG files. There is no need to manually keep track of which variants to retain or discard, deleting the unwanted variants from the directory is sufficient. 

```bash
pyhton STEP3_make_IGV_snapshots_tumor_wbc.py \
--input_variants /path/to/ch_mutations.csv \
--PATH_batch /path/to/save/igv_batch_script.bat \
--DIR_snapshots /path/to/save/snapshots \
--prefix "some_string" \
--suffix "some_string" \
--given_range 150 \
--WBC_only

```
- `input_variants` is the absolute path to the csv file produced by `STEP2_combine_variant_callers.R`.
- `PATH_batch` is the full path (including filename) for the IGV batch script that will be generated.
- `DIR_snapshots` is the directory where IGV will save all the PNG snapshot files.
- `prefix` and `suffix` are optional. If provided, they will be added to the beginning or end of sample names when constructing paths. If not provided, sample names will be used as-is. This useful to apply this to various csv files that may have the sample names for bam files, but may not have the full path available. `suffix` can be used to append the `.bam` to the sample name.
- `given_range` sets the number of base pairs flanking the variant locus to be visualized in IGV. For example, a value of 150 will display 150bp upstream and downstream of the mutation. `200` usually works well.
- `WBC_only` can be passed if the dataset only contains WBC samples with no matched cfDNA. Leave this flag out if both WBC and cfDNA are present.
- If running in `WBC_only` mode required columns for the `input_variants` are `"Patient_id", "Protein_annotation", "Gene", "Chrom", "Position", "Path_bam_n"`. If running in cfDNA-WBC matched mode required colnames are `"Patient_id", "Protein_annotation", "Gene", "Chrom", "Position", "Sample_name_t", "Path_bam_t", "Path_bam_n"`.


### STEP4_curate_mutations.py

After manually reviewing IGV snapshots produced by `STEP3_make_IGV_snapshots_tumor_wbc.py`, this script detects which mutations were removed (i.e. images not present in the curated screenshots directory) and automatically adds them to a blacklist. Retained mutations are saved to a new file for downstream analyses. This is the final mutations file.

```bash
python STEP4_curate_mutations.py \
--mutation_type somatic \
--DIR_working /path/to/main/pipeline/folder \
--DIR_curated_screenshots /path/to/screenshots \
--do_misc_filtering
```

- `mutation_type` should be set to `chip` or `somatic` depending on the type of variants being curated.
- `DIR_working` is the main working directory of the pipeline (not the results folder).
- `DIR_curated_screenshots` is the folder containing the PNG files you have kept after manual IGV review. Variants corresponding to deleted images will be considered excluded. 
- `do_misc_filtering` is optional. If passed, the script performs additional filtering, removing mutations where the cfDNA:WBC VAF ratio exceeds 20 and excluding mutations in LPAR6.

This script will generate two output files:

```bash
resources/validated_variants/{mutation_type}_to_exclude_IGV.csv
results/variant_calling/{mutation_type}_SSCS2_curated.csv

{mutation_type}_to_exclude_IGV.csv contains variants that were removed after IGV review.
{mutation_type}_SSCS2_curated.csv contains the final curated set of retained variants.
```
