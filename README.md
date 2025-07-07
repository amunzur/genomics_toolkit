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
