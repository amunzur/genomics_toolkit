#!/bin/bash
# Computes background profile that can used later on when performing GC correction

# conda activate snakemake

PATH_HLA_genome="/groups/wyattgrp/users/amunzur/ref/abc_complete.2bit"
PATH_HLA_genome_GC_bias_background="/groups/wyattgrp/users/amunzur/ref/abc_complete_GC_bias_background"
bam_file="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/sorted_hla_filtered/GU-17-314_cfDNA-Baseline-2017Dec12.bam"

computeGCBias_background \
    -b ${bam_file} \
    -g ${PATH_HLA_genome} \
    -i \
    -p 12 \
    --output ${PATH_HLA_genome_GC_bias_background}
