#!/bin/python

import pandas as pd
import numpy as np
import os
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from itertools import cycle
# from scipy.stats import ttest_ind
import seaborn as sns
import matplotlib.gridspec as gridspec
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu
# import statsmodels.api as sm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import argparse


mpl.rcParams['font.size'] = 10
mpl.rcParams['text.color'] = 'k'
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.handletextpad'] = '0.8'
mpl.rcParams['legend.labelspacing'] = '0.4'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['axes.labelsize'] = 10

source_functions = "/groups/wyattgrp/users/amunzur/toolkit/plotting_UTILITIES.py"
with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

"""
python /groups/wyattgrp/users/amunzur/toolkit/Make_basic_plots_one_group.py \
    --path_muts /groups/wyattgrp/users/amunzur/ironman_ch/results/variant_calling/CHIP_SSCS2_curated.csv \
    --vafcolname VAF \
    --path_sample_information /groups/wyattgrp/users/amunzur/ironman_ch/resources/sample_lists/sample_information.tsv \
    --dir_figures /groups/wyattgrp/users/amunzur/ironman_ch/results/figures \
    --barcolor deepskyblue \
    --keyword IRONMAN
"""

def main():
    parser = argparse.ArgumentParser(description="Make exploratory plots for CH ina given cohort.")
    parser.add_argument("--path_muts", required=True, help="Path to curated mutation list.")
    parser.add_argument("--vafcolname", help="The name of the column in the curated mutation list that has the VAF information.")
    parser.add_argument("--path_sample_information", required=True, help="Path to sample information. Can accomodate with or without header.")
    parser.add_argument("--dir_figures", required=True, help="DIR where the figures will be saved.")
    parser.add_argument("--barcolor", default="deepskyblue", help="The color of the bars (default: deepskyblue).")
    parser.add_argument("--keyword", required=True, help="Keyword, usually will be the name of the cohort, will be used in figure names.")
    args = parser.parse_args()
    
    path_muts=args.path_muts
    vafcolname=args.vafcolname
    path_sample_information=args.path_sample_information
    dir_figures=args.dir_figures
    barcolor=args.barcolor
    keyword=args.keyword
    
    muts=pd.read_csv(path_muts)
    if muts[vafcolname].min()<0.25:
        muts[vafcolname]=muts[vafcolname]*100
    if "VAF_n" in muts.columns and muts["VAF_n"].min()<0.25:
        muts["VAF_n"]=muts["VAF_n"]*100
    if "VAF_t" in muts.columns and muts["VAF_t"].min()<0.25:
        muts["VAF_t"]=muts["VAF_t"]*100
    
    sample_info = pd.read_csv(path_sample_information, sep = "\t")
    ntotal_pts = sample_info.shape[0]
    
    fig = plt.figure(figsize=(8, 9))
    fig.text(0.5, 0.95, "Baseline CH characteristics", ha='center', fontsize=12, fontweight='bold')
    outer_gs = gridspec.GridSpec(3, 1, height_ratios=[0.1, 1.7, 1], hspace = 0, wspace = 0) # outer most gs with 3 rows
    inner_gs1 = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=outer_gs[1], wspace=0.15, hspace = 0.15)
    inner_gs3 = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=outer_gs[2], wspace=0.35, hspace = 0.3)
    inner_inner_left_gs1 = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1, 1], subplot_spec=inner_gs1[0], wspace=0.5, hspace = 0.5)
    inner_inner_left_top_gs1 = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=inner_inner_left_gs1[0], wspace=0.5, hspace = 0.3)
    
    # 1. VAF CORRELATION
    if "VAF_t" in muts.columns and "VAF_n" in muts.columns:
        ax2 = plt.subplot(inner_inner_left_top_gs1[1])
        ax2 = plot_vaf_scatter(muts, ax2, annotate_genes = False, vafcolname=vafcolname)
        ax2.set_xlim((0, 100))
        ax2.set_ylim((0, 100))
        ax2.set_xticks([0, 25, 50, 75, 100])
        ax2.set_yticks([0, 25, 50, 75, 100])
        ax2.set_xticklabels([0, 25, 50, 75, 100])
        ax2.set_yticklabels([0, 25, 50, 75, 100])
        outer_gs.tight_layout(fig)
    # fig.savefig("/groups/wyattgrp/users/amunzur/lu_chip/results/figures/exploration/therap_figure1.png")
    
    # 2. CH PRESENCE ABSENCE PERCENTAGE
    ax0 = plt.subplot(inner_inner_left_top_gs1[0])
    ax0 = plot_ch_presence_absence_bars_UNGROUPED(muts, ax0, ntotal_pts=ntotal_pts, vafcolname=vafcolname, barcolor=barcolor)
    
    # 3. NUMBER OF MUTATIONS
    ax3 = plt.subplot(inner_inner_left_gs1[1])
    ax3, ax3_twin=plot_per_patient_counts_grouped_bar_chart_UNGROUPED(muts=muts, ntotal_pts=ntotal_pts, ax=ax3, barcolor=barcolor, vafcolname=vafcolname)
    
    # 4. GENE COUNTS 
    gs_genes = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1, 1], subplot_spec=inner_gs1[1], wspace=0.4)
    ax_genes_main = plt.subplot(gs_genes[0])
    ax_genes_main.set_yticks([])
    ax_genes_main.set_yticklabels([])
    ax_genes_main.spines[["top", "left", "right", "bottom"]].set_visible(False)
    ax_genes_main.tick_params(axis='x', labeltop = False)
    ax_genes_main.tick_params(axis='y', labelleft=False, left = False, right = False, labelright = False)
    ax_genes = ax_genes_main.twinx()
    ax_swarm = plt.subplot(gs_genes[1], sharey = ax_genes)
    
    # gene_list = ["FLT3", "TERT", 'SH2B3', 'BRCC3', 'STAG2', 'KDM6A', 'KMT2D', 'GNAS', 'CHEK2', 'ATM', 'TP53', 'ASXL1', 'PPM1D', 'TET2', 'DNMT3A']
    gene_list = ["FLT3", "TERT", 'SH2B3', 'SRSF2', 'SF3B1', 'BRCC3', 'STAG2', 'KDM6A', 'KMT2D', 'GNAS', 'CHEK2', 'ATM', 'TP53', 'ASXL1', 'PPM1D', 'TET2', 'DNMT3A']
    
    ax_genes, gene_order=plot_gene_counts_UNGROUPED_single_ax(muts=muts, ax=ax_genes, gene_list = gene_list, ntotal_pts = ntotal_pts, use_entire_denom = True, barcolor=barcolor)
    ax_swarm=plot_vafs_swarm_UNGROUPED(muts, ax_swarm, vafcolname, gene_list = gene_list, gene_order = gene_order)
    
    outer_gs.tight_layout(fig)
    fig.savefig(os.path.join(dir_figures, f"{keyword}_CH_figures.png"))
    fig.savefig(os.path.join(dir_figures, f"{keyword}_CH_figures.pdf"))
    
    # Number of mutations per patient per gene
    fig, ax = plt.subplots(figsize = (4, 4))
    lighter_blue=(0/255, 191/255, 255/255, 0.2)
    ax = plot_mutcounts_per_gene(muts, ntotal_pts, gene_list=gene_order["Gene"], ax=ax, color_dict={1: lighter_blue, 2: "deepskyblue"})
    ax.set_title(f"Mutation multiplicity per gene in {keyword}")
    fig.tight_layout()
    fig.savefig(os.path.join(dir_figures, "baseline_gene_counts_nmuts_per_patient.png"))
    fig.savefig(os.path.join(dir_figures, "baseline_gene_counts_nmuts_per_patient.pdf"))

if __name__ == "__main__":
    main()

# path_muts="/groups/wyattgrp/users/amunzur/prince_ch/results/variant_calling/CHIP_SSCS2_curated.csv"
# vafcolname="VAF_n"
# path_sample_information="/groups/wyattgrp/users/amunzur/prince_ch/resources/sample_lists/sample_information.tsv"
# dir_figures="/groups/wyattgrp/users/amunzur/prince_ch/results/figures"
# barcolor="deepskyblue"
# keyword="PRINCE"