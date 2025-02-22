#!/bin/python
import os
import pandas as pd
import numpy as np
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.font_manager import FontProperties
from matplotlib.colors import Normalize, to_rgba
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import seaborn as sns
import matplotlib.colors as mcolors
import argparse

mut_dict = {
    "missense": '#79B443',
    "stopgain": '#BD4398',
    "frameshift indel": '#FFC907',
    "nonframeshift indel": '#a9a9a9',
    "splicing": "darkorange", 
    "startloss": "red", 
    "stoploss": "deepskyblue"}

mpl.rcParams['font.size'] = 8
mpl.rcParams['text.color'] = 'k'
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.handletextpad'] = '0.8'
mpl.rcParams['legend.labelspacing'] = '0.4'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 7
mpl.rcParams['axes.labelsize'] = 10
mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2

mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['legend.handletextpad'] = '0.7'
mpl.rcParams['legend.labelspacing'] = '0.2'
plt.rcParams['legend.handlelength'] = 0.4
plt.rcParams['legend.handleheight'] = 0.70

# # LOAD DATASETS
# DIR_working = "/groups/wyattgrp/users/amunzur/ironman_ch"
# path_muts="/groups/wyattgrp/users/amunzur/ironman_ch/results/variant_calling/IRONMAN_CH_calls.csv"
# PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")

# source functions
source_functions = "/groups/wyattgrp/users/amunzur/toolkit/UTILITIES_make_OP.py"
with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

def main():
    parser = argparse.ArgumentParser(description="Makes oncoprint.")
    parser.add_argument("--path_muts", required=True, help="Path to curated mutation list.")
    parser.add_argument("--dir_figures", required=True, help="DIR where the figures will be saved.")
    args = parser.parse_args()
    
    path_muts=args.path_muts
    dir_figures=args.dir_figures
    
    # CH mutations
    chip_main = pd.read_csv(path_muts)
    chip_main["Timepoint"]="Baseline"
    chip_main["Status"]="CHIP"
    # chip_main=chip_main[chip_main["Function"]!="upstream"]

    if "VAF_t" in chip_main.columns and chip_main["VAF_t"].min()<0.25:
        chip_main["VAF_t"]=chip_main["VAF_t"]*100
    
    if "VAF_n" in chip_main.columns and chip_main["VAF_n"].min()<0.25:
        chip_main["VAF_n"]=chip_main["VAF_n"]*100
    
    if "VAF" in chip_main.columns and chip_main["VAF"].min()<0.25:
        chip_main["VAF"]=chip_main["VAF"]*100

    # chip_main=chip_main.rename(columns={"VAF": "VAF_t"})
    # chip_main.loc[chip_main["Gene"]=="TERT", "Gene"]="TERT p."

    # all_vars_chip =co all_vars_chip[(all_vars_chip["Dependent"] == False) & (all_vars_chip["Timepoint"] == "Baseline")]
    # chip_main["Consequence"] = chip_main["Consequence"].replace({"Frameshift deletion":"Frameshift indel", "Frameshift insertion":"Frameshift indel", "frameshift_insertion": "Frameshift indel", "Nonframeshift deletion":"Nonframeshift indel", "Nonframeshift insertion":"Nonframeshift indel"})
    chip_main["Consequence"] = chip_main["Consequence"].str.lower()
    chip_main["Consequence"] = chip_main["Consequence"].replace(
        {"frameshift_deletion":"frameshift indel", 
         "frameshift insertion":"frameshift indel", 
         "nonframeshift deletion":"nonframeshift indel", 
         "nonframeshift insertion":"nonframeshift indel", 
         "nonframeshift_deletion": "nonframeshift indel",
         "frameshift_insertion": "nonframeshift indel",
         ".": "splicing", 
         "nonsynonymous_snv": "missense"})
        
    # Compute gene ordering
    path_gene_groups="/groups/wyattgrp/users/amunzur/ironman_ch/resources/gene_groups.tsv"
    genes_to_include = chip_main["Gene"].unique()
    genes_to_include=["DNMT3A", "TET2", "ATM", "CHEK2", "PPM1D", "TP53", "GNAS", "KDM6A", "TERT", "FLT3", "NRAS", "BRCC3", "KMT2D", "MYD88", "GATA2", "CEBPA", "KRAS", "ASXL1", "SETDB1", "RUNX1", "SRSF2", "CUX1", "BCOR", "SH2B3", "BCORL1", "SETD2", "NF1"]
    # gene_order_df=assign_mutation_values_singledf(genes = genes_to_include, muts_df=chip_main, bar_height=1, omit_missing_row = True, path_gene_groups=path_gene_groups)
    gene_order_df=assign_mutation_values_singledf_from_tsv(path_gene_groups)
    ch_plotting_df=calculate_mutation_burden_per_gene(chip_main, samplenamecol="Patient_id").merge(gene_order_df, how = "inner")
    
    # CH dfs
    # ch_plotting_df = generate_plotting_df(chip_main, mut_dict, gene_order_df) # makes up bulk of the plot
    ch_vaf = get_mutation_counts_per_gene(chip_main)[1].merge(gene_order_df[["Gene", "CH position"]], how = "left").rename(columns = {"CH position": "Order on oncoprint"}).set_index("Order on oncoprint").drop("Gene", axis = 1)
    ch_mut_counts  = get_mutation_counts_per_gene(chip_main, make_a_seperate_category_for_multiple=False)[0].merge(gene_order_df[["Gene", "CH position"]], on="Gene", how = "inner").rename(columns = {"CH position": "Order on oncoprint"}).set_index("Order on oncoprint").drop("Gene", axis = 1)
    
    ch_pt_counts = calculate_patient_counts(chip_main)
    
    samples_enumerated = ch_plotting_df.merge(ch_pt_counts)
    samples_enumerated = samples_enumerated[[
        "Patient_id", 
        'Gene',
        'Consequence',
        'Count',
    ]]
    
    samples_enumerated["Gene"] = pd.Categorical(samples_enumerated["Gene"], categories=genes_to_include, ordered=True)
    samples_enumerated["Consequence"] = pd.Categorical(samples_enumerated["Consequence"], categories=mut_dict.keys(), ordered=True)
    samples_enumerated = samples_enumerated.sort_values(
        by=["Gene", "Count"], 
        ascending=[True, False]).reset_index(drop = True).reset_index(drop = True)
    # samples_enumerated = samples_enumerated.sort_values(by="Count", ascending=False).reset_index(drop = True)
    
    samples_enumerated=samples_enumerated.drop_duplicates("Patient_id")["Patient_id"].reset_index(drop = True).reset_index()
    samples_enumerated = samples_enumerated[["Patient_id", "index"]].rename(columns = {"index": "Samples_enumerated"})
    
    ch_plotting_df = ch_plotting_df.merge(samples_enumerated, how = "left")
    ch_plotting_df["Mutation_color"]=ch_plotting_df["Consequence"].map(mut_dict)
    ch_pt_counts = ch_pt_counts.merge(samples_enumerated, how = "left")
    
    # PLOTTING    
    fig_width = 4
    fig_height = 8
    tc_height = 2
    mutcounts_height = 2
    sex_height = 0.5
    age_height = 0.5
    diagnosis_height = 0.5
    main_height = len(genes_to_include)*2
    
    fig = plt.figure(figsize = (fig_width, fig_height))
    gs = mpl.gridspec.GridSpec(ncols = 2, nrows = 9, height_ratios = [tc_height, tc_height, mutcounts_height, sex_height, age_height, diagnosis_height, main_height, 2, 3], width_ratios = [1, 0.1], hspace = 0, wspace = 0)
    
    # Add and space out subplots
    # ax_ctDNA_vaf = fig.add_subplot(gs[0, 0]) #Plotting the highest VAF per patient
    ax_CH_vaf = fig.add_subplot(gs[1, 0]) #Plotting the highest VAF per patient
    ax_mut_counts_patient = fig.add_subplot(gs[2, 0], sharex=ax_CH_vaf) 
    ax_main = fig.add_subplot(gs[6, 0], sharex = ax_mut_counts_patient)
    ax_mut_counts = fig.add_subplot(gs[6, 1], sharey = ax_main) 
    ax_legend = fig.add_subplot(gs[8, 0], sharex = ax_main) 
    
    bar_height = 0.6
    bar_width = 0.6
    
    # ax_ctDNA_vaf = plot_highest_VAF_per_patient(ax_ctDNA_vaf, ctDNA_main, samples_enumerated, bar_width = bar_width)
    ax_CH_vaf = plot_highest_VAF_per_patient(ax_CH_vaf, chip_main, samples_enumerated, bar_width = bar_width, vaf_col="VAF_t")
    ax_main= plot_oncoprint_single_dataset(ax_main, gene_order_df, ch_plotting_df, bar_height, bar_width, colnamepos="CH position", size_scatter=2)
    ax_main = beautify_OP_ax_single_dataset(ax_main, gene_order_df, samples_enumerated, xlabel = "", fontsize = 6, bar_height=bar_height)    
    [ax_mut_counts, ax_secondary] = plot_mut_counts_per_gene(ch_mut_counts, ch_vaf, ax_mut_counts, mut_dict, bar_height = bar_height, logscale = True, vaf_col = "VAF_t", show_median=False) # Plot CH mutations
    
    ax_mut_counts_patient = plot_mut_counts_per_patient_single_df(ax_mut_counts_patient, ch_pt_counts, bar_width)
    ax_legend = add_legend_simple(mut_dict, ax_legend)
    
    fig.savefig(os.path.join(dir_figures, "OP.pdf"), transparent=True)
    fig.savefig(os.path.join(dir_figures, "OP.png"), transparent=True)
    
    print(f"SAVED TO {dir_figures}")

if __name__ == "__main__":
    main()

# path_muts="/groups/wyattgrp/users/amunzur/prince_ch/results/variant_calling/CHIP_SSCS2_curated.csv"

"""
python /groups/wyattgrp/users/amunzur/toolkit/make_OP_show_all_muts.py \
    --path_muts /groups/wyattgrp/users/amunzur/prince_ch/results/variant_calling/CHIP_SSCS2_curated.csv \
    --dir_figures /groups/wyattgrp/users/amunzur/prince_ch/results/figures/patient_profiles \
"""

path_muts="/groups/wyattgrp/users/amunzur/prince_ch/results/variant_calling/CHIP_SSCS2_curated.csv"
dir_figures="/groups/wyattgrp/users/amunzur/prince_ch/results/figures/patient_profiles"

