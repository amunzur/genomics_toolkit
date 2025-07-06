#!/home/amunzur/anaconda3/envs/snakemake/bin/python

# Make sure to run in base env
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
from matplotlib.lines import Line2D

mut_dict = {
    "missense": '#79B443',
    "stopgain": '#BD4398',
    "frameshift indel": '#FFC907',
    "nonframeshift indel": '#a9a9a9',
    "splicing": "darkorange", 
    "startloss": "red", 
    "stoploss": "deepskyblue", 
    "multiple": "cyan"}

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

# source functions
source_functions = "UTILITIES_make_OP.py" # in the toolkit dir, too.
with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

def assign_custom_index(path_gene_groups, spacing_within=0.3, spacing_between=0.4):
    """
    Assigns a custom index that increases normally within the same group and adds extra space when transitioning.
    """
    df=pd.read_csv(path_gene_groups, sep="\t", header=None, names=["Gene", "Group"])
    df = df.copy()  # Avoid modifying the original DataFrame
    col_group = df.columns[1]  # Automatically detect group column
    new_index = []
    current_idx = 0
    prev_group = None
    
    previous_group=None
    current_idx=0
    new_index = []
    for i, row in df.iterrows():
        row_group=row["Group"]
        if previous_group is None or row_group==previous_group:
            current_idx+=spacing_within
        else:
            current_idx+=spacing_between
        previous_group=row_group
        new_index.append(current_idx)
    
    df["CH position"] = new_index
    return df

def generated_samples_enumerated_dict(ch_plotting_df, ch_pt_counts):
    """
    
    """   
    samples_enumerated = ch_plotting_df.merge(ch_pt_counts)
    samples_enumerated = samples_enumerated[[
        "Patient_id",
        "Arm" ,
        'Gene',
        'Consequence',
        'Count',
    ]]
    
    genes_to_include=["DNMT3A", "ASXL1", "TET2", "PPM1D", "CHEK2", "ATM", "TP53"]
    gene_sort_order=["PPM1D", "CHEK2", "ATM", "TP53", "DNMT3A", "TET2", "ASXL1"]
    samples_enumerated["Gene"] = pd.Categorical(samples_enumerated["Gene"], categories=gene_sort_order, ordered=True)
    samples_enumerated["Consequence"] = pd.Categorical(samples_enumerated["Consequence"], categories=mut_dict.keys(), ordered=True)
    samples_enumerated=samples_enumerated.dropna()
    samples_enumerated_lu=samples_enumerated[samples_enumerated["Arm"]=="LuPSMA"]
    samples_enumerated_caba=samples_enumerated[samples_enumerated["Arm"]=="Cabazitaxel"]
    
    samples_enumerated_lu=samples_enumerated_lu.sort_values(by=["Gene", "Consequence"], ascending=[True, False]).reset_index(drop = True).reset_index(drop = True)
    samples_enumerated_lu=samples_enumerated_lu.drop_duplicates("Patient_id")[["Patient_id", "Arm"]].reset_index(drop = True).reset_index()
    samples_enumerated_lu=samples_enumerated_lu[["Patient_id", "Arm", "index"]].rename(columns = {"index": "Samples_enumerated"})
    
    samples_enumerated_caba=samples_enumerated_caba.sort_values(by=["Gene", "Consequence"], ascending=[True, False]).reset_index(drop = True).reset_index(drop = True)
    samples_enumerated_caba=samples_enumerated_caba.drop_duplicates("Patient_id")[["Patient_id", "Arm"]].reset_index(drop = True).reset_index()
    samples_enumerated_caba=samples_enumerated_caba[["Patient_id", "Arm", "index"]].rename(columns = {"index": "Samples_enumerated"})
    
    samples_enumerated_dict={"LuPSMA": samples_enumerated_lu, "Cabazitaxel": samples_enumerated_caba}
    
    return(samples_enumerated_dict)

def plot_oncoprint_compressed(ax_main, muts_df, bar_height, bar_width):
    for i, row in muts_df.iterrows():
        ypos=row["CH position"]
        xpos=row["Samples_enumerated"]
        bar_color=row["Mutation_color"]
        ax_main.bar(x=xpos, bottom=ypos, height = bar_height, width=bar_width, capstyle='butt', color=bar_color, edgecolor = "none", linewidth = 0.25, zorder = 10)
    return(ax_main)

path_muts="/path/to/CHIP_progression_plasma_curated.csv"
vaf_colname="VAF%"

# CH mutations
chip_main = pd.read_csv(path_muts)
chip_main["Status"]="CHIP"

if chip_main[vaf_colname].min()<0.25:
    chip_main[vaf_colname]=chip_main[vaf_colname]*100

# Renaming etc.
chip_main=chip_main[chip_main["Consequence"]!="synonymous_SNV"]
chip_main["Consequence"] = chip_main["Consequence"].str.lower()
chip_main["Consequence"] = chip_main["Consequence"].replace(
    {"frameshift_deletion":"frameshift indel", 
     "frameshift insertion":"frameshift indel", 
     "nonframeshift deletion":"nonframeshift indel", 
     "nonframeshift insertion":"nonframeshift indel", 
     "nonframeshift_deletion": "nonframeshift indel",
     "nonframeshift_insertion": "nonframeshift indel",
     "frameshift_insertion": "nonframeshift indel",
     "frameshift_substitution": "frameshift indel",
     "nonframeshift_substitution": "nonframeshift indel",
     ".": "splicing", 
     "nonsynonymous_snv": "missense"})

# Here we collapse the mutations
per_patient_collapsed_mutations_dict = {}
for (pt, gene, arm), group in chip_main.groupby(["Patient_id", "Gene", "Arm"]):
    n_mut_types = group["Consequence"].nunique()
    per_patient_collapsed_mutations_dict[(pt, arm, gene)] = "multiple" if n_mut_types > 1 else group["Consequence"].iloc[0]

ch_plotting_df = pd.DataFrame([(pt, arm, gene, consequence) for (pt, arm, gene), consequence in per_patient_collapsed_mutations_dict.items()], columns=["Patient_id", "Arm", "Gene", "Consequence"])

# Compute gene ordering
path_gene_groups="/groups/wyattgrp/users/amunzur/lu_chip/resources/clinical_tables/gene_groups_op.tsv"
gene_pos_df=assign_custom_index(path_gene_groups)
chip_main=chip_main.merge(gene_pos_df, how="inner")
    
# CH dfs
ch_vaf = get_mutation_counts_per_gene(chip_main, vaf_colname=vaf_colname)[1].merge(gene_pos_df[["Gene", "CH position"]], how = "left").rename(columns = {"CH position": "Order on oncoprint"}).set_index("Order on oncoprint").drop("Gene", axis = 1)
ch_mut_counts  = get_mutation_counts_per_gene(chip_main, make_a_seperate_category_for_multiple=True, vaf_colname=vaf_colname)[0].merge(gene_pos_df[["Gene", "CH position"]], on="Gene", how = "inner").rename(columns = {"CH position": "Order on oncoprint"}).set_index("Order on oncoprint").drop("Gene", axis = 1)

ch_pt_counts = calculate_patient_counts(chip_main)
ch_pt_counts.columns=["Patient_id", "Count"]

# Generate samples enumerated files
samples_enumerated_dict=generated_samples_enumerated_dict(ch_plotting_df, ch_pt_counts)
samples_enumerated=pd.concat([samples_enumerated_dict["LuPSMA"], samples_enumerated_dict["Cabazitaxel"]], ignore_index=True)

ch_plotting_df = ch_plotting_df.merge(samples_enumerated, how = "left").merge(gene_pos_df, how="left")
ch_plotting_df["Mutation_color"]=ch_plotting_df["Consequence"].map(mut_dict)
ch_pt_counts = ch_pt_counts.merge(samples_enumerated, how = "left")

# PLOTTING   
n_lupsma=ch_plotting_df[ch_plotting_df["Arm"]=="LuPSMA"]["Patient_id"].unique().shape[0]
n_caba=ch_plotting_df[ch_plotting_df["Arm"]=="Cabazitaxel"]["Patient_id"].unique().shape[0]

fig_width=5
fig_height=2.5
tc_height = 2
mutcounts_height = 2
main_height = len(genes_to_include)*2

bar_height = 0.25
bar_width = 0.6

fig = plt.figure(figsize = (fig_width, fig_height))
gs = mpl.gridspec.GridSpec(ncols = 3, nrows = 4, height_ratios = [tc_height, mutcounts_height, main_height, 2], width_ratios = [n_lupsma, n_caba, n_caba/2], hspace = 0, wspace = 0)
ax_mut_counts=fig.add_subplot(gs[2, 2])
ax_legend = fig.add_subplot(gs[3, 0]) 

for i, arm in enumerate(["LuPSMA", "Cabazitaxel"]):
    ax_vaf=fig.add_subplot(gs[0, i])
    ax_mut_counts_patient=fig.add_subplot(gs[1, i], sharex=ax_vaf)
    ax_main=fig.add_subplot(gs[2, i], sharex=ax_vaf, sharey=ax_mut_counts)
    
    samples_enumerated_arm=samples_enumerated_dict[arm]
    chip_main_arm=chip_main[chip_main["Arm"]==arm]
    ch_plotting_df_arm=ch_plotting_df[ch_plotting_df["Arm"]==arm]
    ch_pt_counts_arm=ch_pt_counts[ch_pt_counts["Arm"]==arm]
    
    # Plot grid
    pts = samples_enumerated_arm["Samples_enumerated"].unique()
    genes = ch_plotting_df["CH position"].unique()
    
    for i in genes: 
        for j in pts: 
            ax_main.bar(x = j, bottom = i, height = bar_height, width=bar_width, capstyle='butt', color = "whitesmoke", edgecolor = "none", linewidth = 0.25, zorder = 10)
    
    ax_vaf = plot_highest_VAF_per_patient(ax_vaf, chip_main_arm, samples_enumerated, bar_width = bar_width, vaf_col=vaf_colname)
    ax_main=plot_oncoprint_compressed_tiny(ax_main, ch_plotting_df_arm, bar_height, bar_width)
    remove_frame(ax_main, ["top", "right", "bottom", "left"])
    ax_vaf.set_ylim((0, 12))
    ax_vaf.set_yticks([0, 5, 10])
    ax_vaf.set_yticklabels(["0", "5", "10"])
    ax_main.set_xticks([])
    ax_main.set_xticklabels([])
    ax_vaf.spines['bottom'].set_visible(True)
    ax_mut_counts_patient.spines['bottom'].set_visible(True)
        
    [ax_mut_counts, ax_secondary] = plot_mut_counts_per_gene(ch_mut_counts, ch_vaf, ax_mut_counts, mut_dict, bar_height=bar_height, logscale=True, vaf_col=vaf_colname, show_median=False) # Plot CH mutations
    ax_mut_counts_patient = plot_mut_counts_per_patient_single_df(ax_mut_counts_patient, ch_pt_counts_arm, bar_width)
    
    if arm=="Cabazitaxel":
        ax_main.tick_params("y", left=False, labelleft=False)
        ax_mut_counts_patient.tick_params("y", labelleft=False)
        ax_vaf.tick_params("y", labelleft=False)
        ax_vaf.set_ylabel("")
        ax_mut_counts_patient.set_ylabel("")
    else:
        ax_main.tick_params("y", left=False)
        ax_main.set_yticks(gene_pos_df["CH position"]+bar_height/2)
        ax_main.set_yticklabels(gene_pos_df["Gene"])

ax_main.set_ylim(gene_pos_df["CH position"].min(), gene_pos_df["CH position"].max()+bar_height*1.2)
ax_legend = add_legend_simple(mut_dict, ax_legend)
ax_legend.axis("off")

fig.savefig("/path/to/figures/oncoprint.png")
fig.savefig("/path/to/figures/oncoprint.pdf", transparent=True)