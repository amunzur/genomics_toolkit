#!/bin/python

"""
Generates per-patient plots of CH mutations.

For each patient:
- Plots mutation VAFs by gene
- Shades gene backgrounds and annotates protein changes
- Outputs PNG and PDF figures

Good to run this for every cohort I work on.

Inputs:
- --path_muts: Path to curated mutation list
- --vafcolname: Column name for VAF values
- --path_sample_information: Sample metadata with timepoints
- --dir_figures: Output directory for the figures

Example:
python per_patient_plots_JUST_CHIP.py \
    --path_muts CHIP_SSCS2_curated.csv \
    --vafcolname VAF_n \
    --path_sample_information sample_information.tsv \
    --dir_figures ./patient_profiles
"""

import pandas as pd
import numpy as np
import os
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import matplotlib.patches as patches
import argparse

"""
Visualizes CH mutations for each patient.
"""

def plot_gene_shading(mutations_df, ax, color_dict, rect_color="whitesmoke"):
    """
    """
    ymin = 0  # Set ymin to 0
    ymax=ax.get_ylim()[1]
    rect_midx_point_min=0.5
    rect_midx_point = rect_midx_point_min
    rect_width = 1  
    for i, gene in enumerate(color_dict.keys()):
        rect_width = (len(color_dict[gene]) / 3)*0.8  # Define the width of the rectangle
        rect_midx_point = rect_midx_point_min + i  # Calculate the midpoint on x
        rect_x = rect_midx_point - rect_width / 2  # Calculate the bottom-left x-coordinate
        # Create the rectangle
        rect = patches.Rectangle(
            (rect_x, ymin),
            rect_width,
            ymax - ymin,
            edgecolor='None',
            facecolor=rect_color,
            zorder=5
        )
        # Add the rectangle to the axis
        ax.add_patch(rect)
        ax.text(rect_midx_point, ymax*0.95, gene, ha='center', va='bottom', fontsize=7, zorder=10, fontstyle="italic")
    
    ax.set_xlim(0, rect_midx_point+rect_width/2)
    
    return(ax)


def plot_patient_mutations(mutations_df, name_vaf_column, ax, color_dict, xval, vafcolname="VAF%", add_legend=True):
    """
    Plots mutation vafs.
    """
    # ax.set_xlim(0.5, 1.5)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_ylabel(vafcolname)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for i, gene in enumerate(color_dict.keys()):
        gene_df=mutations_df[mutations_df["Gene"]==gene]
        for j, row in gene_df.iterrows():
            color = color_dict.get(gene, "black")  # Use black as the default color if the gene is not in color_dict
            jitter=np.random.uniform(-0.2, 0.2, 1)
            ax.scatter(i+0.5, row[name_vaf_column], s=20, edgecolor=None, color="black", zorder=999)
    
    ax.set_ylim(0, ax.get_ylim()[1]*1.1)
    
    # Legend
    if add_legend:
        legend_colors = colors.values()
        legend_labels = colors.keys()
        legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
        ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize = 8, handletextpad=0.05)
    
    return(ax)

def annotate_patient_muts(patient_muts, ax_annotations, gene_order, vafcolname):
    """
    Adds protein alterations as annotations besides gene names.
    """
    ax_annotations.annotate("Protein alterations", (0.5, 0.95), ha='center', va='center', fontsize = 9, xycoords='axes fraction')
    
    ypos_max=0.92
    for i, gene in enumerate(gene_order):
        df_subset=patient_muts[patient_muts["Gene"]==gene]
        for j, row in df_subset.iterrows():
            protein_alteration=row["Protein_annotation"]
            if protein_alteration is np.nan:
                protein_alteration="Splice site"
            vaf=round(row[vafcolname], 2)
            txt=f"{gene} {protein_alteration} {vaf}%"
            ypos=ypos_max-0.07
            ax_annotations.annotate(txt, (0.1, ypos), ha='left', va='center', fontsize = 9, xycoords='axes fraction')
            ypos_max=ypos
    
    # AES
    # ax.set_xlim(0.5, 1.5)
    ax_annotations.set_xticks([])
    ax_annotations.set_yticks([])
    ax_annotations.set_xticklabels([])
    ax_annotations.set_yticklabels([])
    ax_annotations.spines[["top", "right", "left", "bottom"]].set_visible(False)
    return(ax_annotations)


def plot_lines_between_two_scatters(list1, list2, genes, color_dict, ax):
    """
    Connects two given dots with a line.
    """
    for i in range(len(list1)):
        gene = genes[i]
        color = color_dict.get(gene, "black")  # Default to black if gene not in color_dict
        ax.plot([1, 2], [list1[i], list2[i]], color=color, linewidth=2)  # Connect points for each gene
    
    return(ax)

def get_time_diff(patient, path_sample_info):
    """
    Given a patient ID, returns the time difference between baseline and progression in months.
    """
    sample_info = pd.read_csv(path_sample_info, sep = "\t")
    sample_info["Date_collected"] = pd.to_datetime(sample_info["Date_collected"], format="%Y%b%d")
    
    baseline_time = sample_info[(sample_info["Patient_id"] == patient) & (sample_info["Timepoint"] == "Baseline")]["Date_collected"].iloc[0]
    prog_time = sample_info[(sample_info["Patient_id"] == patient) & (sample_info["Timepoint"] == "FirstProgression")]["Date_collected"].iloc[0]
    
    time_diff_in_mo = round((prog_time - baseline_time).days/30.44)
    return(time_diff_in_mo)

def generate_color_palette_for_genes(genes_list):
    """
    Generates a color palette, each gene gets a distinct color.
    """
    palette = sns.color_palette("Set2")
    extended_palette = palette * (len(genes_list) // len(palette) + 1) # Extend the palette if needed
    color_dict = {gene: extended_palette[i] for i, gene in enumerate(genes_list)} # Map each gene to a color
    return(color_dict)

def align_ctdna_dataframes(baseline_df, progression_df):
    """
    Aligns two DataFrames containing baseline and progression ctDNA data by adding missing entries.
    Missing entries will have VAF set to 0, Effect set to None, and Genomic position set to None.
    
    Parameters:
    - baseline_df (pd.DataFrame): DataFrame containing baseline ctDNA data.
    - progression_df (pd.DataFrame): DataFrame containing progression ctDNA data.
    
    Returns:
    - (pd.DataFrame, pd.DataFrame): Aligned baseline and progression DataFrames.
    """
    # Create unique (Gene, Genomic position) pairs from both DataFrames
    baseline_pairs = set(zip(baseline_df["Gene"], baseline_df["Genomic position"]))
    progression_pairs = set(zip(progression_df["Gene"], progression_df["Genomic position"]))
    
    # Union of all (Gene, Genomic position) pairs
    all_pairs = baseline_pairs.union(progression_pairs)
    
    # Add missing entries to baseline
    missing_in_baseline = all_pairs - baseline_pairs
    if missing_in_baseline:
        missing_baseline_rows = pd.DataFrame(
            {
                "Patient": [progression_df["Patient"].iloc[0]] * len(missing_in_baseline),
                "Timepoint": ["Baseline"] * len(missing_in_baseline),
                "Gene": [pair[0] for pair in missing_in_baseline],
                "Genomic position": [pair[1] for pair in missing_in_baseline],
                "Variant Allele Frequency (VAF%)": [0] * len(missing_in_baseline),
                "Effect": [None] * len(missing_in_baseline),
            }
        )
        baseline_df = pd.concat([baseline_df, missing_baseline_rows], ignore_index=True)
    
    # Add missing entries to progression
    missing_in_progression = all_pairs - progression_pairs
    if missing_in_progression:
        missing_progression_rows = pd.DataFrame(
            {
                "Patient": [baseline_df["Patient"].iloc[0]] * len(missing_in_progression),
                "Timepoint": ["Progression"] * len(missing_in_progression),
                "Gene": [pair[0] for pair in missing_in_progression],
                "Genomic position": [pair[1] for pair in missing_in_progression],
                "Variant Allele Frequency (VAF%)": [0] * len(missing_in_progression),
                "Effect": [None] * len(missing_in_progression),
            }
        )
        progression_df = pd.concat([progression_df, missing_progression_rows], ignore_index=True)
    
    # Sort DataFrames (optional, for readability)
    baseline_df = baseline_df.sort_values(by=["Gene", "Genomic position"]).reset_index(drop=True)
    progression_df = progression_df.sort_values(by=["Gene", "Genomic position"]).reset_index(drop=True)
    
    return baseline_df, progression_df

def generate_legend(ax, color_dict):
    """
    Makes a legend by matching colors with gene names.
    """ 
    legend_colors = color_dict.values()
    legend_labels = color_dict.keys()
    legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="lower left", frameon=False, fontsize = 10, handletextpad=0.1)
    
    ax.axis('off')
    
    return(ax)

def main():
    parser = argparse.ArgumentParser(description="Makes individual plots for all patients.")
    parser.add_argument("--path_muts", required=True, help="Path to curated mutation list.")
    parser.add_argument("--vafcolname", help="The name of the column in the curated mutation list that has the VAF information.")
    parser.add_argument("--path_sample_information", required=True, help="Path to sample information. Can accomodate with or without header.")
    parser.add_argument("--dir_figures", required=True, help="DIR where the figures will be saved.")
    args = parser.parse_args()
    
    path_muts=args.path_muts
    vafcolname=args.vafcolname
    path_sample_information=args.path_sample_information
    dir_figures=args.dir_figures
    
    muts=pd.read_csv(path_muts)
    if muts[vafcolname].min()<0.25:
        muts[vafcolname]=muts[vafcolname]*100
    
    sample_info = pd.read_csv(path_sample_information, sep = "\t")
    all_pts=sample_info["Patient_id"].unique()
    
    for patient in all_pts:
        print(f"Working on {patient}.")

        # STEP1. GENERATE RELEVANT PATIENT DFS AND COLOR PALETTES.
        patient_muts = muts[muts["Patient_id"] == patient]
        genes = patient_muts["Gene"].drop_duplicates().tolist()
        colors = generate_color_palette_for_genes(genes)

        # STEP2. PLOTTING SCATTER POINTS.
        fig = plt.figure(figsize=(6.5, 3.5))
        fig.suptitle(patient)
        gs = gridspec.GridSpec(1, 2, width_ratios = [1, 0.8], hspace=0.4)

        ax_muts=plt.subplot(gs[0])
        ax_annotations=plt.subplot(gs[1])

        ax_muts=plot_patient_mutations(mutations_df = patient_muts, name_vaf_column = vafcolname, ax = ax_muts, color_dict = colors, xval = 1, add_legend=False)
        ax_annotations=annotate_patient_muts(patient_muts=patient_muts, ax_annotations=ax_annotations, gene_order=colors.keys(), vafcolname=vafcolname)
        ax_muts=plot_gene_shading(patient_muts, ax=ax_muts, color_dict=colors, rect_color="whitesmoke")

        gs.tight_layout(fig)
        fig.savefig(os.path.join(dir_figures, patient+".png"))
        fig.savefig(os.path.join(dir_figures, patient+".pdf"))

if __name__ == "__main__":
    main()

"""
python /path/to/toolkit/per_patient_plots_JUST_CHIP.py \
    --path_muts /path/to/results/variant_calling/CHIP_SSCS2_curated.csv \
    --vafcolname VAF_n \
    --path_sample_information /path/to/resources/sample_lists/sample_information.tsv \
    --dir_figures /path/to/results/figures/patient_profiles
"""
