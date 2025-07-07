

import pandas as pd
import os
import math

def return_gene_start_end(gene, panel_df):
    """
    Given a gene name returns start and end from the panel target locations.
    """
    panel_df_subset=panel_df[panel_df["gene"].str.contains(gene)]
    start=panel_df_subset["start"].min()
    end=panel_df_subset["start"].max()
    
    return start, end

import pandas as pd

def return_exon_coordinates(gene, path_ref_gene):
    """
    Extracts exon coordinates for the canonical transcript of a given gene.
    """
    columns = [
        "bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
        "exonCount", "exonStarts", "exonEnds", "score", "name2",
        "cdsStartStat", "cdsEndStat", "exonFrames"
    ]
    df = pd.read_csv(path_ref_gene, sep="\t", names=columns)
    
    # Filter for the gene
    gene_df = df[df["name2"] == gene]
    
    if gene_df.empty:
        print(f"No data found for gene: {gene}")
        return None
    
    # Pick the longest transcript (you can modify this logic)
    gene_df["tx_length"] = gene_df["txEnd"] - gene_df["txStart"]
    canonical_transcript = gene_df.loc[gene_df["tx_length"].idxmax()]
    
    # Extract exon coordinates
    exon_starts = list(map(int, canonical_transcript["exonStarts"].strip(",").split(",")))
    exon_ends = list(map(int, canonical_transcript["exonEnds"].strip(",").split(",")))
    
    # Store exon info in DataFrame
    exons = [{"chrom": canonical_transcript["chrom"], "start": start, "end": end, "strand": canonical_transcript["strand"]}
             for start, end in zip(exon_starts, exon_ends)]
    
    exon_df = pd.DataFrame(exons)
    
    return exon_df

def annotate_exons(exons_df, ax_cn):
    """
    Annotates exons as shaded regions on the given Matplotlib axis.
    """
    for i, row in exons_df.iterrows():
        ax_cn.axvspan(row["start"], row["end"], color="deepskyblue", alpha=0.2)
    
    return(ax_cn)

def return_ploidy_purity(sample_name, path_ct_frac):
    """
    
    """
    df=pd.read_csv(path_ct_frac, sep="\t")
    df=df[df["sample"]==sample_name]
    
    dl=df[df["sample"]==sample_name]["diploid level"].values[0]
    ctdna=df[df["sample"]==sample_name]["ctdna"].values[0]
    ploidy=df[df["sample"]==sample_name]["ploidy"].values[0]
    
    result_dict={"dl": dl, "ctdna": ctdna, "ploidy": ploidy}
    return(result_dict)

def get_logr_line_from_cn(cn, ct_frac, diploid_level):
    logr=math.log2((cn*ct_frac)/2+1-ct_frac)+diploid_level
    return(logr)


def plot_cn_logratio_lines(ct_frac, diploid_level, ax, xpos_for_text=3):
    """
    Plots those horizontal lines in the plot for ploidy.
    """
    # STEP 1.Plot CN logratio lines
    color_dict={
        8: "darkred",
        7: "darkred",
        6: "darkred",
        5: "darkred",
        4: "red",
        3: "lightcoral",
        2: "gray",
        1: "cornflowerblue"
    }
    for cn in np.arange(1, 6):
        cn_logr=get_logr_line_from_cn(cn, ct_frac, diploid_level)
        logr_line_color=color_dict[cn]
        ax.axhline(y=cn_logr, color=logr_line_color, linewidth=1, linestyle="dashed")
        ax.text(xpos_for_text, cn_logr, cn, color=logr_line_color, fontsize=8, rotation=0, ha='center')
    
    return(ax)

def main(path_segmentation, path_fix, path_hla_panel, path_ct_frac, path_ref_gene, path_output, sample_name, gene):
    """
    Main plotting function.
    """
    panel_df=pd.read_csv(path_hla_panel, sep="\t", header=None, names=['chrom', 'start', 'end', 'gene'])
    segmentation_main=pd.read_csv(path_segmentation, sep="\t").replace(r"\bSETD1B\b", "SETDB1", regex=True) 
    
    logr_main=pd.read_csv(path_fix, sep="\t")
    logr_main["Midpoint"]=(logr_main["start"]+logr_main["end"])/2
    
    fig = plt.figure(figsize=(8, 6))
    outer_gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace = 0.2, wspace = 0) # outer most gs with 3 rows  
    
    ax_cn=plt.subplot(outer_gs[0])
    ax_snps=plt.subplot(outer_gs[1])    
    
    gene_start, gene_end=return_gene_start_end(gene, panel_df)
    gene_segmentation=segmentation_main[segmentation_main["gene"].str.contains(gene)]
    gene_logr=logr_main[logr_main["gene"].str.contains(gene)]   
    
    # Return ctdna etc
    ploidy_dict=return_ploidy_purity(sample_name, path_ct_frac)
    ploidy=round(ploidy_dict["ploidy"], 2)
    ctdna=round(ploidy_dict["ctdna"]*100, 2)
    dl=ploidy_dict["dl"]
    ax_title=f"{gene}\n{sample_name}\nctDNA%: {ctdna} Diploid level: {dl} Ploidy: {ploidy}"
    
    # Plot segmentation lines
    for i, row in gene_segmentation.iterrows():
        segment_logr=row["log2"]
        segment_start=row["start"]
        segment_end=row["end"]
        ax_cn.hlines(y=segment_logr, xmin=segment_start, xmax=segment_end, color='black')
    
    # Plot logr scatter
    for i, row in gene_logr.iterrows():
        xpos=row["Midpoint"]
        logr=row["log2"]
        ax_cn.scatter(xpos, logr, s=7, color="black", zorder=100)   
    
    # Annotate exons
    exons_df=return_exon_coordinates(gene, path_ref_gene)
    ax_cn=annotate_exons(exons_df, ax_cn)   
    
    ax_cn.set_xlim((gene_start-1000, gene_end+1000))
    ax_cn.set_ylim((-2, 2))
    ax_cn.set_yticks([-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2])
    ax_cn.set_yticklabels(["-2", "-1.5", "-1", "-0.5", "0", "0.5", "1", "1.5", "2"])
    ax_cn.set_ylabel("Coverage logratio")
    ax_cn.set_title(ax_title)
    ax_cn.spines[["top", "right"]].set_visible(False)   
    
    ax_snps.set_ylabel("Het SNP allele frequency")
    ax_snps.spines[["top", "right"]].set_visible(False) 
    
    # Plot CN lines
    ax_cn=plot_cn_logratio_lines(ctdna/100, dl, ax_cn, xpos_for_text=gene_start)
    
    outer_gs.tight_layout(fig)
    fig.savefig(path_output)


path_segmentation="/groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/fix/genes/GU-21-633_FiT_22Oct2021_1.cnn.fix.segment"
path_fix="/groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/fix/genes/GU-21-633_FiT_22Oct2021_1.cnn.fix"
path_hla_panel="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/hla_panel/target_regions.bed"
path_ref_gene="/groups/wyattgrp/users/amunzur/pipeline/resources/references/refGene.txt"
path_ct_frac="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_ploidy.tsv"
path_sample_list="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv"
path_gene_list="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/hla_panel/gene_list.tsv"
dir_figures="/groups/wyattgrp/users/amunzur/hla_pipeline/results/copynum/zoomed_in_plots"

genes = [line.strip() for line in open(path_gene_list)]

sample_df=pd.read_csv(path_sample_list, sep="\t")
tumor_samples=sample_df["Tumor_name"].unique()

for tumor_sample in tumor_samples:
    for gene in genes:
        path_output=os.path.join(dir_figures, f"{tumor_sample}_{gene}")
        main(path_segmentation, path_fix, path_hla_panel, path_ct_frac, path_ref_gene, path_output, tumor_sample, gene)
        print(f"{tumor_sample} / {gene}")