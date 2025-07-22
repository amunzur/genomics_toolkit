#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 6 07:41:52 2024
@author: jbacon
"""

import pandas as pd
import matplotlib.pyplot as plt
import cnvlib
import math
import numpy as np
import os
import seaborn as sns
import re
from collections import defaultdict
from cnvlib import core
from cnvlib import segment


########################
#### User Variables ####
########################

cnrdir = "/groups/wyattgrp/users/amunzur/hla_project/copynum/crpc_panel/fix/snps" #output of cnvkit.fix
vardir_crpc = "/groups/wyattgrp/users/amunzur/hla_project/copynum/crpc_panel/wbc_snps" #Raw VCF output of VarDict
vardir_wes_and_wgs = "/groups/wyattgrp/users/amunzur/hla_project/copynum/wbc_snps_in_snp_grid"
plot_output_dir = "/groups/wyattgrp/users/amunzur/hla_project/copynum/crpc_panel/figures"

##########################
#### Define Functions ####
##########################

# Define the function to apply the logic for copy status to SNP grid
def get_copy_status_backbone(row):
    if row['chromosome'] in ['chrX', 'chrY']:
        return 'N/A'
    if (row['log2'] < -1 and row['baf'] >= 0.6):
        return 'Deep Deletion'
    elif (row['log2'] > 0.7 and row['baf'] >= 0.6):
        return 'Amplification'
    elif (-1.0 <= row['log2'] <= -0.15 and row['baf'] >= 0.6):
        return 'Shallow Deletion'
    elif (0.15 <= row['log2'] <= 0.7 and row['baf'] >= 0.6):
        return 'Gain'
    elif (-0.15 <= row['log2'] <= 0.15 and row['baf'] >= 0.6):
        return 'CN LOH'
    else:
        return 'Neutral'

# Define the function to apply the logic for copy status to violins
def get_copy_status_violin(row):
    if row['Chromosome'] in ['chrX', 'chrY']:
        return 'N/A'
    if (row['log2'] < -1):
        return 'Deep Deletion'
    elif (row['log2'] > 0.7):
        return 'Amplification'
    elif (-1 <= row['log2'] <= -0.3) or (-0.3 <= row['log2'] <= -0.15 and row['baf'] >= 0.6 and row['snp_count'] > 1):
        return 'Shallow Deletion'
    elif (0.3 <= row['log2'] <= 0.7) or (0.15 <= row['log2'] <= 0.3 and row['baf'] >= 0.6 and row['snp_count'] > 1):
        return 'Gain'
    elif (-0.15 <= row['log2'] <= 0.15 and row['baf'] >= 0.6 and row['snp_count'] > 1):
        return 'CN LOH'
    else:
        return 'Neutral'


#####################
#### Import Data ####
#####################

#IMPORT A LIST OF SNPS TO EXCLUDE DUE TO RECURRENT DEVIATION FROM 50% IN WBC DATA
with open("/groups/wyattgrp/users/jbacon/reference/exclude_snps_CN.txt", 'r') as f:
    lines = f.readlines()
    values_to_remove = [(line.split()[0], int(line.split()[1])) for line in lines]

values_set = set(values_to_remove)

#IMPORT ALL FILES FOR PLOTTING. Check SNP only if violin directories are commented out.
if 'cnrdir' in locals() and 'vardir' in locals():
    cnr_files = [os.path.join(cnrdir, f) for f in os.listdir(cnrdir) if f.endswith(".fix")]
    var_files=[os.path.join(vardir, f) for f in os.listdir(cnrdir) if f.endswith(".vcf")]
    common_files = sorted(list(set(cnr_files) & set(var_files)))


# --- Helper ---
def get_patient_id(name):
    return re.split(r'[_-](?:cfDNA|WBC)', os.path.basename(name))[0]

# --- Collect VCF files ---
all_vcf_crpc = [f for f in os.listdir(vardir_crpc) if f.endswith(".vcf")]
all_vcf_weswgs = [f for f in os.listdir(vardir_wes_and_wgs) if f.endswith(".vcf")]

# --- Select best VCF per patient ---
best_var_files = {}
for f in all_vcf_crpc + all_vcf_weswgs:
    pid = get_patient_id(f)
    source = "weswgs" if f in all_vcf_weswgs else "crpc"
    
    if pid not in best_var_files:
        best_var_files[pid] = {"source": source, "files": [f]}
    else:
        if best_var_files[pid]["source"] == "crpc" and source == "weswgs":
            best_var_files[pid] = {"source": "weswgs", "files": [f]}
        elif best_var_files[pid]["source"] == "weswgs" and source == "weswgs":
            best_var_files[pid]["files"].append(f)

# --- Apply WGS preference if needed ---
final_var_files = {}
for pid, entry in best_var_files.items():
    candidates = entry["files"]
    if entry["source"] == "weswgs":
        wgs_files = [f for f in candidates if "WGS" in f]
        chosen = wgs_files[0] if wgs_files else candidates[0]
        path = os.path.join(vardir_wes_and_wgs, chosen)
    else:
        path = os.path.join(vardir_crpc, candidates[0])
    final_var_files[pid] = path

# --- Build dictionary ---
temp_dict = defaultdict(lambda: {'cnv': [], 'var': []})
for f in cnr_files:
    f_full_path=os.path.join(cnrdir, f)
    temp_dict[get_patient_id(f)]['cnv'].append(f_full_path)

for pid, var_path in final_var_files.items():
    temp_dict[pid]['var'].append(var_path)

# --- Expand for multiple CNVs ---
patient_dict = {}
for pid, d in temp_dict.items():
    cnvs = d['cnv']
    vars_ = d['var']
    if len(cnvs) > 1:
        for i, cnv_file in enumerate(cnvs, start=1):
            patient_dict[f"{pid}_{i}"] = {'cnv': [cnv_file], 'var': vars_}
    else:
        patient_dict[pid] = {'cnv': cnvs, 'var': vars_}



pd.DataFrame(patient_dict).T

##########################
#### Begin Processing ####
##########################

for key in patient_dict.keys():
            
    #####################################################
    ######### Load and format the SNP Grid data #########
    #####################################################
    
    if patient_dict[key]["cnv"] and patient_dict[key]["var"]:
        path_cnr = patient_dict[key]["cnv"][0]
        path_var = patient_dict[key]["var"][0]
        
        # Do something with path_cnr and path_var
        print(f"{key}:\n  CNR: {path_cnr}\n  VAR: {path_var}")
    
        cnr=cnvlib.read(path_cnr).drop_low_coverage()
        try:
            var = cnvlib.cmdutil.load_het_snps(path_var, sample_id=None, normal_id=None, min_variant_depth=20)
        # Continue processing var...
        except ValueError as e:
            print(f"Skipping {path_var}: {e}")
            continue
        
        #######VCF FILE FILTERING######
        
        #REMOVE HOMOZYGOTES
        variant_df=var.data
        variant_df=variant_df.loc[variant_df['alt_freq'] < 0.9]
        
        # Remove excluded SNPs
        snp_mask = variant_df.apply(lambda row: (row['chromosome'], row['start']) not in values_set, axis=1)
        variant_df = variant_df[snp_mask]
        
        #REMOVE INDELS
        remove_indel_df = variant_df[(variant_df['ref'].str.len() <= 1) & (variant_df['alt'].str.len() <= 1)]
        
        #REFORMAT ALLELE FREQUENCY PLOT LIKE MATTI DOES
        remove_indel_df['alt_freq']=abs(0.5 - remove_indel_df['alt_freq'])+0.5
        
        #DO NOT ALLOW VARIANTS TO BE EXACTLY 0.5 ELSE CNVLIB PLOTS TWO LINES
        remove_indel_df['alt_freq'] = remove_indel_df['alt_freq'].apply(lambda x: x + 0.001 if x == 0.5 else x)
        
        #REMOVE SNPs WHERE VAF IS MORE THAN 10% OFF FROM AVERAGE OF ADJACENT SNPs
        var_rolling_avg = remove_indel_df['alt_freq'].rolling(window=3, center=True).mean()
        var_outlier_mask = abs(remove_indel_df['alt_freq'] - var_rolling_avg) > 0.10
        
        #APPLY BACK TO VAR OBJECT
        remove_indel_df = remove_indel_df[~var_outlier_mask]
        var.data = remove_indel_df
        
        #######RUN CBS SEGMENTATION######
        
        path_cns_segmentation=path_cnr+".segment"
        cns_df=cnvlib.read(path_cns_segmentation).data
        
        #####
        ###Fix assignment of BAF to each segment here
        #####
        
        # Iterate over each row in cns_df
        for index, cns_row in cns_df.iterrows():
            # Extract the chromosome and start/end positions
            chromosome = cns_row[0]
            start = cns_row[1]
            end = cns_row[2]
            
            # Filter remove_indel_df for rows that overlap with the current segment
            overlapping_rows = remove_indel_df[
                (remove_indel_df['chromosome'] == chromosome) &
                (remove_indel_df['start'] <= end) &
                (remove_indel_df['end'] >= start)
            ]
            
            # Calculate the median alt_freq if there are overlapping rows
            if not overlapping_rows.empty:
                median_alt_freq = overlapping_rows['alt_freq'].median()
            else:
                median_alt_freq = None  # or use np.nan for a numeric representation
            
            # Assign the median_alt_freq to the corresponding row in cns_df
            cns_df.at[index, 'baf'] = median_alt_freq
            
            # cns.data=cns_df
            
            #Assign copy status
            cns_df['copy_status'] = cns_df.apply(get_copy_status_backbone, axis=1)
            
            cns_df.to_csv("/groups/wyattgrp/users/amunzur/COMPOST_BIN/temp.cns", sep="\t", index=False)
            cns = cnvlib.read("/groups/wyattgrp/users/amunzur/COMPOST_BIN/temp.cns")
            
            #Export segmentation to  CSV
            # cns_df=cns_df.drop(labels=['gene'], axis=1)
            # cns_df.insert(0, 'SampleID', ([sample] * len(cns_df)))
            # cns_df.to_csv(segment_output_dir+sample+".tsv", sep='\t', index=False)
            
            #######COPY NUMBER RATIO FILTERING######
            
            #REMOVE PROBES WHERE LOG2 IS MORE THAN 0.3 OFF FROM AVERAGE OF ADJACENT PROBES
            filter_cnr_df=cnr.data
            cnr_rolling_avg = filter_cnr_df['log2'].rolling(window=3, center=True).mean()
            cnr_outlier_mask = abs(filter_cnr_df['log2'] - cnr_rolling_avg) > 0.3
            
            #APPLY BACK TO CNR OBJECT
            filter_cnr_df = filter_cnr_df[~cnr_outlier_mask]
            cnr.data=filter_cnr_df
        
        #Plot SNP Grid scatter first using the cnvlib library from CNVKit
        fig = cnvlib.do_scatter(cnr, segments=cns, segment_color='black', variants=var)
        # fig = cnvlib.do_scatter(cnr, segment_color='black', variants=var)
        fig.set_size_inches(24, 8)
        axes = fig.get_axes() #Extract ax1 (log ratio) and ax2 (BAF) produced by cnvlib.do_scatter
        ax1 = axes[0]
        ax2 = axes[1]
        
        #Color and Coordinate the SNP GRid    
        #set alpha of points    
        for scatter in ax1.collections:
            scatter.set_alpha(0.4)
            
        for scatter in ax2.collections:
            scatter.set_alpha(0.4)
        
        #Color the chromosomes
        # colors = ['#435098', '#86b289', '#8b3730', '#853ea9', '#c4a75c', '#435098', '#86b289', '#8b3730', '#853ea9',
        #           '#c4a75c', '#435098', '#86b289', '#8b3730', '#853ea9', '#c4a75c', '#435098', '#86b289', '#8b3730',
        #           '#853ea9', '#c4a75c', '#435098', '#86b289', '#8b3730', '#853ea9']
        # for i, collection in enumerate(ax1.collections):
        #    collection.set_color(colors[i])
           
        # for i, collection in enumerate(ax2.collections):
        #    collection.set_color(colors[i])
          
        #No attribute for SNP segment line color. Make them black and bring forward here.
        for line in ax2.get_lines():
            x_data = line.get_xdata()
            y_data = line.get_ydata()
            
            line.set_color('black')
            line.set_zorder(10) 
        
        #Set Plot Ylims
        log2_segvals=cns["log2"]
        
        seg_max=log2_segvals.max()
        rounded_max_value = math.ceil(seg_max * 2) / 2
        
        if rounded_max_value > 1.5:
            y_lim_high=rounded_max_value
        else:
            y_lim_high=1.5
        
        ax1.set_ylim(-1.5, y_lim_high)
        
        # Add gridlines
        for y in np.arange(-2, y_lim_high, 0.5):
            ax1.axhline(y = y, color = 'k', linestyle = 'dotted', linewidth=0.75, alpha=0.3)
            
        for y in np.arange(0.5, 1.0, 0.1):
            ax2.axhline(y = y, color = 'k', linestyle = 'dotted', linewidth=0.75, alpha=0.3)
        
        #customize plot   
        ax1.set_ylabel("Log Ratio", fontsize=14)
        ax1.tick_params(axis='both', labelsize=12)
        
        basename = os.path.basename(path_cnr)
        clean_name = re.sub(r'.cnn.*$', '', basename)
        title=f"{clean_name}\nVariants file: {path_var}\nCNV file: {path_cnr}"
        ax1.set_title(title, pad=16, fontsize=18)
        
        ax2.set_ylabel("BAF", fontsize=14, labelpad=16)
        ax2.tick_params(axis='both', labelsize=12)
        ax2.set_ylim(0.5, 1.0)
        
        # Line through HLA regions
        hla_location=29941260           
        
        fig.savefig(os.path.join(plot_output_dir, clean_name+".png"), dpi=300, bbox_inches='tight')
        
