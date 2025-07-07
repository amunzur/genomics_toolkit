# Script to make IGV snapshots with tumor and WBC bams shown together. 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 13:45:38 2021
@author: amunzur
"""
import os
import pandas as pd

def make_igv_batch_script(
	df, 
	PATH_batch, 
	DIR_snapshots, 
	add_prefix, 
	add_suffix, 
	given_range, 
	WBC_only):
	for index, row in df.iterrows(): # iterate through each snv in the maf
		if WBC_only: 
			BAM_wbc = str(add_prefix + str(row["Path_bam_n"])+ add_suffix)		
		else:
			BAM_wbc = str(add_prefix + str(row["Path_bam_n"])+ add_suffix)		
			BAM_tumor = str(add_prefix + str(row["Path_bam_t"]) + add_suffix)
		position = int(row["Position"])
		if given_range: # if the user wants to see a given range, consider the end position of the range as well
			start_position = str(int(position - given_range/2))
			end_position = str(int(position + given_range/2))	
		# os.remove(IGV_script)
		os.makedirs(DIR_snapshots, exist_ok = True) # make the snapshost dir if it doesnt exist already
		if WBC_only:
			output_file_name = row["Gene"] + "_" + str(row["Protein_annotation"]) + "_" + row["Chrom"] + "_" + str(row["Position"]) + "_" + row["Sample_name"] + ".png" # one snapshot for each variant		
		else:
			output_file_name = row["Gene"] + "_" + str(row["Protein_annotation"]) + "_" + row["Chrom"] + "_" + str(row["Position"]) + "_" + row["Sample_name_t"] + ".png" # one snapshot for each variant		
		# Begin compiling the batch file
		with open(PATH_batch, 'a') as the_file:
			the_file.write('new\n')
			if index == 0: the_file.write('genome hg38\n') # only load the genome if it is the first time we are starting IGV
			if WBC_only:
				the_file.write(str("load " + BAM_wbc + '\n'))
			else:
				the_file.write(str("load " + BAM_tumor + '\n'))
				the_file.write(str("load " + BAM_wbc + '\n'))
			the_file.write(str("snapshotDirectory " + DIR_snapshots + '\n'))
			chrom = str(row["Chrom"])
			if given_range: the_file.write(str('goto ' + chrom + ":" + start_position + "-" + end_position + "\n"))
			else: the_file.write(str('goto ' + chrom + ":" + str(position) + "\n"))
			the_file.write('sort start location\n')
			the_file.write('sort base\n')
			the_file.write('maxPanelHeight -1\n')
			the_file.write(str('snapshot ' + output_file_name + '\n'))
			the_file.write("\n")
	with open(PATH_batch, 'a') as the_file: the_file.write('exit') # append an exist statement at the end so that IGV closes on its own
	print("/home/amunzur/IGV_Linux_2.11.3/igv.sh --batch ", PATH_batch) # print the exact command needed to turn IGV to terminal

PATH_variants_df = "/groups/wyattgrp/users/amunzur/hla_pipeline/results/variant_calling/ctdna_mutations.csv" # variants csv file
PATH_batch = "/groups/wyattgrp/users/amunzur/hla_pipeline/workflow/scripts/batch_scripts/IGV_batch_ctdna_mutations.txt" # where the batch script is saved
DIR_snapshots = "/groups/wyattgrp/users/amunzur/hla_pipeline/results/figures/IGV_snapshots" # dir where snapshots will be saved
given_range = 100
add_prefix = "" # if .bam suffix needs to be added 
add_suffix = ""

df = pd.read_csv(PATH_variants_df)[["Patient_id", "Protein_annotation", "Gene", "Chrom", "Position", "Sample_name_t", "Path_bam_t", "Path_bam_n"]] # assuming these cols exist in the variants file
# df["Path_bam_n"]="/groups/wyattgrp/users/amunzur/lu_chip/results/bam/SSCS2_final/"+df["Sample_name"]+".bam"

try:
	os.remove(PATH_batch)
except:
	print("Error while deleting batch files")

make_igv_batch_script(df, PATH_batch, DIR_snapshots, add_prefix, add_suffix, given_range, WBC_only = False)

# For progression cfDNA only mutations
path="/groups/wyattgrp/users/amunzur/ironman_ch/results/variant_calling/CHIP_SSCS2.csv"
df = pd.read_csv(path)
df = progression_only[["Patient_id", "Protein_annotation", "Gene", "Chrom", "Position", "Sample_name", "bam_path"]] # assuming these cols exist in the variants file
df.columns = ["Patient_id", "Protein_annotation", "Gene", "Chrom", "Position", "Sample_name", "Path_bam_n"]


