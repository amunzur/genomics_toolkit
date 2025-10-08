
# Script to make IGV snapshots with tumor and WBC bams shown together. 

# For tumor + WBC BAMs:
# python make_igv_snapshots.py \
#   --path_variants /path/to/variants.csv \
#   --PATH_batch /path/to/output_batch.txt \
#   --DIR_snapshots /path/to/snapshots \
#   --prefix "" \
#   --suffix ".bam" \
#   --given_range 100

# For unpaired BAMs:
# python make_igv_snapshots.py \
#   --path_variants /path/to/wbc_variants.csv \
#   --PATH_batch /path/to/output_batch_wbc.txt \
#   --DIR_snapshots /path/to/snapshots_wbc \
#   --unpaired \
#   --given_range 150 \

"""
Created on Thu May 27 13:45:38 2021
@author: amunzur
"""
import os
import pandas as pd
import argparse

def apply_prefix_suffix(path, prefix, suffix):
	"""
	Apply optional prefix and/or suffix to a given file path.
	"""
	if prefix:
		path = str(prefix) + str(path)
	if suffix:
		path = str(path) + str(suffix)
	return path

def make_igv_batch_script(df, PATH_batch, DIR_snapshots, prefix, suffix, given_range, unpaired):
	"""
	Generate an IGV batch script to create snapshots of tumor and/or WBC BAM files.
	"""
	for index, row in df.iterrows(): # iterate through each snv in the variants file
		
		if "Path_bam" in df.columns and unpaired:
			BAM_wbc = apply_prefix_suffix(row["Path_bam"], prefix, suffix)
		elif "Path_bam" not in df.columns and unpaired:
			BAM_wbc = apply_prefix_suffix(row["Sample_name"], prefix, suffix) # Can be tumor or WBC
		elif "Path_bam" not in df.columns and not unpaired and "Sample_name_n" in df.columns and "Sample_name_n" in df.columns:
			BAM_wbc = apply_prefix_suffix(row["Sample_name_n"], prefix, suffix)
			BAM_tumor = apply_prefix_suffix(row["Sample_name_t"], prefix, suffix)
		else:
			BAM_wbc = apply_prefix_suffix(row["Path_bam_n"], prefix, suffix)
			BAM_tumor = apply_prefix_suffix(row["Path_bam_t"], prefix, suffix)
				
		position = int(row["Position"])
		
		if given_range: # if the user wants to see a given range, consider the end position of the range as well
			start_position = str(int(position - given_range/2))
			end_position = str(int(position + given_range/2))	
		
		# os.remove(IGV_script)
		os.makedirs(DIR_snapshots, exist_ok = True) # make the snapshost dir if it doesnt exist already
		if unpaired:
			output_file_name = row["Gene"] + "_" + str(row["Protein_annotation"]) + "_" + row["Chrom"] + "_" + str(row["Position"]) + "_" + row["Sample_name"] + ".png" # one snapshot for each variant		
		else:
			output_file_name = row["Gene"] + "_" + str(row["Protein_annotation"]) + "_" + row["Chrom"] + "_" + str(row["Position"]) + "_" + row["Sample_name_t"] + ".png" # one snapshot for each variant		
		# Begin compiling the batch file
		with open(PATH_batch, 'a') as the_file:
			the_file.write('new\n')
			if index == 0: the_file.write('genome hg38\n') # only load the genome if it is the first time we are starting IGV
			if unpaired:
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
	with open(PATH_batch, 'a') as the_file: the_file.write('exit') # append an exit statement at the end so that IGV closes on its own
	print("/home/amunzur/IGV_Linux_2.11.3/igv.sh --batch ", PATH_batch) # print the exact command needed to turn IGV to terminal

def main():
	parser = argparse.ArgumentParser(description="Generate IGV batch script for tumor + WBC BAM visualization")
	parser.add_argument("--path_variants", required=True, help="Path to variants CSV file")
	parser.add_argument("--PATH_batch", required=True, help="Output path for IGV batch script")
	parser.add_argument("--DIR_snapshots", required=True, help="Directory where snapshots will be saved")
	parser.add_argument("--prefix", help="Prefix to add to BAM paths (optional)")
	parser.add_argument("--suffix", help="Suffix to add to BAM paths (optional)")
	parser.add_argument("--given_range",type=int,default=200,help="Genomic range around the position to display")
	parser.add_argument("--unpaired",action="store_true",help="If set, only load WBC BAMs (no tumor BAMs)")
	args = parser.parse_args()
		
	# Remove old batch file if it exists
	try:
		os.remove(args.PATH_batch)
	except FileNotFoundError:
		pass
	
	# Load dataframe
	df = pd.read_csv(args.path_variants)

    # Ensure required columns exist
	# required_cols = (
	# 	["Patient_id", "Protein_annotation", "Gene", "Chrom", "Position", "Sample_name"]
	# 	if args.unpaired
	# 	else ["Patient_id", "Protein_annotation", "Gene", "Chrom", "Position", "Sample_name_t", "Path_bam_t", "Path_bam_n",
	# 	]
	# )
	# df = df[required_cols]

	make_igv_batch_script(
		df,
		args.PATH_batch,
		args.DIR_snapshots,
		args.prefix,
		args.suffix,
		args.given_range,
		args.unpaired,
	)

if __name__ == "__main__":
	main()
