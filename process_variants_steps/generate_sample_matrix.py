"""
Report read support for all variants in all samples. 
Requires mpileups to be computed first.
"""

import pandas as pd
import os
import subprocess
import time
import re

path_muts="/groups/wyattgrp/users/amunzur/ironman_ch/results/variant_calling/IRONMAN_CH_calls.csv"
dir_filtered_mpileups="/groups/wyattgrp/users/amunzur/ironman_ch/results/metrics/mpileup_for_matrix"
dir_mpileups="/groups/wyattgrp/users/amunzur/ironman_ch/results/metrics/mpileup/SSCS2"
path_sample_information="/groups/wyattgrp/users/amunzur/ironman_ch/resources/sample_lists/sample_information.tsv"
dir_bams="/groups/wyattgrp/users/amunzur/ironman_ch/results/data/bam/SSCS2_final"
path_matrix="/groups/wyattgrp/users/amunzur/ironman_ch/results/variant_calling/ironman_mutation_sample_matrix.csv"

def generate_unique_mutation_list(path_muts):
    """
    Generates a list of unique mutations to run grep on these locations.
    """
    muts=pd.read_csv(path_muts)
    unique_muts=muts[["Gene", "Chrom", "Position", "Alt", "Type", "Protein_annotation"]].drop_duplicates()
    
    return(unique_muts)

def run_commands_in_batches(command_list, batch_size=50, delay_between_batches=1):
    """
    Runs a list of shell commands in batches, waiting for all commands in a batch to finish 
    before starting the next batch.
    
    Parameters:
    - command_list (list of str): List of shell commands to run.
    - batch_size (int): Number of commands to run simultaneously in each batch.
    - delay_between_batches (int): Time to wait (in seconds) between batches.
    
    Returns:
    - None
    """
    # Split commands into batches
    batches = [command_list[i:i + batch_size] for i in range(0, len(command_list), batch_size)]
    
    for batch_index, batch in enumerate(batches, start=1):
        processes = []
        print(f"Running batch {batch_index}/{len(batches)}...")
        # Start all commands in the current batch
        for command in batch:
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            processes.append(process)
        
        # Wait for all processes in the batch to complete
        for process in processes:
            process.wait()  # Wait for the process to finish
        print(f"Batch {batch_index} completed.")
        # Delay between batches
        if delay_between_batches > 0:
            time.sleep(delay_between_batches)
    print("All commands have been executed.")

def check_if_mutation_was_called(sample, chrom, pos, protein_annot, alt, path_muts):
    """
    Given info about a mutation, checks whether it was called in the given sample.
    """
    muts=pd.read_csv(path_muts)
    muts["Position"]=muts["Position"].astype(int)
    pos=int(pos)
    
    muts_filtered=muts[(muts["Patient_id"]==sample) & (muts["Position"]==pos) & (muts["Chrom"]==chrom) & (muts["Protein_annotation"]==protein_annot)]
    if muts_filtered.empty:
        mut_called= False
    else:
        mut_called=True
    
    return mut_called

def parse_insertions(consensus_read_bases, alt):
    # Normalize ALT to lowercase to match against lowercase insertions
    alt = alt.lower()[1:] # To allow flexible match without the initial base
    # Regex pattern to match insertions of the form +<len><sequence>
    insertion_pattern = re.compile(r'\+(\d+)([acgtnACGTN]+)', re.IGNORECASE)
    matches = insertion_pattern.findall(consensus_read_bases)
    count = 0
    for length_str, sequence in matches:
        if sequence.lower() == alt:
            count += 1
    
    return count

def count_deletions(pileup_seq, expected_del=None):
    """
    Count deletions in pileup consensus string.
    Args:
        pileup_seq (str): Pileup consensus string.
        expected_del (str): Optional. If provided, only count this deletion sequence.
    Returns:
        int: Number of deletion-supporting reads.
    """
    i = 0
    count = 0
    while i < len(pileup_seq):
        if pileup_seq[i] == '*':
            count += 1
            i += 1
        elif pileup_seq[i] == '-':
            i += 1
            length_str = ''
            while i < len(pileup_seq) and pileup_seq[i].isdigit():
                length_str += pileup_seq[i]
                i += 1
            if not length_str:
                continue
            length = int(length_str)
            del_seq = pileup_seq[i:i + length]
            i += length
            if expected_del is None or del_seq.upper() == expected_del.upper():
                count += 1
        else:
            i += 1
    
    return(count)

def return_depth_and_n_altered_from_filtered_mpileup(path_filtered_mpileup, chrom, pos, mut_type, alt):
    """
    Return the depth and the number of alt reads from grepped mpileup output.
    """
    df=pd.read_csv(path_filtered_mpileup, sep="\t", header=None, names=["chrom", "position", "ref", "Coverage", "Consensus_read_bases", "MapQ"])
    df = df[(df["chrom"] == chrom) & (df["position"] == int(pos))]
    
    if df.empty:
        return {"n_alt": 0, "depth": 0, "VAF": 0.0}
    
    depth = int(df.iloc[0]["Coverage"]) # Extract depth
    consensus_read_bases = df.iloc[0]["Consensus_read_bases"] # Extract the consensus read bases
    
    # Determine n_alt based on mutation type
    if mut_type.upper() == "SNV":
        n_alt = consensus_read_bases.lower().count(alt.lower())
    elif mut_type.upper() == "INSERTION":
        n_alt = parse_insertions(consensus_read_bases, alt)
        # Match insertions with the format +<length><alt>
        # insertion_pattern = rf"\+{len(alt)}{alt}"
        # n_alt = len(re.findall(insertion_pattern, consensus_read_bases))
    elif mut_type.upper() == "DELETION":
        n_alt=count_deletions(consensus_read_bases, expected_del=alt)
        # Match deletions with the format -<length><ref>
        # deletion_pattern = rf"-{len(alt)}{alt}"
        # n_alt = len(re.findall(deletion_pattern, consensus_read_bases))
    else:
        raise ValueError("Invalid mutation type. Use 'SNV', 'INSERTION', or 'DELETION'.")
    
    # Calculate VAF
    vaf = (n_alt / depth) * 100 if depth > 0 else 0.0
    
    # Return results as a dictionary
    result_dict = {"n_alt": n_alt, "depth": depth, "VAF": vaf}
    return result_dict

sample_list=pd.read_csv(path_sample_information, sep="\t", header=None, names=["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
sample_list["Patient_id"]=sample_list["Patient_id"].str.replace("-", "_")

# STEP 1. RUN GREP ON ALL SAMPLES TO RETURN READ SUPPORT FOR ALL MUTATIONS
unique_muts=generate_unique_mutation_list(path_muts)
grep_command_list=[]
for i, row in unique_muts.iterrows():
    chrom = row["Chrom"]
    pos = row["Position"]
    for sample in sample_list["Patient_id"]:
        path_mpileup = os.path.join(dir_mpileups, sample + ".mpileup")
        path_grepped_output = os.path.join(dir_filtered_mpileups, f"{sample}_{chrom}_{pos}.tsv")
        
        # Check if file does not exist or has size 0
        if not os.path.isfile(path_grepped_output) or os.path.getsize(path_grepped_output) == 0:
            grep_command = f"grep {pos} {path_mpileup} > {path_grepped_output}"
            grep_command_list.append(grep_command)

run_commands_in_batches(grep_command_list, batch_size=200)

# STEP 2. GO THROUGH GREPPED OUTPUTS FROM MPILEUP AND START BUILDING THE SAMPLE/MUTATION MATRIX, 
# PAYING ATTENTION TO WHETHER THE MUTATION WAS CALLED OR NOT IN A GIVEN SAMPLE.

# Initialize an empty dictionary to store the data
mut_sample_matrix = {}

for i, row in unique_muts.iterrows():
    chrom = row["Chrom"]
    pos = str(row["Position"])
    mut_type = row["Type"]
    alt = row["Alt"]
    
    # if mut_type=="Insertion":
    #     print(mut_type)
    #     break
        
    if pd.isna(row["Protein_annotation"]):
        annot="splice"
    else:
        annot=row["Protein_annotation"]
    mut_key = row["Gene"] + " " + annot + " " + chrom + ":" + pos  # How we refer to this mutation in the matrix
    
    # Initialize a dictionary to store VAFs for the current mutation across all samples
    sample_vafs = {}
    
    for idx, sample in enumerate(sample_list["Patient_id"], start=1):
        # Construct the path to the filtered mpileup file
        path_filtered_mpileup = os.path.join(dir_filtered_mpileups, f"{sample}_{chrom}_{pos}.tsv")
        
        # Retrieve depth, alt reads, and VAF
        try:
            result = return_depth_and_n_altered_from_filtered_mpileup(path_filtered_mpileup, chrom, pos, mut_type, alt)
            mutation_called=check_if_mutation_was_called(sample.replace("_", "-"), chrom, pos, row["Protein_annotation"], alt, path_muts)
            n_alt=result["n_alt"]
            depth=result["depth"]
            if mutation_called:
                sample_vafs[sample] = f"*{n_alt}/{depth}"
            else:
                sample_vafs[sample] = f"{n_alt}/{depth}"
        
        except Exception as e:
            # Handle missing files or errors by assigning NaN
            sample_vafs[sample] = float("nan")
            
        # Print progress message
    print(f"Processing {mut_key}")
    
    # Add the sample VAFs for this mutation to the main dictionary
    mut_sample_matrix[mut_key] = sample_vafs

# Convert the dictionary to a DataFrame
mut_sample_df = pd.DataFrame.from_dict(mut_sample_matrix, orient="index")

# Reset index to make mutation IDs a column
mut_sample_df = mut_sample_df.reset_index()
mut_sample_df.columns = ["Mutation_ID"] + list(mut_sample_df.columns[1:])

# Display the resulting DataFrame
print(mut_sample_df)


mut_sample_df.to_csv("/groups/wyattgrp/users/amunzur/ironman_ch/results/variant_calling/ironman_mutation_sample_matrix_TEST.csv", index=False)