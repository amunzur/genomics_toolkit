import os 
import re
import subprocess
import pandas as pd
import argparse
import sys

'''
Rename FastQ files based on barcode-to-sample ID mapping.

This script:
- Scans a directory for `.fastq.gz` or `.fq.gz` files.
- Matches barcodes in filenames to sequencing IDs using a provided barcode sheet.
- Generates a bash script to rename and move the files to a specified directory.
- Optionally reverse-complements barcodes or runs in dry-run mode.

Usage:
python /path/to/toolkit/rename_fastq.py \
--input_dir /path/to/results/data/fastq/new \
--barcode_sheet /path/to/resources/sequencing/sequencing_sheet.csv \
--output_dir /path/to/results/data/fastq/raw \
--script_path /path/to/workflow/batch_scripts/rename_fastq/rename_dec27_run.bash \
--dry_run

Barcode sheet can be a tsv or csv file. It must contain these three colums: Sequencing_ID, Barcode, GSC Pool ID
- Sequencing_ID will be the new sample name.
- Barcode is GSC assigned sequencing barcode.
- GSC Pool ID: This is the pool identifier, usually something like PX3713. The raw sample names also have this identifier at the beginning of the file name.

'''

def find_fastq_files(DIR_data): 
    '''
    Given a directory, find all files ending with either fq.gz or fastq.gz. Save their absolute paths to a list.
    '''
    file_names = []
    for root, dirs, files in os.walk(DIR_data):
        for file in files: 
            if file.endswith(("fastq.gz", "fq.gz")):
                file_names.append(os.path.join(root, file))
    return file_names

def find_barcode(str_list):
    '''
    Given a file path to a fastq file with the barcode somewhere in its path, return the barcode.
    '''
    letters = ["A", "C", "G", "T", "-"]
    barcode = [x for x in str_list if all(letter in letters for letter in x)][0]
    return(barcode)

def get_reverse_compliment(some_fastq):
    '''
    Sometimes the barcode includes a portion that is the reverse compliment. Check sample sheet and determine if this is the case. 
    This function identifies the reverse compliment of the second part of the barcode and returns a complete new file name as string.
    '''
    base_name = os.path.basename(some_fastq)
    seq_id = base_name.split("_")[1].split("-")[1] # Extract the part you need for reverse complimenting, usually i5 (second index after -) needs to be reverse complimented.
    id_reverse = seq_id[::-1] # reverse
    mydict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'} # dict to take the reverse 
    id_reverse_compliment = []
    
    for char in  id_reverse:
        for pattern, repl in mydict.items():
            if char == pattern: # only continue if the pattern and the character we are at match
                s = re.sub(pattern, repl, char)
                id_reverse_compliment.append(s)
    
    id_reverse_compliment = "".join(id_reverse_compliment)
    
    # replace the id in the fastq name with the reverse compliment 
    base_name_modified = base_name.replace(seq_id, id_reverse_compliment)
    reverse_compliment_barcode=base_name_modified.split("_")[1]
    return(reverse_compliment_barcode)

def get_sequencing_id(some_fastq, barcode_sheet, reverse_compliment = False):
    '''
    Given a file path to a fastq file with barcodes, scan the barcodes sheet to find its sequencing ID and return it as a string.
    barcode_sheet could be a csv or tsv file.
    '''
    if barcode_sheet.endswith(".tsv"):
        df_barcodes = pd.read_csv(barcode_sheet, sep = "\t")
    elif barcode_sheet.endswith(".csv"): 
        df_barcodes = pd.read_csv(barcode_sheet)
    
    # Sometimes if the whole lane is used or something like that the file naming is slightly different.
    if not some_fastq.endswith("merge.fastq.gz"):
        pool = "_".join(os.path.basename(some_fastq).split("_")[:2])
        read = os.path.basename(some_fastq).split("_")[2]
        suffix = f"_{read}.fq.gz"
    else:
        pool = os.path.basename(some_fastq).split("_")[0]
        suffix = "_1.fq.gz" if "_1_" in os.path.basename(some_fastq) else "_2.fq.gz"
    
    if reverse_compliment: 
        barcode = get_reverse_compliment(some_fastq)
    else:
        barcode = find_barcode(some_fastq.split("_"))
    
    seq_id_list = df_barcodes[(df_barcodes["Barcode"] == barcode) & (df_barcodes["GSC Pool ID"] == pool)]["Sequencing_ID"]
    
    if len(seq_id_list) > 1: 
        raise ValueError(f"The barcode {barcode} from {os.path.basename(some_fastq)} matched to more than one sample. Specify pool.")
    elif len(seq_id_list) == 0:
        raise ValueError(f"The barcode {barcode} from {os.path.basename(some_fastq)} didn't match any sample.")
    else:
        seq_id = seq_id_list.values[0]
    
    # Add suffix for R1 or R2
    seq_id += suffix
    
    return {'FastQ_path': some_fastq, 'barcode': barcode, 'seq_id': seq_id}

def write_script(file_names, DIR_output, PATH_script):
    '''
    Writes a bash script to rename move the fastq files into a new location. You would then run that script to rename the fastqs.
    file_names: the list of dictionaries, outputted by the get_sequencing_id function.
    DIR_output: the dir to save the renamed fastq files.
    PATH_script: The bash script that will be generated to rename and move files.
    '''
    if not os.path.exists(os.path.dirname(PATH_script)):
        os.makedirs(os.path.dirname(PATH_script))
    
    if os.path.exists(PATH_script):
        try:
            os.remove(PATH_script)
        except OSError as e:
            print(f"Error: {e.filename} - {e.strerror}")            
    for my_dict in file_names:
        source=my_dict["FastQ_path"]
        destination=os.path.join(DIR_output, my_dict["seq_id"])
        command="mv " + source + " " + destination
        with open(PATH_script, 'a') as file:
            file.write(command)
            file.write("\n")    

def execute_script(PATH_script, dry_run=False):
    if dry_run:
        print(f"Dry run: the following bash script would be executed:\n")
        with open(PATH_script, 'r') as file:
            print(file.read())
    else:
        try:
            subprocess.run(f"bash {PATH_script}", check=True, shell=True)
            print(f"Bash script executed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error executing bash script: {e}")

# Main function to parse arguments and run the renaming process
def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Rename FastQ files based on barcode matching")
    parser.add_argument("--input_dir", required=True, help="Directory containing the FastQ files")
    parser.add_argument("--barcode_sheet", required=True, help="CSV file containing barcode mapping")
    parser.add_argument("--output_dir", required=True, help="Directory where renamed FastQ files will be saved")
    parser.add_argument("--script_path", required=True, help="Path to save the generated bash script")
    parser.add_argument("--reverse_compliment", action='store_true', help="Flag to reverse complement barcodes")
    parser.add_argument("--dry_run", action='store_true', help="Flag to show the bash script without executing it")
    
    args = parser.parse_args()
    
    print(f"input_dir='{args.input_dir}'")
    print(f"barcode_sheet='{args.barcode_sheet}'")
    print(f"output_dir='{args.output_dir}'")
    print(f"script_path='{args.script_path}'")
    print(f"reverse_compliment={args.reverse_compliment}")
    print(f"dry_run={args.dry_run}")
    
    # If user asks for help, print it and exit
    if "--help" in sys.argv:
        print_help()
        sys.exit()
    
    # Find all FastQ files in the input directory
    file_names = find_fastq_files(args.input_dir)
    
    renamed_files = []
    for file in file_names:
        try:
            file_info = get_sequencing_id(file, args.barcode_sheet, args.reverse_compliment)
            renamed_files.append(file_info)
        except ValueError as e:
            print(f"Error: {e}")
    
    # Write bash script to rename and move files
    write_script(renamed_files, args.output_dir, args.script_path)
    print(f"Bash script has been written to: {args.script_path}")
    
    # Execute the bash script if not a dry run
    execute_script(args.script_path, dry_run=args.dry_run)

if __name__ == "__main__":
    main()
