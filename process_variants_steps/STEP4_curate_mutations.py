import numpy as np
import pandas as pd 
import os
import re
import argparse

"""
After curating mutations by going through IGV snapshots, this script detects which mutations have been removed, and then adds them to the blacklist.
"""

def do_misc_filtering_function(df):
    """
    Performs additional filtering by removing variants based on ratio between cfDNA and WBC VAF,
    and removed mutations in LPAR6.
    """
    df = df[df["tumor_wbc_vaf_ratio"] < 20]
    df = df[df["Gene"] != "LPAR6"]
    
    df = df.reset_index(drop = True)
    
    return(df)

def curate(PATH_muts, DIR_curated_screenshots, path_to_keep, path_to_exclude, mut_type, do_misc_filtering): 
    
    # sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])[["Patient_id", "Diagnosis"]].drop_duplicates()
    all_muts = pd.read_csv(PATH_muts)
    # all_muts.loc[all_muts["Diagnosis"] == "Kidney During treatment", "Diagnosis"] = "Kidney"
    # all_muts = add_timepoint(all_muts, PATH_sample_information)
    
    # subset to pts
    # pts = sample_info["Patient_id"].str.replace("GU-", "").unique()
    # all_muts = all_muts[all_muts["Patient_id"].isin(pts)]
    
    if "Sample_name_t" in all_muts.columns:
        all_muts["IGV_screenshot_name"] = all_muts.apply(lambda row: "_".join(map(str, row[['Gene', 'Protein_annotation', 'Chrom', 'Position', 'Sample_name_t']])), axis=1) + ".png" # all called muts
    else:
        all_muts["IGV_screenshot_name"] = all_muts.apply(lambda row: "_".join(map(str, row[['Gene', 'Protein_annotation', 'Chrom', 'Position', 'Sample_name']])), axis=1) + ".png" # all called muts
    
    curated_muts = pd.DataFrame({"IGV_status": "keep", "IGV_screenshot_name": os.listdir(DIR_curated_screenshots)}) # muts we are keeping
    merged = all_muts.merge(curated_muts, how = "left")
    
    muts_to_keep = merged[merged["IGV_status"] == "keep"]
    muts_to_exclude = merged[pd.isnull(merged["IGV_status"])]
    
    # i am excluding this one - possible germline
    # idx = muts_to_keep[~(muts_to_keep["Gene"] == "FLT3") & (muts_to_keep["Sample_name_t"] == "GU-20-305_cfDNA-Baseline-IDT-2021Mar26")]
    muts_to_keep.loc[muts_to_keep["Gene"].str.contains("U2AF"), "Gene"] = "U2AF1"
    
    # exclude oxidative damage patients 
    # chip_oxidative_patients = ["21-307", "21-377"]
    # somatic_oxidative_patients = ["20-313", "21-184", "21-430"]
    # if mut_type == "chip":
    #     muts_to_keep = muts_to_keep[~muts_to_keep["Patient_id"].isin(chip_oxidative_patients)]
    # else:
    #     muts_to_keep = muts_to_keep[~muts_to_keep["Patient_id"].isin(somatic_oxidative_patients)]    
    del muts_to_keep["IGV_status"]
    del muts_to_exclude["IGV_status"]
    
    # Only keep muts present in the patient id list
    # muts_to_keep = sample_info.merge(muts_to_keep, how = "inner")
    
    if do_misc_filtering: 
        muts_to_keep = do_misc_filtering_function(muts_to_keep)
    
    muts_to_exclude.to_csv(path_to_exclude, index = False)
    muts_to_keep.to_csv(path_to_keep, index = False)
    
    print(f"Curated variants saved to {path_to_keep}")
    print(f"Excluded variants saved to {path_to_exclude}")

def main():
    parser = argparse.ArgumentParser(description="Curate mutations after IGV review.")
    parser.add_argument("--mutation_type", required=True, choices=["chip", "somatic"], help="Specify mutation type: chip or somatic.")
    parser.add_argument("--do_misc_filtering", action="store_true", help="Perform additional filtering.")
    parser.add_argument("--DIR_working", required=True, help="")
    parser.add_argument("--DIR_curated_screenshots", required=True, help="")
    args = parser.parse_args()
    
    DIR_working = args.DIR_working
    mutation_type = args.mutation_type.lower()
    do_misc_filtering = args.do_misc_filtering
    DIR_curated_screenshots=args.DIR_curated_screenshots
    
    PATH_muts = os.path.join(DIR_working, f"results/variant_calling/{mutation_type}_SSCS2.csv")
    PATH_retained_variants = os.path.join(DIR_working, f"results/variant_calling/{mutation_type}_SSCS2_curated.csv")
    PATH_excluded_variants = os.path.join(DIR_working, f"resources/validated_variants/{mutation_type}_to_exclude_IGV.csv")
        
    curate(
        PATH_muts=PATH_muts,
        DIR_curated_screenshots=DIR_curated_screenshots,
        path_to_keep=PATH_retained_variants,
        path_to_exclude=PATH_excluded_variants,
        mut_type=mutation_type,
        do_misc_filtering=do_misc_filtering,
    )

if __name__ == "__main__":
    main()

"""
Run example:

python /groups/wyattgrp/users/amunzur/toolkit/STEP4_curate_mutations.py \
    --mutation_type somatic \
    --DIR_working /groups/wyattgrp/users/amunzur/hla_pipeline \
    --DIR_curated_screenshots /groups/wyattgrp/users/amunzur/hla_pipeline/results/figures/IGV_screenshots/non_hla_ctdna_mutations

"""