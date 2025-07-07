#!/usr/bin/env python3

# Given a segmentation file calculates the segment copy numbers, saves to a new file.
# Make sure to run this in the base environment: conda activate base
# Per sample diploid levels and ctDNA fractions should be available

# /groups/wyattgrp/users/amunzur/toolkit/calculate_CN_segment.py \
#     --path_segmentation /path/to/ecdna_project/igv_tracks/sample_segments.tsv \
#     --path_output /groups/wyattgrp/users/amunzur/ecdna_project/igv_tracks/sample_segments_CN.tsv \
#     --diploid_level -0.707 \
#     --ctdna_fraction 0.691 

import sys
import pandas as pd
import argparse

def calculate_CN(path_segmentation, diploid_level, ctdna_fraction):
    df = pd.read_csv(path_segmentation, sep = "\t")
    df = df[~df["CHR"].isin(["chrX", "chrY", "chrM"])]
    df["CN"] = df.apply(lambda row:  (2 * (ctdna_fraction + 2**(row["LOGRATIO"] - diploid_level) - 1)) / ctdna_fraction, axis=1)
    df = df[['CHR', 'START', 'END', 'LOGRATIO', 'CN']]
    return(df)

def main():
    # Argument parsing
    # DL and cfDNA fraction must be precomputed.
    parser = argparse.ArgumentParser(description="Calculate raw copy numbers from segmentation file.")
    parser.add_argument("--path_segmentation", required=True, help="Absolute path to the segmentation file")
    parser.add_argument("--path_output", help="Optional, if given saves here.")
    parser.add_argument("--diploid_level", required=True, help="Supply the sample's diploid level")
    parser.add_argument("--ctdna_fraction", required=True, help="Supply the sample's ctdna fraction as an integer, not as a fraction.")
    
    args = parser.parse_args()
    
    if "--help" in sys.argv:
        parser.print_help()
        sys.exit()
    
    if args.path_output is None:
        path_output = args.path_segmentation.replace(".tsv", "_CN.tsv")
    else: 
        path_output = args.path_segmentation.replace(".tsv", "_CN.tsv")
    
    # print(f"path_segmentation='{args.path_segmentation}'")
    # print(f"path_output='{args.path_output}'")
    # print(f"diploid_level='{args.diploid_level}'")
    # print(f"ctdna_fraction='{args.ctdna_fraction}'")
    
    if float(args.ctdna_fraction) > 1:
        raise ValueError("Please provide the ctDNA fraction as an integer, not fraction.")
    elif float(args.ctdna_fraction) == 0:
        raise ValueError("Provided ctDNA fraction is 0.")
    
    df = calculate_CN(args.path_segmentation, float(args.diploid_level), float(args.ctdna_fraction))
    df.to_csv(path_output, index = None, header = None, sep = "\t")
    print(f"Output saved to {path_output}")

if __name__ == "__main__":
    main()