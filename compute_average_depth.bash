#!/bin/bash

# Usage: 
# /path/to/toolkit/compute_average_depth.bash \
# /path/to/results/metrics/depth/SSCS2 \
# /path/to/results/metrics/averaged_depth.tsv

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <DIR_depth> <PATH_output>"
    exit 1
fi

# Assign command-line arguments to variables
DIR_depth="$1"
PATH_output="$2"

# Ensure output file is empty before writing
rm "$PATH_output"
touch "$PATH_output"

echo -e "Sample_name\tMean_depth\tMedian_depth" > "$PATH_output"


for f in "$DIR_depth"/*.txt; do
    sample_name=$(basename "$f")
    echo "$sample_name"
    
    avg_depth=$(awk '{ total += $3; count++ } END { if (count > 0) print total/count; else print 0 }' "$f")
    median_depth=$(cut -f3 "$f" | sort -n | awk ' { a[i++]=$1; } END { if (i>0) print a[int(i/2)]; else print 0 }')
    
    echo -e "${sample_name%.*}\t$avg_depth\t$median_depth" >> "$PATH_output"
done
