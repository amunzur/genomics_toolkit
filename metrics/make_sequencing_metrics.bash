#!/bin/bash

PATH_hs_metrics="$1"
PATH_output="$2"
DIR_raw_read_counts="$3"
# DIR_markdup_perc="$4"

get_picard_field() {
    local file="$1"
    local field="$2"
    
    awk -v target="$field" '
        BEGIN { FS=OFS="\t" }
        # Line 6 is where header is
        NR == 7 {
            # Find correct col number by column name
            for (i = 1; i <= NF; i++) {
                if ($i == target) {
                    col = i
                    break
                }
            }
            next
        }
        NR == 8 {
            # Print the value if present, otherwise print NA
            if (col) {
                if ($col == "" || $col == "-") print "NA"
                else print $col
            } else {
                # Column missing entirely
                print "NA"
            }
        }
    ' "$file"
}


# Ensure output file is empty before writing
mkdir -p "$(dirname "$PATH_output")"

# Tab delim
: > "$PATH_output"
echo -e "Sample_name\tMedian_depth\tOn_target_percentage\tNumber_of_reads(M)\tDuplicate_percentage" > "$PATH_output"

echo "Looking in: $PATH_hs_metrics, matching files:"
ls "$PATH_hs_metrics"/*.HS_metrics

for f in "$PATH_hs_metrics"/*.HS_metrics; do
    sample_name=$(basename "$f" .HS_metrics)
    echo "$sample_name"

    # --- here Parse PICARD METRICS ---
    MEDIAN_TARGET_COVERAGE=$(get_picard_field "$f" "MEDIAN_TARGET_COVERAGE")
    PCT_SELECTED_BASES=$(get_picard_field "$f" "PCT_SELECTED_BASES")
    TOTAL_READS_FROM_HS=$(get_picard_field "$f" "TOTAL_READS")

    if [ "$TOTAL_READS_FROM_HS" = "NA" ]; then
        TOTAL_READS_FROM_HS=$(get_picard_field "$f" "PF_READS")
    fi
    # ---------------------------------------

    # --- Prefer raw read counts if available ---
    if [ -z "$DIR_raw_read_counts" ]; then
        TOTAL_READS=$(printf "%.2f" "$(echo "$TOTAL_READS_FROM_HS" / 1e6 | bc -l)") # convert to megabases
    else
        path_read_counts="${DIR_raw_read_counts}/${sample_name}.txt"
        if [ -f "$path_read_counts" ]; then
            TOTAL_READS=$(awk '{printf "%.2f", $2 / 1e6}' "$path_read_counts") # convert to megabases
        else
            echo "Warning: $path_read_counts not found, using HS metrics instead."
            TOTAL_READS=$(printf "%.2f" "$(echo "$TOTAL_READS_FROM_HS" / 1e6 | bc -l)")
        fi
    fi

    echo -e "${sample_name}\t${MEDIAN_TARGET_COVERAGE}\t${PCT_SELECTED_BASES}\t${TOTAL_READS}" >> "$PATH_output"
done
