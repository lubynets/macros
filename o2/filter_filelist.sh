#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <original_filelist> <selected_runs>" >&2
    exit 1
fi

filelist="$1"
selected_runs="$2"
output_file="filelist_filtered.txt"

if [[ ! -f "$filelist" ]]; then
    echo "Cannot find filelist: $filelist" >&2
    exit 1
fi

if [[ ! -f "$selected_runs" ]]; then
    echo "Cannot find selected runs: $selected_runs" >&2
    exit 1
fi

# Read selected runs into a set
declare -A runs

while IFS= read -r run || [[ -n "$run" ]]; do
    run="${run%$'\r'}"
    [[ -z "$run" ]] && continue
    runs["$run"]=1
done < "$selected_runs"

# Create output file
: > "$output_file"

n_total=0
n_selected=0

while IFS= read -r file_path || [[ -n "$file_path" ]]; do
    file_path="${file_path%$'\r'}"
    [[ -z "$file_path" ]] && continue

    ((n_total+=1))

    directory="${file_path%/*}"
    runnumber_file="$directory/runnumber.txt"

    if [[ ! -f "$runnumber_file" ]]; then
        echo "Warning: No runnumber.txt found in $directory" >&2
        continue
    fi

    # Read run numbers, including comma-separated entries
    run_numbers=$(<"$runnumber_file")
    run_numbers="${run_numbers//,/ }"

    for run in $run_numbers; do
        if [[ -v runs[$run] ]]; then
            printf '%s\n' "$file_path" >> "$output_file"
            ((n_selected+=1))
            break
        fi
    done

done < "$filelist"

echo "Total entries: $n_total"
echo "Selected entries: $n_selected"
echo "Output file: $output_file"
