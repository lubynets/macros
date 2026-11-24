#!/usr/bin/env bash

set -euo pipefail

# Check arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_file> <path_prefix>" >&2
    exit 1
fi

input_file="$1"
prefix="${2%/}"

if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file not found: $input_file" >&2
    exit 1
fi

# Output file in the directory from which the script is called
output_file="$PWD/runnumbers.txt"

# Clear/create the output file
: > "$output_file"

# Process each line in the input file
while IFS= read -r file_path || [[ -n "$file_path" ]]; do

    # Remove possible Windows carriage return
    file_path="${file_path%$'\r'}"

    # Replace the original filename with metrics_summary.json
    directory="${file_path%/*}"
    metrics_file="$directory/metrics_summary.json"
    runnumber_file="$directory/runnumber.txt"

    if [[ ! -f "$metrics_file" ]]; then
        result="No metrics_summary.json found"
        echo "$result" >&2
    else
        # Extract unique 6-digit numbers after the specified prefix
        result=$(
            grep -oE \
                "${prefix}/[^/]+/[0-9]{6}/" \
                "$metrics_file" |
            grep -oE '[0-9]{6}/$' |
            tr -d '/' |
            sort -u |
            paste -sd ',' - |
            sed 's/,/, /g'
        )

        # Handle files without matching run numbers
        if [[ -z "$result" ]]; then
            result="No run numbers found"
        fi

        # Write result next to metrics_summary.json
        printf '%s\n' "$result" > "$runnumber_file"
    fi

    # Append result to the global output file
    printf '%s\n' "$result" >> "$output_file"

done < "$input_file"

echo "Done. Results saved to $output_file"
