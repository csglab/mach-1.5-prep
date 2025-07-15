#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 input.csv.gz"
    exit 1
fi

input_file="$1"
output_file="${input_file%.csv.gz}.jsonl.zst"

# Check if output file already exists
if [ -f "$output_file" ]; then
    echo "Output file $output_file already exists, skipping..."
    exit 0
fi

zcat "$input_file" | \
    awk -F, 'NR>1 { # Skip header row
        # Create the meta string
        meta = sprintf("transcript_id=%s;gene_id=%s;chr=%s;seq_len=%s;contains_unk=%s", 
                      $1, $2, $3, $5, $6)
        
        # Output JSON format
        printf "{\"text\":\"%s\",\"meta\":\"%s\"}\n", $4, meta
    }' | zstd > "$output_file"

echo "Conversion complete: $output_file created"
