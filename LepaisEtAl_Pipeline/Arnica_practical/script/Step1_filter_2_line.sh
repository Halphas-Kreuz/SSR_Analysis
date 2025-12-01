#!/bin/bash

# Read the number of lines from the first command-line argument ($1)
# If no argument is provided (${1:-2}), default to 2.
LINES_TO_KEEP=${1:-2}

echo "Filtering to keep the top $LINES_TO_KEEP lines for each locus."

filter_files() {
  local input_folder="$1"
  local output_folder="../filtered_results/$(basename "$2")"
  # Assign the third argument to a local variable
  local lines_to_keep="$3" 
  local names_file="../nSSR_LocusList.txt"

  # Ensure the output folder exists
  rm -rf "$output_folder"
  mkdir -p "$output_folder"

  # Process each file in the input folder
  for data_file in "$input_folder"/*.txt; do
    # Get the base name of the file (e.g., "A.txt" -> "A")
    base_name=$(basename "$data_file" .txt)
    base_name=${base_name%_tssv}
    output_file="$output_folder/${base_name}.txt"

    # Clear the output file if it exists, to avoid appending
    > "$output_file"

    # Process each name in the names file
    while IFS= read -r name || [[ -n "$name" ]]; do # now it does not ignore the last line
      # Extract the first N lines (using the new variable)
      grep "^$name" "$data_file" | head -n "$lines_to_keep" >> "$output_file"
    done < "$names_file"
  done
}
# automatically convert everything 
convert_txt_to_csv() {
  local folder="../filtered_results/$(basename "$1")"
  for file in "$folder"/*.txt; do
    mv "$file" "${file%.txt}.csv"
  done
}


filter_files "../tssvResults" "../filtered_tssvResults" "$LINES_TO_KEEP"
filter_files "../stuttermark" "../filtered_stuttermark" "$LINES_TO_KEEP"
filter_files "../tssvReports" "../filtered_tssvReports" "$LINES_TO_KEEP"

convert_txt_to_csv "../filtered_tssvResults"
convert_txt_to_csv "../filtered_stuttermark"
convert_txt_to_csv "../filtered_tssvReports"

