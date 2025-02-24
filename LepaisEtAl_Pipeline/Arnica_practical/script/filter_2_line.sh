#!/bin/bash

filter_files() {
  local names_file="../nSSR_LocusList.txt"
  local input_folder="$1"
  local output_folder="../filtered_results/$(basename "$2")"

  # Ensure the output folder exists
  rm -rf "$output_folder"
  mkdir -p "$output_folder"

  # Process each file in the input folder
  for data_file in "$input_folder"/*.txt; do
    # Get the base name of the file (e.g., "A.txt" -> "A")
    base_name=$(basename "$data_file" .txt)
    new_name=$(echo "$base_name" | sed -E 's/.*AM(.*)/AM\1/')
    
    # Create the output file name with "_filtered" appended
    output_file="$output_folder/${new_name}_filtered.txt"

    # Clear the output file if it exists, to avoid appending
    > "$output_file"

    # Process each name in the names file
    while IFS= read -r name || [[ -n "$name" ]]; do # now it does not ignore the last line
      # Extract the first two lines for each name and append to the output file
      grep "^$name" "$data_file" | head -n 2 >> "$output_file"
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


# Hardcoded folder paths
filter_files "../tssvResults" "../filtered_tssvResults"
filter_files "../stuttermark" "../filtered_stuttermark"
filter_files "../tssvReports" "../filtered_tssvReports"

convert_txt_to_csv "../filtered_tssvResults"
convert_txt_to_csv "../filtered_stuttermark"
convert_txt_to_csv "../filtered_tssvReports"
