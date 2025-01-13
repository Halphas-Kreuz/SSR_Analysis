#!/bin/bash

# Input folder containing the data files
input_folder="../stuttermark"

# File containing the names to search for
names_file="../nSSR_LocusList.txt"

# Output folder for the filtered files
output_folder="../filtered_2_line"

# Ensure the output folder exists
rm -rf "$output_folder"
mkdir  "$output_folder"

# Process each file in the input folder
for data_file in "$input_folder"/*.txt; do
  # Get the base name of the file (e.g., "A.txt" -> "A")
  base_name=$(basename "$data_file" .txt)
  
  # Create the output file name with "_filtered" appended
  output_file="$output_folder/${base_name}_filtered.txt"

  # Clear the output file if it exists, to avoid appending
  > "$output_file"

  # Process each name in the names file
  while IFS= read -r name; do
    # Extract the first two lines for each name and append to the output file
    grep "^$name" "$data_file" | head -n 2 >> "$output_file"
  done < "$names_file"
done


#usage : ./filter_folder_files.sh input_folder names_file.txt output_folder