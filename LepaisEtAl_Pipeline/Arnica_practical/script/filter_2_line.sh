#!/bin/bash

#16.01 : I am trying to make a whole script which can do everything. 
#later change : for the tssvResults and stuttermark, two loops needed  

# File containing the names to search for
names_file="../nSSR_LocusList.txt"

# Input folder containing the data files
input_folder="../tssvReports"

# Output folder for the filtered files
 output_folder="../filtered_tssvReports"

# Ensure the output folder exists
rm -rf "$output_folder"
mkdir  "$output_folder"

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
  while IFS= read -r name; do
    # Extract the first two lines for each name and append to the output file
    grep "^$name" "$data_file" | head -n 2 >> "$output_file"
  done < "$names_file"
done

# be careful, there is a nc sample which need to be correct manually.