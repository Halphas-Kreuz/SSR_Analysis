#!/bin/bash

# Folder containing the files
input_folder="../tssvReports"

# Loop through all files in the folder
for file in "$input_folder"/*; do
    # Extract the base name of the file
    base_name=$(basename "$file")

    # Extract the part after "AM"
    new_name=$(echo "$base_name" | sed -E 's/.*AM(.*)/AM\1/')

    # Rename the file
    mv "$file" "$input_folder/$new_name"
    echo "Renamed $file to $new_name"
done
