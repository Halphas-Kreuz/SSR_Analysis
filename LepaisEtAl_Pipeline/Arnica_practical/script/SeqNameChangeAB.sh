#!/bin/bash

# Specify the directory to operate on
directory="$1"

# Check if the directory is provided and exists
if [ -z "$directory" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

if [ ! -d "$directory" ]; then
    echo "Error: Directory '$directory' does not exist."
    exit 1
fi

# Find and rename files
find "$directory" -type f | while read -r file; do
    # Compute the new file name
    new_name=$(echo "$file" | sed 's/_R1/_A/g; s/_R2/_B/g')
    
    # Rename the file if the name has changed
    if [ "$file" != "$new_name" ]; then
        mv -v "$file" "$new_name"
    fi
done

echo "File names updated successfully."
