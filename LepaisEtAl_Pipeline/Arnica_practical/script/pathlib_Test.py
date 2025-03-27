from pathlib import Path
import re

def extract_names_from_folder(folder_path):
    folder = Path(folder_path)
    extracted_names = []
    pattern = re.compile(r"(.*_S\d+_[AB])_.*\.csv$")  

    for file in folder.rglob("*.csv"):  # Recursively find all CSV files
        match = pattern.match(file.name)
        if match:
            extracted_names.append(match.group(1))

    return extracted_names

# Example usage
folder_path = '../filtered_results/filtered_tssvResults' # Change this to your actual folder path
names_list = extract_names_from_folder(folder_path)
print(names_list)