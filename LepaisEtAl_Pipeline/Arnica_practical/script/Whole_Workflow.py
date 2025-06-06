# Import necessary modules
import csv
import os
import json
from pathlib import Path
import re

# Constants for repeated values
NA = "N/A"
OTHER_SEQUENCE = "other sequence"

# Parameters
folder_path = '../filtered_results/filtered_tssvResults'
loci_list = ["Arm01", "Arm03", "Arm04", "Arm06", "Arm07", "Arm08", "Arm11"]
locus_coverage_file = "../LocusCoverageperIndividual_nSSR_FullLength.csv"
allele_info_csv_filename = '../new_output/AlleleInfo.csv'
alternative_output_csv_filename = '../new_output/Alternative_AlleleInfo.csv'

# Utility Functions
def extract_names_from_folder(folder_path):
    """Extract sample names from the folder."""
    folder = Path(folder_path)
    extracted_names = []
    pattern = re.compile(r"(.*_S\d+_[AB])_.*\.csv$")

    for file in folder.rglob("*.csv"):  # Recursively find all CSV files
        match = pattern.match(file.name)
        if match:
            extracted_names.append(match.group(1))

    return extracted_names

def save_csv_as_matrix(csv_file_path):
    """Load a CSV file into a matrix."""
    with open(csv_file_path, mode='r') as file:
        reader = csv.reader(file)
        matrix = [row for row in reader]
    return matrix

def write_csv(output_filename, sample_list, process_sample_func):
    """Writes processed sample data to a CSV file."""
    with open(output_filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Serial", "Sample", "Loci", "Allele1", "Allele2", "Read_count", "RC_Allele1", "RC_Allele2", "Frequency1", "Frequency2"])
        
        serial = 1
        for sample_name in sample_list:
            data = process_sample_func(sample_name)
            
            if isinstance(data, str):  # Handle error messages
                print(data)
                continue
            
            for loci, candidates in data.items():
                if isinstance(candidates, str):  # Handle insufficient loci data
                    writer.writerow([serial, sample_name, loci, NA, NA, NA, NA, NA, NA, NA])
                else:
                    allele1_details = candidates.get("first_candidate", [NA, NA, NA, NA, NA])
                    allele2_details = candidates.get("second_candidate", [OTHER_SEQUENCE, len(OTHER_SEQUENCE), 1, NA, NA])
                    
                    writer.writerow([
                        serial, sample_name, loci,
                        allele1_details[0], allele2_details[0],
                        allele1_details[2], allele1_details[3], allele2_details[3],
                        allele1_details[4], allele2_details[4]
                    ])
                serial += 1
        print(f"CSV file '{output_filename}' created successfully!")

# Core Processing Functions
def process_sample(sample_name):
    """Process a single sample and extract allele information."""
    def search_folder(folder_path, keyword):
        for file_name in os.listdir(folder_path):
            if keyword in file_name and file_name.endswith('.csv'):
                file_path = os.path.join(folder_path, file_name)
                with open(file_path, 'r') as file:
                    reader = csv.reader(file, delimiter='\t')
                    return [row for row in reader] 
        return None

    sample_file = search_folder(folder_path, sample_name)
    if sample_file is None:
        return f"Sample '{sample_name}' not found in folder."

    def open_file(file_path):
        with open(file_path, newline='') as file:
            reader = csv.reader(file, delimiter='\t')
            return [row for row in reader]

    score_matrix = open_file(locus_coverage_file)

    def search_column(matrix, word):
        if word in matrix[0]:
            col_index = matrix[0].index(word)
            return [row[col_index] for row in matrix]
        return None

    score_column = search_column(score_matrix, sample_name)
    if not score_column or len(score_column) < 2 or score_column[1] == '0':
        return f"Invalid score column data for {sample_name}"

    def delete_lines_with_keyword(matrix, keyword="Other sequence"):
        return [row for row in matrix if all(keyword not in element for element in row)]

    cleaned_sample_matrix = delete_lines_with_keyword(sample_file)

    def select_loci(sample_file, loci): 
        loci_rows = [row for row in sample_file if any(loci in element for element in row)]
        return [(row[1], row[2]) for row in loci_rows[:2]]  # Extract first two candidates

    def PercentNumber(sample_name, loci):
        try:
            get_locus_index = lambda locus_name: loci_list.index(locus_name) + 1
            loci_data = select_loci(cleaned_sample_matrix, loci)

            if len(loci_data) == 0:
                return {loci: "Insufficient loci data"}

            locus_index = get_locus_index(loci)
            if locus_index >= len(score_column) or not score_column[locus_index].isdigit():
                return {loci: "Invalid score column data"}

            locus_score = int(score_column[locus_index])
            if locus_score == 0:
                locus_score = 1

            candidate1 = loci_data[0]
            read_count_1 = int(candidate1[1]) if len(candidate1) > 1 and candidate1[1].isdigit() else NA
            first_percent = round(read_count_1 / locus_score, 3) if read_count_1 != NA else NA

            if len(loci_data) == 1:
                return {
                    loci: {
                        "first_candidate": [candidate1[0], len(candidate1[0]), locus_score, read_count_1, first_percent]
                    }
                }

            candidate2 = loci_data[1]
            read_count_2 = int(candidate2[1]) if len(candidate2) > 1 and candidate2[1].isdigit() else NA
            second_percent = round(read_count_2 / locus_score, 3) if read_count_2 != NA else NA

            return {
                loci: {
                    "first_candidate": [candidate1[0], len(candidate1[0]), locus_score, read_count_1, first_percent],
                    "second_candidate": [candidate2[0], len(candidate2[0]), locus_score, read_count_2, second_percent],
                }
            }
        except (ValueError, IndexError) as e:
            return {loci: f"Error processing loci data: {str(e)}"}

    results = {}
    for loci in loci_list:
        results.update(PercentNumber(sample_name, loci))
    return results

# Main Execution
if __name__ == "__main__":
    # Extract sample names
    sample_list = extract_names_from_folder(folder_path)

    # Write AlleleInfo.csv
    write_csv(allele_info_csv_filename, sample_list, process_sample)

    # Write Alternative_AlleleInfo.csv
    write_csv(alternative_output_csv_filename, sample_list, process_sample)

    # Load AlleleInfo.csv as a matrix
    allele_info_matrix = save_csv_as_matrix(allele_info_csv_filename)
    print("Matrix loaded from AlleleInfo.csv:")
    for row in allele_info_matrix:
        print(row)
