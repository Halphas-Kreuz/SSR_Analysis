import csv
import os
import json
from pathlib import Path
import re
import sys  # Import sys to read command-line arguments

# --- Parameter Part ---

# Get 'n' from command line, default to 2
try:
    # sys.argv[0] is the script name, sys.argv[1] is the first arg
    N_LINES_TO_PROCESS = int(sys.argv[1])
except (IndexError, ValueError):
    N_LINES_TO_PROCESS = 2 # Default to 2 if no argument is given

MIN_SEQUENCE_LENGTH = 8   
MIN_FREQUENCY_RATIO = 0.1 # manuel changeable
print(f"✅ Processing the top {N_LINES_TO_PROCESS} candidates for each locus.")

folder_path = '../filtered_results/filtered_tssvResults'
with open ('../nSSR_LocusList.txt', 'r') as f:
    loci_list = [line.strip() for line in f if line.strip()]
locus_coverage_file = "../LocusCoverageperIndividual_nSSR_FullLength.txt"

#recursively extract all the data name from the folder
def extract_names_from_folder(folder_path):
    folder = Path(folder_path)
    extracted_names = []
    pattern = re.compile(r"^(.*)\.csv$")

    for file in folder.rglob("*.csv"):  # Recursively find all CSV files
        match = pattern.match(file.name)
        if match:
            extracted_names.append(match.group(1))

    return sorted(extracted_names)

sample_list = extract_names_from_folder(folder_path)


def process_sample(extrande, num_lines_to_keep):
    if extrande.endswith('_tssv_filtered'):
        extrande = extrande.replace('_tssv_filtered', '') 
    def search_folder(folder_path, keyword):
        for file_name in os.listdir(folder_path):
            if keyword in file_name and file_name.endswith('.csv'):
                file_path = os.path.join(folder_path, file_name)
                with open(file_path, 'r') as file:
                    reader = csv.reader(file, delimiter='\t')
                    return [row for row in reader] 
        return None
        

    sample_file = search_folder(folder_path, extrande)
    if sample_file is None:
        return f"Sample '{extrande}' not found in folder."

    score_matrix = []
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

    score_column = search_column(score_matrix, extrande)
    if not score_column or len(score_column) < 2:
        return f"Invalid score column data for {extrande}"
    
    def select_loci(sample_file, loci, n): 
        loci_rows = [row for row in sample_file if row and row[0] == loci]
        # Use the 'n' parameter here instead of hardcoded :2
        return [(row[1], row[2]) for row in loci_rows[:n]]

    def PercentNumber(sample_name, loci):
        try:
            get_locus_index = lambda locus_name: loci_list.index(locus_name) + 1

            # Pass the number of lines to select_loci
            loci_data = select_loci(sample_file, loci, num_lines_to_keep)

            if len(loci_data) == 0:
                return {loci: "Insufficient loci data"}

            locus_index = get_locus_index(loci)
            if locus_index >= len(score_column) or not score_column[locus_index].strip().isdigit():
                return {loci: "Invalid score column data"}

            locus_score = int(score_column[locus_index])
            if locus_score == 0:
                locus_score = 1

            # --- Process ALL candidates found (up to n) ---
            candidates_list = []
            for candidate in loci_data:
                try:
                    percent = round(int(candidate[1]) / locus_score, 3)
                    seq_len = len(candidate[0])
                    if seq_len >= MIN_SEQUENCE_LENGTH and percent >= MIN_FREQUENCY_RATIO:                      
                        candidates_list.append((candidate[0], seq_len, percent))
                    # Store as a tuple: (sequence, length, frequency)
                except ValueError:
                    # Handle case where candidate[1] is not a number
                    candidates_list.append((candidate[0], len(candidate[0]), "N/A"))
            
            # Return the whole list of candidates
            return {loci: candidates_list}

        except (ValueError, IndexError) as e:
            return {loci: f"Error processing loci data: {str(e)}"}
    
    results = {}
    for loci in loci_list:
        results.update(PercentNumber(extrande, loci))
        
    return results



output_csv_filename = '../new_output/Step2_AlleleInfo.csv'
output_path = Path(output_csv_filename)
output_path.parent.mkdir(parents=True, exist_ok=True)

with open(output_csv_filename, mode='w', newline='') as file:
    writer = csv.writer(file)

    # --- (MODIFICATION 1) ---
    # --- Write DYNAMIC grouped header ---
    header = ["Serial", "Sample", "Loci"]
    # First, add all Allele columns
    for i in range(1, N_LINES_TO_PROCESS + 1):
        header.append(f"Allele{i}")
    # Second, add all Frequency columns
    for i in range(1, N_LINES_TO_PROCESS + 1):
        header.append(f"Frequency{i}")
    writer.writerow(header)
    # --- (END MODIFICATION 1) ---

    serial = 1  # Serial number starts at 1

    for i in range(len(sample_list)):
        sample_name = sample_list[i]
        # Pass N_LINES_TO_PROCESS to the function
        data = process_sample(sample_name, N_LINES_TO_PROCESS)

        if isinstance(data, str):  # Handle error messages
            print(data)
            continue

        for loci, candidates in data.items():
            # Base row
            row_data = [serial, sample_name, loci]
            
            alleles_list = []
            freqs_list = []

            if isinstance(candidates, str):  # Handle insufficient loci data or errors
                # Fill all allele/freq columns with N/A
                alleles_list = ["N/A"] * N_LINES_TO_PROCESS
                freqs_list = ["N/A"] * N_LINES_TO_PROCESS
            else:
                # 'candidates' is now a list of tuples: (seq, length, freq)
                for candidate_tuple in candidates:
                    alleles_list.append(candidate_tuple[0]) # Allele Seq
                    freqs_list.append(candidate_tuple[2])   # Frequency
                
                # Pad with "N/A" if fewer than N_LINES candidates were found
                while len(alleles_list) < N_LINES_TO_PROCESS:
                    alleles_list.append("N/A")
                    freqs_list.append("N/A")

            # --- (MODIFICATION 2) ---
            # Add all alleles first
            row_data.extend(alleles_list)
            # Then add all frequencies
            row_data.extend(freqs_list)
            # --- (END MODIFICATION 2) ---
            
            # Write the complete dynamic row
            writer.writerow(row_data)
            serial += 1  # Increment serial number

    print(f"CSV file '{output_csv_filename}' created successfully!")