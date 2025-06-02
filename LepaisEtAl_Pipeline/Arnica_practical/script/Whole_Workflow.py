#import module please insert here 
import csv
import os
import json
from pathlib import Path
import re


# Table 1 part 
# Originally from test.py
# Input: whole filter_stuttermark folder; LocusCoverage_File; Loci list 
# Output: Table 1 in the following format (actual name is AlleleInfo.csv)
# [Serial, Sample, Loci, Allele1, Allele2, Frequency1, Frequency2] 

#parameter part 

folder_path = '../filtered_results/filtered_tssvResults'
# sample_name = 'AM30-09_S509_A' # this is only as test , later it will be applied 
loci_list = ["Arm01", "Arm03", "Arm04", "Arm06", "Arm07", "Arm08", "Arm11"]  
locus_coverage_file = "../LocusCoverageperIndividual_nSSR_FullLength.csv"
AlleleInformation = "../AlleleInformationFile_nSSR_FullLength_ParameterSet2_sa70_sb10_m10_n20.csv"

#recursively extract all the data name from the folder 

def extract_names_from_folder(folder_path):
    folder = Path(folder_path)
    extracted_names = []
    pattern = re.compile(r"(.*_S\d+_[AB])_.*\.csv$")  

    for file in folder.rglob("*.csv"):  # Recursively find all CSV files
        match = pattern.match(file.name)
        if match:
            extracted_names.append(match.group(1))

    return extracted_names

sample_list = extract_names_from_folder(folder_path)


def process_sample(extrande):
    # still keep as filter or not ? 
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
    #if the count is < 0.25, then drop it 
    if not score_column or len(score_column) < 2 or score_column[1] == '0':
        return f"Invalid score column data for {extrande}"

    def delete_lines_with_keyword(matrix, keyword="Other sequence"):
        return [row for row in matrix if all(keyword not in element for element in row)]
    
    cleaned_sample_matrix = delete_lines_with_keyword(sample_file)
    
    def select_loci(sample_file, loci): 
        loci_rows = [row for row in sample_file if any(loci in element for element in row)]
        return [(row[1], row[2]) for row in loci_rows[:2]]  # Extract first two candidates

    def PercentNumber(sample_name, loci):
        try:
            # Map locus name to index in loci_list
            get_locus_index = lambda locus_name: loci_list.index(locus_name) + 1

            # Select loci data
            loci_data = select_loci(cleaned_sample_matrix, loci)

            if len(loci_data) == 0:
                return {loci: "Insufficient loci data"}

            # Get the score column value for the locus
            locus_index = get_locus_index(loci)
            if locus_index >= len(score_column) or not score_column[locus_index].isdigit():
                return {loci: "Invalid score column data"}

            locus_score = int(score_column[locus_index])
            if locus_score == 0:
                locus_score = 1

            # Process first candidate
            candidate1 = loci_data[0]
            first_percent = round(int(candidate1[1]) / locus_score, 3)

            if len(loci_data) == 1:
                return {loci: {"first_candidate": (candidate1[0], len(candidate1[0]), locus_score, first_percent)}}

            # Process second candidate
            candidate2 = loci_data[1]
            second_percent = round(int(candidate2[1]) / locus_score, 3)

            return {
                loci: {
                    "first_candidate": (candidate1[0], len(candidate1[0]), locus_score, int(candidate1[1]), first_percent),
                    "second_candidate": (candidate2[0], len(candidate2[0]), locus_score, int(candidate2[1]) second_percent),
                }
            }
        except (ValueError, IndexError) as e:
            return {loci: f"Error processing loci data: {str(e)}"}
    
    results = {}
    for loci in loci_list:
        results.update(PercentNumber(extrande, loci))
    
    # return json.dumps(results, indent=4)
    return results



output_csv_filename = '../new_output/AlleleInfo.csv'

with open(output_csv_filename, mode='w', newline='') as file:
    writer = csv.writer(file)

    # Write header
    writer.writerow(["Serial", "Sample", "Loci", "Allele1", "Allele2", "Read_count", "Frequency1", "Frequency2"])

    serial = 1  # Serial number starts at 1

    for i in range(len(sample_list)):
        sample_name = sample_list[i]
        data = process_sample(sample_name)

        if isinstance(data, str):  # Handle error messages
            print(data)
            continue

        for loci, candidates in data.items():
            if isinstance(candidates, str):  # Handle insufficient loci data
                writer.writerow([serial, sample_name, loci, "N/A", "N/A", "N/A", "N/A", "N/A"])
            else:
                # Extract allele 1 details
                allele1_seq, allele1_length, Read_count, allele1_freq = candidates["first_candidate"]

                # Extract allele 2 details if present
                if "second_candidate" in candidates:
                    allele2_seq, allele2_length, Read_count, allele2_freq = candidates["second_candidate"]
                else:
                    allele2_seq, allele2_length, Read_count, allele2_freq = "N/A", "N/A", "N/A", "N/A"

                # Write row to CSV
                writer.writerow([serial, sample_name, loci, allele1_seq, allele2_seq, Read_count, allele1_freq, allele2_freq])

            serial += 1  # Increment serial number

    print(f"CSV file '{output_csv_filename}' created successfully!")