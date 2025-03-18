import csv
import os
import json

folder_path = '../filtered_results/filtered_tssvResults'
sample_name = 'AM30-01_S540_A'
loci_list = ["Arm01", "Arm03", "Arm04", "Arm06", "Arm07", "Arm08", "Arm11"]  
locus_coverage_file = "../LocusCoverageperIndividual_nSSR_FullLength.csv"
AlleleInformation = "../AlleleInformationFile_nSSR_FullLength_ParameterSet2_sa70_sb10_m10_n20.csv"
def process_sample():
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

    score_column = search_column(score_matrix, sample_name)
    #if the count is < 0.25, then drop it 
    if not score_column or len(score_column) < 2 or score_column[1] == '0':
        return f"Invalid score column data for {sample_name}"

    def delete_lines_with_keyword(matrix, keyword="Other sequence"):
        return [row for row in matrix if all(keyword not in element for element in row)]
    
    cleaned_sample_matrix = delete_lines_with_keyword(sample_file)
    
    def select_loci(sample_file, loci): 
        loci_rows = [row for row in sample_file if any(loci in element for element in row)]
        return [(row[1], row[2]) for row in loci_rows[:2]]  # Extract first two candidates

    def PercentNumber(sample_name, loci):
        get_locus_index = lambda locus_name: loci_list.index(locus_name) + 1
        loci_data = select_loci(cleaned_sample_matrix, loci)
        
        if len(loci_data) == 0:
            return {loci: "Insufficient loci data"}
        
        candidate1 = loci_data[0]
        first_percent = round(int(candidate1[1]) / int(score_column[get_locus_index(loci)]), 3)
        
        if len(loci_data) == 1:
            return {loci: {"first_candidate": (candidate1[0], len(candidate1[0]), first_percent)}}
        
        candidate2 = loci_data[1]
        second_percent = round(int(candidate2[1]) / int(score_column[get_locus_index(loci)]), 3)
        
        return {loci: {"first_candidate": (candidate1[0], len(candidate1[0]),  first_percent), "second_candidate": (candidate2[0], len(candidate2[0]), second_percent)}}
    
    results = {}
    for loci in loci_list:
        results.update(PercentNumber(sample_name, loci))
    
    return json.dumps(results, indent=4)

result = process_sample()
print(result)