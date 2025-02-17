import csv
import os


#step 1 : find out the sample to be searched (if found, then read it )

test_sample_name = 'AM30-01_S540_A'
loci_list = ["Arm01", "Arm03", "Arm04", "Arm06", "Arm07", "Arm08"]  

def search_folder(folder_path, keyword):
    for file_name in os.listdir(folder_path):
        if keyword in file_name and file_name.endswith('.csv'):
            file_path = os.path.join(folder_path, file_name)
            with open(file_path, 'r') as file:
                reader = csv.reader(file, delimiter='\t')  # Use tab as the delimiter
                matrix = [row for row in reader] 
            # print(f"Contents of {file_name} as matrix:")
            # for row in matrix:
            #     print(row)
            # print("-" * 50)
            return matrix 
        
sample_file= search_folder('../filtered_part/filtered_2_line_stuttermarked', test_sample_name)


# step 2 : construct the score matrix

#only need to initial once : there is only a score matrix!
score_matrix = []
def open_file(file_path):
    with open( file= file_path, newline='',) as file:
        reader = csv.reader(file,delimiter = '\t')
        for row in reader: 
            score_matrix.append(row)

#step 3 : with the keyword, search for the right column 
def search_column(matrix, word):
    if word in matrix[0]:
        col_index = matrix[0].index(word)
        column = [row[col_index] for row in matrix]
        print(f"Column '{word}' found at index {col_index}")
        print(column)
        return column
    else:
        print(f"Column '{word}' not found")
        return 0

open_file('../LocusCoverageperIndividual_nSSR_FullLength.csv')
search_column(score_matrix,test_sample_name) #first print : from the score matrix , give every loci count 

#step 4 : from the filtered_2_line document, reclean the data 

#4-1 : from the sample file, clean the line with "Other sequence"
def delete_lines_with_keyword(matrix, keyword="Other sequence"):
    filtered_matrix = [row for row in matrix if all(keyword not in element for element in row)]
    return filtered_matrix

cleaned_empty_sample_matrix = delete_lines_with_keyword(sample_file)

#4-2 : from the sample file (also the stuttermark), extract the count for the loci 
def select_loci(sample_file, loci='Arm01'): 
    candidate1 = None
    candidate2 = None

    loci_rows = [row for row in sample_file if any(loci in element for element in row)]  # Find matching rows

    for i, row in enumerate(loci_rows):
        # print(row)  # Print matching rows

        if i == 0:
            candidate1 = row[2]
        elif i == 1:
            candidate2 = row[2]
    # print(candidate1, candidate2)
    return (candidate1, candidate2)  # Return the selected candidates (two values )

select_loci(cleaned_empty_sample_matrix)

def PercentNumber(sample_name, loci="Arm01"):
    score_column = search_column(score_matrix, sample_name)  
    loci_count = select_loci(cleaned_empty_sample_matrix, loci)  # Dynamically select loci

    if score_column and len(score_column) > 1 and score_column[1] != 0:  # Avoid division by zero
        first_candidate = round(int(loci_count[0]) / int(score_column[1]), 3) if loci_count[0] is not None else None
        second_candidate = round(int(loci_count[1]) / int(score_column[1]), 3) if loci_count[1] is not None else None

        print(f"First candidate for {sample_name} {loci} is {first_candidate}")
        print(f"Second candidate for {sample_name} {loci} is {second_candidate}")
    else:
        print(f"Invalid score column data for {sample_name}")

PercentNumber(test_sample_name,"Arm03")



