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
        # print(f"Column '{word}' found at index {col_index}")
        # print(column)
        return column
    else:
        # print(f"Column '{word}' not found")
        return 0

open_file('../LocusCoverageperIndividual_nSSR_FullLength.csv')
search_column(score_matrix,test_sample_name) #first print : from the score matrix , give every loci count 
#example output : ['AM30-01_S540_A', '2109', '1764', '3352', '1439', '264', '3260']

#step 4 : from the filtered_2_line document, reclean the data 

#4-1 : from the sample file, clean the line with "Other sequence"
def delete_lines_with_keyword(matrix, keyword="Other sequence"):
    filtered_matrix = [row for row in matrix if all(keyword not in element for element in row)]
    return filtered_matrix

cleaned_empty_sample_matrix = delete_lines_with_keyword(sample_file)
#Example output :(including all loci)
# Arm01	TGTGTGTCTATATATC(1)CATA(14)CACATGTATATATAT(1)	1782	1782	0	
# Arm01	TGTGTGTCTATATATC(1)CATA(13)CACATGTATATATAT(1)	159	159	0	STUTTER:1247.4x1(2-1)

#4-2 : from the sample file (also the stuttermark), extract the name and count for the loci 
def select_loci(sample_file, loci='Arm01'): 
    candidate1_Name = None
    candidate2_Name = None
    candidate1_Count = None
    candidate2_Count = None

    loci_rows = [row for row in sample_file if any(loci in element for element in row)]  # Find matching rows

    for i, row in enumerate(loci_rows):
        # print(row)  # Print matching rows

        if i == 0:
            candidate1_Name = row[1]
            candidate1_Count = row[2]
        elif i == 1:
            candidate2_Name = row[1]
            candidate2_Count = row[2]
    # print(candidate1, candidate2)
    print (candidate1_Name,candidate1_Count,candidate2_Name,candidate2_Count) 
    return (candidate1_Name,candidate1_Count,candidate2_Name,candidate2_Count)  # Return the selected candidates (two values )

select_loci(cleaned_empty_sample_matrix)

def PercentNumber(sample_name, loci="Arm01"):
    score_column = search_column(score_matrix, sample_name)  #load the sum of loci count 
    loci_line = select_loci(cleaned_empty_sample_matrix, loci)  #load the loci name and count (beware, this tuple has four values )

    if score_column and len(score_column) > 1 and score_column[1] != 0:  # Avoid division by zero
        first_candidate_percent = round(int(loci_line[1]) / int(score_column[1]), 3) if loci_line[1] is not None else None
        second_candidate_percent = round(int(loci_line[3]) / int(score_column[1]), 3) if loci_line[3] is not None else None

        print(f"First candidate for {sample_name} {loci} is {loci_line[0]} {first_candidate_percent}")
        print(f"Second candidate for {sample_name} {loci} is {loci_line[2]} {second_candidate_percent}")
    else:
        print(f"Invalid score column data for {sample_name}")

# PercentNumber(test_sample_name,"Arm03")

for i in loci_list:
    PercentNumber(test_sample_name,i)



