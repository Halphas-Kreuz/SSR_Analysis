import csv
import os


#step 1 : find out the sample to be searched (if found, then read it )

sample_name = 'AM30-01_S540_A'

def search_folder(folder_path, keyword):
    for file_name in os.listdir(folder_path):
        if keyword in file_name and file_name.endswith('.csv'):
            file_path = os.path.join(folder_path, file_name)
            with open(file_path, 'r') as file:
                reader = csv.reader(file, delimiter='\t')  # Use tab as the delimiter
                matrix = [row for row in reader] 
            print(f"Contents of {file_name} as matrix:")
            for row in matrix:
                print(row)
            print("-" * 50)
            return matrix 
        
sample_file= search_folder('../filtered_part/filtered_2_line_stuttermarked', sample_name)


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
search_column(score_matrix,sample_name)

#step 4 : from the filtered_2_line document, reclean the data 

#4-1 : clean the line with "Other sequence"
def delete_lines_with_keyword(matrix, keyword="Other sequence"):

    # filtered_content = [line for line in content if keyword not in line]
    # return filtered_content
    filtered_matrix = [row for row in matrix if all(keyword not in element for element in row)]
    return filtered_matrix

# def select_loci(sample_matrix, score_matrix): #try for Arm01
    
    
# print(delete_lines_with_keyword(sample_file))

    




