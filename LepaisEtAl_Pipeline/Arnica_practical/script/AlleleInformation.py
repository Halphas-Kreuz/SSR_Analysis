import csv
import os

Allele = '../AlleleInformationFile_nSSR_FullLength_ParameterSet2_sa70_sb10_m10_n20.csv'

def read_csv_file(file_path):
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')  # Use tab as the delimiter
        matrix = [row for row in reader]  # Convert the reader object to a list of lists
    return matrix

matrix = read_csv_file(Allele)
print(matrix[1])

#step 1 : check the locus 

