import csv

matrix = []

with open('../LocusCoverageperIndividual_nSSR_FullLength.csv', newline='',) as file:
    reader = csv.reader(file,delimiter = '\t')
    for row in reader: 
        matrix.append(row)

def search_column(matrix, word):
    if word in matrix[0]:
        col_index = matrix[0].index(word)
        


