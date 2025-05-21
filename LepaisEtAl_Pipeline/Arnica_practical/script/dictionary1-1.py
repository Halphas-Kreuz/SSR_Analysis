#this script read the sorted_frequency_wuith_first_occurrence.csv file and Alleles.csv file 
#Main focus : both of the file contain the main seq
#thus : a translation table is created

# Read the first CSV file into a list
import csv

with open('../new_output/sorted_frequency_with_first_occurrence.csv', 'r') as file1:
    reader1 = csv.reader(file1)
    table1 = [row for row in reader1]

# Read the second CSV file into a list
with open('../new_output/AlleleInformationFile_nSSR_FullLength_ParameterSet2_sa70_sb10_m10_n20.csv', 'r') as file2:
    reader2 = csv.reader(file2)
    table2 = [row for row in reader2]

# Extract the 2nd and 5th elements from each row and create a dictionary
# the dictionary need to be in both direction, let me think how to do this .
# Assuming `table2` is already populated with the CSV data
abbr_to_fullSeqDict = {row[1]: row[4] for row in table2 if len(row) > 4}
reversed_abbr_to_fullSeqDict = {row[4]: row[1] for row in table2 if len(row) > 4}

def get_counterpart(forward, reverse, query):
    if query in forward:
        return forward[query]
    elif query in reverse:
        return reverse[query]
    else:
        return None

#try for a key 
word_to_search = "TGTGTGTCTATATATC(1)CATA(13)CACATGTATATATAT(1)"  # Replace with the word you want to search for
reverseSearch = "CATGGCTGCCACACCCACTTCTTAAACCATCATCATTTTCATGTTGTATGTATGTATGTATGATTGATGATGCATGAGATA"

result1= get_counterpart(abbr_to_fullSeqDict, reversed_abbr_to_fullSeqDict, word_to_search)
result2= get_counterpart(abbr_to_fullSeqDict, reversed_abbr_to_fullSeqDict, reverseSearch)
print("Result for word_to_search:", result1)
print("Result for reverseSearch:", result2)