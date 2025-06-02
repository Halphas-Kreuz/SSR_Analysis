from dictionary_1 import search_in_matrix
import csv

with open('../filtered_results/filtered_tssvResults/AM30-01_S540_A_tssv_filtered.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    header = next(reader)  # Skip the header row
    data = [row for row in reader]

# Process each row in the matrix
processed_data = []
for row in data:
    # Join the row into a single string, split by '\t', and remove the split elements
    split_row = "\t".join(row).split("\t")
    processed_row = [element for element in split_row if element.strip()]  # Remove empty elements
    processed_data.append(processed_row)

# Print the processed data
for row in processed_data:
    print(row)