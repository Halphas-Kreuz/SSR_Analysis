import csv

def load_data():
    # Read the first CSV file into a list
    with open('../new_output/sorted_frequency_with_first_occurrence.csv', 'r') as file1:
        reader1 = csv.reader(file1)
        table1 = [row for row in reader1]

    # Read the second CSV file into a list
    with open('../new_output/AlleleInformationFile_nSSR_FullLength_ParameterSet2_sa70_sb10_m10_n20.csv', 'r') as file2:
        reader2 = csv.reader(file2)
        table2 = [row for row in reader2]
        data_rows = table2[1:]  # Skip the header row

    matrix = []
    for row in data_rows:
        if len(row) >= 5:  # Ensure sufficient columns
            fused_index = row[0] + "_" + row[2]
            triplet = [fused_index, row[1], row[4], len(row[4])]  # 2nd, 3rd, and 5th columns
            matrix.append(triplet)
    return matrix

def search_in_matrix(matrix, search_string):
    for row in matrix:
        if len(row) >= 3:  # Ensure the row has at least 3 elements
            if row[1] == search_string:  # Check if the second element matches
                return row[2]  # Return the third element
            elif row[2] == search_string:  # Check if the third element matches
                return row[1]  # Return the second element
    return None  # Return None if no match is found

# Protect executable code
if __name__ == "__main__":
    matrix = load_data()  # Load data only when the script is run directly
    search_string = "TGTGTGTCTATATATC(1)CATA(13)CACATGTATATATAT(1)"  # Replace with the string you want to search for
    result = search_in_matrix(matrix, search_string)
    if result:
        print(f"Found match! The corresponding element is: {result}")
    else:
        print("No match found.")