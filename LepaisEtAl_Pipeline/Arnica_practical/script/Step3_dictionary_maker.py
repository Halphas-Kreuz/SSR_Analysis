import csv
import os

def txt_to_csv(txt_path, csv_path):
    if not os.path.exists(csv_path):
        with open(txt_path, 'r') as infile, open(csv_path, 'w', newline='') as outfile:
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile)
            for row in reader:
                writer.writerow(row)
        print(f"Converted {txt_path} to {csv_path}")

def load_data(csv_path):
    with open(csv_path, 'r') as file2:
        reader2 = csv.reader(file2)
        table2 = [row for row in reader2]
        data_rows = table2[1:]  # Skip the header row

    matrix = []
    for row in data_rows:
        if len(row) >= 5:  # Ensure sufficient columns
            fused_index = row[0] + "_" + row[2]
            triplet = [fused_index, row[1], row[4], len(row[4]), row[0]]  # Append original row[0]
            matrix.append(triplet)
    return matrix

def search_in_matrix(matrix, search_string):
    for row in matrix:
        if len(row) >= 3:
            if row[1] == search_string:
                return row[2]
            elif row[2] == search_string:
                return row[1]
    return None

if __name__ == "__main__":
    txt_path = '../AlleleInformationFile_nSSR_FullLength_ParameterSet1_sa50_sb10_m15_n20.txt'
    csv_path = '../AlleleInformationFile_nSSR_FullLength_ParameterSet1_sa50_sb10_m15_n20.csv'
    output_path = '../new_output/AlleleInfo_dictionary.csv'

    # Convert TXT to CSV if needed
    txt_to_csv(txt_path, csv_path)

    # Load data from CSV
    matrix = load_data(csv_path)

    # Example search
    search_string = "TGTGTGTCTATATATC(1)CATA(13)CACATGTATATATAT(1)"  # Replace as needed
    result = search_in_matrix(matrix, search_string)
    if result:
        print(f"Found match! The corresponding element is: {result}")
    else:
        print("No match found.")

    # Write the columns to a new CSV file, including the original row[0] at the end
    with open(output_path, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['Fused_Index', 'Column2', 'Column5', 'Len_Column5', 'Original_Column0'])  # Header
        for row in matrix:
            writer.writerow(row)

    print(f"Triplet CSV file saved to {output_path}")