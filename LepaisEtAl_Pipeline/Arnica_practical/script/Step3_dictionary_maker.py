import csv
import os

def txt_to_csv(txt_path, csv_path):
    """Converts a tab-delimited text file to a CSV file if it doesn't already exist."""
    if not os.path.exists(csv_path):
        with open(txt_path, 'r') as infile, open(csv_path, 'w', newline='') as outfile:
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile)
            for row in reader:
                writer.writerow(row)
        print(f"Converted {txt_path} to {csv_path}")

def load_data(csv_path):
    """Loads data from CSV, transforms it by adding a new fused ID column, and ensures its uniqueness."""
    with open(csv_path, 'r') as file2:
        reader2 = csv.reader(file2)
        table2 = [row for row in reader2]
        data_rows = table2[1:]  # Skip the header row

    matrix = []
    # Dictionary to track the count of each ID to handle duplicates
    id_counts = {}

    for row in data_rows:
        if len(row) >= 5:  # Ensure the original row has sufficient columns
            # Original fused index from column 1 and 3 of the input file
            fused_index = f"{row[0]}_{row[2]}"

            # --- NEW FEATURE LOGIC ---
            # Calculate the length of the 5th column (index 4)
            len_col5 = len(row[4])
            # Create the base for the new fused ID from the 1st column (index 0) and the calculated length
            base_new_fused_id = f"{row[0]}_{len_col5}"

            # Check for duplicates and create a unique ID
            if base_new_fused_id in id_counts:
                # If the ID already exists, increment its count
                id_counts[base_new_fused_id] += 1
                # Append the count to make the ID unique (e.g., Arm03_131_2)
                unique_new_fused_id = f"{base_new_fused_id}_{id_counts[base_new_fused_id]}"
            else:
                # If it's the first time seeing this ID, initialize its count to 1
                id_counts[base_new_fused_id] = 1
                # The first instance of the ID does not get a suffix
                unique_new_fused_id = base_new_fused_id
            # --- END NEW FEATURE LOGIC ---

            # Assemble the new row, inserting the unique_new_fused_id as the second element
            output_row = [
                fused_index,         # 1st column: e.g., Arm03_103
                unique_new_fused_id, # 2nd column (NEW & UNIQUE): e.g., Arm03_131 or Arm03_131_2
                row[1],              # 3rd column: The long TGTGT... string
                row[4],              # 4th column: The other long TGTGT... string
                len_col5,            # 5th column: The length, e.g., 131
                row[0]               # 6th column: The original first part, e.g., Arm03
            ]
            matrix.append(output_row)
    return matrix

def search_in_matrix(matrix, search_string):
    """
    Searches for a string. NOTE: Indices are shifted due to the new column.
    The search now compares against the 3rd and 4th columns of the output.
    """
    for row in matrix:
        if len(row) >= 4:
            # Check against the new column indices for the searchable strings
            if row[2] == search_string:
                return row[3]
            elif row[3] == search_string:
                return row[1]
    return None

if __name__ == "__main__":
    txt_path = '../AlleleInformationFile_nSSR_FullLength_ParameterSet1_sa50_sb10_m15_n20.txt'
    csv_path = '../AlleleInformationFile_nSSR_FullLength_ParameterSet1_sa50_sb10_m15_n20.csv'
    output_path = '../new_output/AlleleInfo_dictionary.csv'

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Convert TXT to CSV if needed
    txt_to_csv(txt_path, csv_path)

    # Load data from CSV with the new transformation logic
    matrix = load_data(csv_path)

    # Example search
    # search_string = "TGTGTGTCTATATATC(1)CATA(13)CACATGTATATATAT(1)"  # Replace as needed
    # result = search_in_matrix(matrix, search_string)
    # if result:
    #     print(f"Found match! The corresponding element is: {result}")
    # else:
    #     print("No match found.")

    # Write the new structure to a CSV file with an updated header
    with open(output_path, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        # Update header to include the new column
        writer.writerow(['Fused_Index', 'Arm_Length_ID', 'Column2', 'Column5', 'Len_Column5', 'Original_Column0'])
        writer.writerows(matrix)

    print(f"Modified dictionary CSV file saved to {output_path}")

