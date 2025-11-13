import csv
import sys  # Import sys to read command-line arguments

# --- Parameter Part ---
# Get 'n' from command line, default to 2
try:
    # sys.argv[0] is the script name, sys.argv[1] is the first arg
    N_ALLELES = int(sys.argv[1])
except (IndexError, ValueError):
    N_ALLELES = 2  # Default to 2 if no argument is given

print(f"✅ Processing {N_ALLELES} allele columns for naming.")
# --- End Parameter Part ---


def load_dictionary(dict_path):
    dictionary = []
    try:
        with open(dict_path, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # skip header
            for row in reader:
                if len(row) >= 4:  # Ensure row has at least 4 columns
                    dictionary.append(row)
    except FileNotFoundError:
        print(f"Error: Dictionary file not found at {dict_path}")
        sys.exit(1) # Exit the script if dictionary is missing
    return dictionary


def search_in_dictionary(dictionary, search_string):
    # Optimization: Don't search for "N/A"
    if search_string == "N/A":
        return "N/A"
        
    for row in dictionary:
        # Assumes sequence is in 4th column (index 3) and name is in 2nd (index 1)
        if row[3] == search_string:
            return row[1]  # Return the name
    return search_string  # If not found, keep original sequence


alleleinfo_path = '../new_output/AlleleInfo.csv'
dict_path = '../new_output/AlleleInfo_dictionary.csv'
output_path = '../new_output/AlleleInfo_named.csv'

# Load dictionary
dictionary = load_dictionary(dict_path)
if not dictionary:
    print("Dictionary is empty or could not be loaded. Exiting.")
    sys.exit(1)

# Read, convert, and write to a new file
try:
    with open(alleleinfo_path, 'r') as infile, open(output_path, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        # Read and write header
        try:
            header = next(reader)
        except StopIteration:
            print(f"Error: Input file {alleleinfo_path} is empty.")
            sys.exit(1)
            
        writer.writerow(header)

        # --- (MODIFICATION) ---
        # Define the start and end index for allele columns
        allele_start_index = 3  # Corresponds to Allele1 (row[3])
        # The loop will go up to, but not include, this index
        allele_end_index = allele_start_index + N_ALLELES 
        # --- (END MODIFICATION) ---

        for row in reader:
            # Ensure row is not empty
            if not row:
                continue

            # --- (MODIFICATION) ---
            # Loop from Allele1 (index 3) up to AlleleN
            for i in range(allele_start_index, allele_end_index):
                # Check if the row has this many columns (avoids errors on short rows)
                if i < len(row):
                    # Pass the value (e.g., row[3]) to the search function
                    row[i] = search_in_dictionary(dictionary, row[i])
            # --- (END MODIFICATION) ---
            
            writer.writerow(row)

    print(f"AlleleInfo with names saved to {output_path}")

except FileNotFoundError:
    print(f"Error: Input file not found at {alleleinfo_path}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")