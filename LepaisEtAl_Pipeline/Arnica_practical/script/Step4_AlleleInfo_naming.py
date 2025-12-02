import csv
import sys  # Import sys to read command-line arguments
import os

# --- Parameter Part ---
# Usage examples:
# 1. Standard: python Step4_AlleleInfo_naming.py 2 human
# 2. Shortcut: python Step4_AlleleInfo_naming.py human (Defaults to N=2)
# 3. Default:  python Step4_AlleleInfo_naming.py (Defaults to N=2, Mode=genalex)

# Defaults
N_ALLELES = 2
MODE_INPUT = 'human'

# Map modes to Column Indices based on your dictionary
# Index 1: Human_ID, Index 2: GenAlEx_ID, Index 3: Original_AlleleSeqCode
MODE_MAP = {
    'human': 1,
    'genalex': 2,
    'original': 3
}

# Argument Parsing Logic
if len(sys.argv) > 1:
    first_arg = sys.argv[1]
    
    # Case A: First argument is a number (Ploidy/N_ALLELES)
    if first_arg.isdigit():
        N_ALLELES = int(first_arg)
        
        # Check for second argument (Mode)
        if len(sys.argv) > 2:
            MODE_INPUT = sys.argv[2].lower()
            
    # Case B: First argument is NOT a number (Likely Mode)
    else:
        potential_mode = first_arg.lower() # Handle case-insensitivity
        if potential_mode in MODE_MAP:
            MODE_INPUT = potential_mode
            N_ALLELES = 2  # Default assumption for shortcut
            print(f"Detected mode '{MODE_INPUT}' as first argument. Setting N_ALLELES to default (2).")
        else:
            print(f"⚠️ Warning: Argument '{first_arg}' is not a valid number or mode.")
            print("   Using defaults: N_ALLELES=2, MODE='genalex'")

# Final Validation
if MODE_INPUT not in MODE_MAP:
    print(f"⚠️ Warning: Unknown mode '{MODE_INPUT}'. Using default 'genalex'.")
    MODE_INPUT = 'genalex'

TARGET_INDEX = MODE_MAP[MODE_INPUT]

print(f"✅ Processing {N_ALLELES} allele columns.")
print(f"✅ Output Mode: {MODE_INPUT} (Column Index {TARGET_INDEX})")

# --- End Parameter Part ---


def load_dictionary(dict_path):
    dictionary = []
    try:
        with open(dict_path, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # skip header
            for row in reader:
                # We need at least 5 columns to reach Full_Sequence (Index 4)
                if len(row) >= 5:
                    dictionary.append(row)
    except FileNotFoundError:
        print(f"Error: Dictionary file not found at {dict_path}")
        sys.exit(1)
    return dictionary


def search_in_dictionary(dictionary, search_string, target_col_idx):
    # Optimization: Don't search for "N/A"
    if search_string in ["N/A", ""]:
        return "N/A"
        
    for row in dictionary:
        # Full_Sequence is at Index 4 (Column 5)
        # Target Name is at target_col_idx
        if row[4] == search_string:
            return row[target_col_idx]
            
    return search_string  # If not found, keep original sequence


alleleinfo_path = '../new_output/Step2_AlleleInfo.csv'

ORIGINAL_DICT_PATH = '../new_output/Step3_AlleleInfo_dictionary.csv'
PATCHED_DICT_PATH = '../new_output/Step3b_Patched_AlleleInfo_dictionary.csv'

# 2. Decide which dictionary to use
if os.path.exists(PATCHED_DICT_PATH):
    dict_path = PATCHED_DICT_PATH
    print("Step 4: use patched dictionary (Step 3b output).")
else:
    dict_path = ORIGINAL_DICT_PATH
    print("Step 4: use original dictionary (Step 3 output).")

# 3. Define output path (also update naming convention)
output_path = f'../new_output/Step4_AlleleInfo_named_{MODE_INPUT}.csv'

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

        # Define the start and end index for allele columns
        allele_start_index = 3  # Corresponds to Allele1 (row[3])
        allele_end_index = allele_start_index + N_ALLELES 

        for row in reader:
            if not row:
                continue

            # Loop from Allele1 up to AlleleN
            for i in range(allele_start_index, allele_end_index):
                if i < len(row):
                    # Pass the target_index to the search function
                    row[i] = search_in_dictionary(dictionary, row[i], TARGET_INDEX)
            
            writer.writerow(row)

    print(f"✅ Finished! Saved to: {output_path}")

except FileNotFoundError:
    print(f"Error: Input file not found at {alleleinfo_path}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")