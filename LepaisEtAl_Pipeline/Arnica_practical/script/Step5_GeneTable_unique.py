import os  # ADDED: For checking file existence
import csv
import sys

# --- Parameter Part ---
# Usage: python Step5_GeneTable_unique.py [MODE]
# Modes: 'human' (default), 'genalex', 'original'

# Defaults
MODE_INPUT = 'human' 

# Map modes to Dictionary Column Indices (From Step 3 Dictionary)
# Index 1: Human_ID (e.g. Arm01_83) 
# Index 2: GenAlEx_ID (e.g. 83)
# Index 3: Original_AlleleSeqCode (e.g. 101)
MODE_MAP = {
    'human': 1,
    'genalex': 2,
    'original': 3
}

# Argument Parsing
if len(sys.argv) > 1:
    potential_mode = sys.argv[1].lower()
    if potential_mode in MODE_MAP:
        MODE_INPUT = potential_mode
    else:
        print(f"⚠️ Warning: Unknown mode '{sys.argv[1]}'. Using default 'human'.")

TARGET_INDEX = MODE_MAP[MODE_INPUT]

print(f"✅ Processing Genotype Table.")
print(f"✅ Output Mode: {MODE_INPUT} (Using Dictionary Column {TARGET_INDEX})")
# --- End Parameter Part ---


def load_lookup_dictionary(dict_path, target_idx):
    """
    Loads the dictionary CSV from Step 3b output.
    Key: Fused_Index (Column 0, e.g. Arm01_101)
    Value: The Target ID defined by MODE (Human/GenAlEx/Original)
    """
    lookup_dict = {}
    try:
        with open(dict_path, 'r', newline='') as f:
            reader = csv.reader(f)
            next(reader)  # Skip the header row
            for row in reader:
                # Ensure row is long enough for the target index
                if len(row) > target_idx:
                    fused_index = row[0]       # Key: Arm01_101
                    target_id = row[target_idx] # Value: Target Name
                    lookup_dict[fused_index] = target_id
    except FileNotFoundError:
        # NOTE: File check is now handled by os.path.exists outside this function, 
        # so this block is for internal function safety only.
        print(f"Error: Dictionary file not found at '{dict_path}'")
        sys.exit(1) 
    return lookup_dict

# --- File Paths and MANDATORY Dependency Check ---
# INPUT 1: The Raw Genotypic Table 
genotype_input_file = "../GenotypicTable_nSSR_FullLength_ParameterSet1_sa50_sb10_m15_n20.txt"

# INPUT 2: The Patched Dictionary from Step 3b (Mandatory)
PATCHED_DICT_PATH = '../new_output/Step3b_Patched_AlleleInfo_dictionary.csv'

# Check for mandatory dependency (Step 3b output)
if os.path.exists(PATCHED_DICT_PATH):
    dictionary_file = PATCHED_DICT_PATH
    print("Step 5: Using MANDATORY patched dictionary.")
else:
    # If the patched file is missing, fail immediately and clearly.
    print(f"[FAILED] Step 5 Error: Cannot run.")
    print(f"Solution: Please ensure Step 3b (Dictionary Patcher) ran successfully and created: {PATCHED_DICT_PATH}")
    sys.exit(1) # Stop script

# OUTPUT: New fused table with step number and mode in filename
# RENAMING: Step5_GenotypicTable_fused_[MODE].txt
output_file = f"../new_output/Step5_GenotypicTable_fused_{MODE_INPUT}.txt"

# --- Main Script Logic ---

# 1. Load the lookup dictionary
allele_lookup = load_lookup_dictionary(dictionary_file, TARGET_INDEX)

# 2. Process the main genotypic table
try:
    with open(genotype_input_file, "r") as fin, open(output_file, "w", newline='') as fout:
        lines = fin.readlines()
        
        if not lines:
             print(f"❌ Error: Input file '{genotype_input_file}' is empty.")
             sys.exit(1)

        headers = lines[0].strip().split("\t")

        writer = csv.writer(fout, delimiter='\t')
        writer.writerow(headers)

        for line in lines[1:]:
            cells = line.strip().split("\t")
            # Pad the row with 'NA' if it's shorter than the header
            if len(cells) < len(headers):
                cells += ["NA"] * (len(headers) - len(cells))

            new_cells = [cells[0]]  # Keep the first column (Sample Name)

            # Iterate through the remaining columns
            for i in range(1, len(headers)):
                # Determine the correct header for the allele (e.g., 'Arm01')
                pair_header_idx = i if i % 2 == 1 else i - 1
                pair_header = headers[pair_header_idx]
                value = cells[i]

                if value not in ["NA", ""]:
                    # Create the fused key to search (e.g., 'Arm01_109')
                    lookup_key = f"{pair_header}_{value}"

                    # Look up the key using our Mode-specific dictionary
                    final_allele_id = allele_lookup.get(lookup_key, "NA")
                    new_cells.append(final_allele_id)
                else:
                    new_cells.append("NA")

            writer.writerow(new_cells)

    print(f"✅ Transformation complete. Final table written to: {output_file}")

except FileNotFoundError:
    print(f"Error: Raw Genotype Input file not found at '{genotype_input_file}'")
    print("Solution: This file is usually the output of the initial SSRseq extraction. Please check the file path.")
    sys.exit(1)
except Exception as e:
    print(f"An unexpected error occurred: {e}")