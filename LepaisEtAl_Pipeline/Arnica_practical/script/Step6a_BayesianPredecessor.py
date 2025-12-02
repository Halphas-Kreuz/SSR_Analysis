import csv
import re
import sys
import os

# --- Parameter Part ---
# Usage: python Step6a_BayesianPredecessor.py [MODE]
# Modes: 'human' (default), 'genalex', 'original'

# Defaults
MODE_INPUT = 'human'

# Argument Parsing
if len(sys.argv) > 1:
    potential_mode = sys.argv[1].lower()
    # Basic check to avoid filenames being treated as modes blindly, 
    # though here we trust the user input mainly.
    MODE_INPUT = potential_mode

print(f"✅ Processing Bayesian Predecessor Step.")
print(f"✅ Mode: {MODE_INPUT}")

# --- File Paths ---
# INPUT 1: The Named Allele File from Step 4
# We look for the file corresponding to the requested mode
alleleinfo_path = f'../new_output/Step4_AlleleInfo_named_{MODE_INPUT}.csv'

# INPUT 2: The Comparison Table (Confirmed stable)
comparison_path = '../new_output/table/Step6_All_SequenceComparisons.csv'

# OUTPUT: The Fused File
output_path = f'../new_output/Step6a_AlleleInfo_named_fused_{MODE_INPUT}.csv'

# --- Main Script Logic ---

# 1. Load All_SequenceComparisons.csv into a lookup dictionary
comparison_dict = {}
try:
    with open(comparison_path, newline='') as compfile:
        reader = csv.DictReader(compfile)
        # Verify required columns exist
        if not {'Allele1', 'Allele2', 'Value'}.issubset(reader.fieldnames):
             print(f"❌ Error: {comparison_path} is missing required columns (Allele1, Allele2, Value).")
             sys.exit(1)

        for row in reader:
            # Create a tuple key: (Allele1, Allele2)
            key = (row['Allele1'], row['Allele2'])
            comparison_dict[key] = row['Value']
    print(f"✅ Loaded {len(comparison_dict)} comparisons from reference table.")

except FileNotFoundError:
    print(f"❌ Error: Comparison file not found at {comparison_path}")
    print("👉 Solution: Please ensure the 'table' folder and comparison file exist.")
    sys.exit(1)

# A single, unified set of all invalid tokens
INVALID_SEQUENCES = {'no data', 'other sequence', 'other sequences', 'n/a', 'na', ''}

def get_comparison_value(valA, valB):
    """
    Look up the distance/value between valA and valB.
    """
    value_to_add = None 
    
    # Check validity of both alleles
    is_valA_valid = valA and valA.lower() not in INVALID_SEQUENCES
    is_valB_valid = valB and valB.lower() not in INVALID_SEQUENCES
    
    # 1. Prioritize direct comparison only if BOTH columns have valid sequences
    if is_valA_valid and is_valB_valid:
        key = (valA, valB)
        value_to_add = comparison_dict.get(key)
        
        # If not found directly, try reverse order (if the table isn't symmetric)
        if value_to_add is None:
             key_reverse = (valB, valA)
             value_to_add = comparison_dict.get(key_reverse)
    
    # 2. If no comparison was found (or one was invalid), apply fallback
    if value_to_add is None:
        # If valA is valid but valB is invalid -> '0'
        # If both invalid -> '0'
        # Basically, default to '0' if lookup fails
        value_to_add = '0'

    return str(value_to_add)

# --- Processing Input File ---
try:
    with open(alleleinfo_path, newline='') as infile, open(output_path, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        # Read the header
        try:
            header = next(reader)
        except StopIteration:
            print(f"❌ Error: Input file {alleleinfo_path} is empty.")
            sys.exit(1)
        
        # Find all "Allele" columns (Allele1, Allele2...)
        allele_cols = [] 
        allele_pattern = re.compile(r'^Allele(\d+)$') 
        
        for i, col_name in enumerate(header):
            match = allele_pattern.match(col_name)
            if match:
                allele_cols.append((i, col_name))
        
        if not allele_cols:
            print("❌ Error: No 'Allele#' columns found in the header.")
            sys.exit(1)
            
        print(f"Found {len(allele_cols)} Allele columns. Generating all-pairs comparisons...")

        # Create new header columns for ALL pairs
        new_comp_cols = []
        for idx_i, name_i in allele_cols:
            num_i = allele_pattern.match(name_i).group(1)
            for idx_j, name_j in allele_cols:
                if idx_i == idx_j: 
                    continue
                num_j = allele_pattern.match(name_j).group(1)
                new_comp_cols.append(f"Comp{num_i}_{num_j}")
            
        writer.writerow(header + new_comp_cols)
        
        # Process each row
        row_count = 0
        for row in reader:
            if not row: continue
            
            comparison_values = []
            
            # Loop through all pairs
            for idx_i, _ in allele_cols:
                for idx_j, _ in allele_cols:
                    if idx_i == idx_j: continue
                    
                    val_i = row[idx_i].strip() if idx_i < len(row) else ''
                    val_j = row[idx_j].strip() if idx_j < len(row) else ''
                    
                    comp_val = get_comparison_value(val_i, val_j)
                    comparison_values.append(comp_val)
            
            writer.writerow(row + comparison_values)
            row_count += 1

    print(f"✅ Finished! Processed {row_count} rows.")
    print(f"✅ Fused file written to: {output_path}")

except FileNotFoundError:
    print(f"❌ Error: Input file not found at {alleleinfo_path}")
    print(f"👉 Solution: Please run Step 4 first with mode '{MODE_INPUT}'.")
    sys.exit(1)
except Exception as e:
    print(f"An unexpected error occurred: {e}")