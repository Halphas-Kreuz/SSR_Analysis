
import csv
import re
import sys # Added sys for error handling

alleleinfo_path = '../new_output/AlleleInfo_named.csv'
comparison_path = '../new_output/table/All_SequenceComparisons.csv'
output_path = '../new_output/AlleleInfo_named_fused.csv'

# Load All_SequenceComparisons.csv into a lookup dictionary
comparison_dict = {}
try:
    with open(comparison_path, newline='') as compfile:
        reader = csv.DictReader(compfile)
        for row in reader:
            key = (row['Allele1'], row['Allele2'])
            comparison_dict[key] = row['Value']
except FileNotFoundError:
    print(f"Error: Comparison file not found at {comparison_path}")
    sys.exit(1) # Use sys.exit

# A single, unified set of all invalid tokens
INVALID_SEQUENCES = {'no data', 'other sequence', 'other sequences', 'n/a'}

def get_comparison_value(valA, valB):
    """
    Applies the comparison and fallback logic from the original script
    to any two allele values (valA, valB).
    """
    value_to_add = None # Use None to track if a value has been set
    
    # Check validity of both alleles
    is_valA_valid = valA and valA.lower() not in INVALID_SEQUENCES
    is_valB_valid = valB and valB.lower() not in INVALID_SEQUENCES
    
    # --- LOGIC from original file START ---

    # 1. Prioritize direct comparison only if BOTH columns have valid sequences
    if is_valA_valid and is_valB_valid:
        key = (valA, valB)
        value_to_add = comparison_dict.get(key)
    
    # 2. If no comparison was found, apply the fallback rule
    if value_to_add is None:
        # Check if valB contains an invalid sequence
        is_valB_invalid = not valB or valB.lower() in INVALID_SEQUENCES

        # Rule: Assign '1' ONLY if valA is valid AND valB is invalid -> no more needed
        if is_valA_valid and is_valB_invalid:
            value_to_add = '0'
        # For ALL other fallback cases, assign '0'
        else:
            value_to_add = '0'

    # --- LOGIC from original file END ---
    
    return str(value_to_add)

# --- Main Processing ---
try:
    with open(alleleinfo_path, newline='') as infile, open(output_path, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        # Read the header
        header = next(reader)
        
        # --- (MODIFICATION: All-Pairs) ---
        # Find all "Allele" columns and store their (index, name)
        allele_cols = [] # Will store tuples of (index, header_name)
        allele_pattern = re.compile(r'^Allele(\d+)$') # Captures the number
        for i, col_name in enumerate(header):
            match = allele_pattern.match(col_name)
            if match:
                allele_cols.append((i, col_name))
        
        if not allele_cols:
            print("Error: No 'Allele#' columns found in the header.")
            sys.exit(1)
            
        print(f"Found {len(allele_cols)} Allele columns. Generating all-pairs comparisons...")

        # Create new header columns for ALL pairs (e.g., Comp1_2, Comp1_3, Comp2_1, ...)
        new_comp_cols = []
        # Use nested loops to get all pairs
        for idx_i, name_i in allele_cols:
            num_i = allele_pattern.match(name_i).group(1) # Get "1" from "Allele1"
            
            for idx_j, name_j in allele_cols:
                if idx_i == idx_j: # Don't compare an allele to itself
                    continue
                
                num_j = allele_pattern.match(name_j).group(1) # Get "2" from "Allele2"
                new_comp_cols.append(f"Comp{num_i}_{num_j}")
            
        writer.writerow(header + new_comp_cols)
        # --- (END MODIFICATION) ---
        
        # Process each row
        for row in reader:
            if not row: # Skip empty rows
                continue

            comparison_values = []
            
            # --- (MODIFICATION: All-Pairs) ---
            # Loop through all pairs of alleles again to get values
            for idx_i, _ in allele_cols:
                for idx_j, _ in allele_cols:
                    if idx_i == idx_j: # Skip self-comparison
                        continue
                    
                    # Safely get values from columns, stripping whitespace
                    val_i = row[idx_i].strip() if idx_i < len(row) else ''
                    val_j = row[idx_j].strip() if idx_j < len(row) else ''
                    
                    # Get the comparison value using the refactored function
                    comp_val = get_comparison_value(val_i, val_j)
                    comparison_values.append(comp_val)
            # --- (END MODIFICATION) ---
            
            # Write the original row plus all the new comparison values
            writer.writerow(row + comparison_values)

    print(f'Fused file written to {output_path}')

except FileNotFoundError:
    print(f"Error: Input file not found at {alleleinfo_path}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")