import csv

alleleinfo_path = '../new_output/AlleleInfo_named.csv'
comparison_path = '../new_output/table/All_SequenceComparisons.csv'
output_path = '../new_output/AlleleInfo_named_fused.csv'

# Load All_SequenceComparisons.csv into a lookup dictionary
comparison_dict = {}
with open(comparison_path, newline='') as compfile:
    reader = csv.DictReader(compfile)
    for row in reader:
        key = (row['Allele1'], row['Allele2'])
        comparison_dict[key] = row['Value']

# A single, unified set of all invalid tokens to be applied to all columns.
INVALID_SEQUENCES = {'no data', 'other sequence', 'other sequences', 'n/a'}

# Process AlleleInfo_named.csv and fuse values
with open(alleleinfo_path, newline='') as infile, open(output_path, 'w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)
    
    header = next(reader)
    writer.writerow(header + ['ComparisonValue'])
    
    for row in reader:
        # Safely get values from columns, stripping whitespace.
        val4 = row[3].strip() if len(row) > 3 else ''
        val5 = row[4].strip() if len(row) > 4 else ''
        
        value_to_add = None # Use None to track if a value has been set.
        
        # --- FINAL LOGIC START ---

        # 1. Prioritize direct comparison only if BOTH columns have valid sequences.
        if val4 and val5 and val4.lower() not in INVALID_SEQUENCES and val5.lower() not in INVALID_SEQUENCES:
            key = (val4, val5)
            value_to_add = comparison_dict.get(key)
        
        # 2. If no comparison was found, apply the fallback rule.
        if value_to_add is None:
            # Check if column 4 contains a valid sequence.
            is_val4_valid = val4 and val4.lower() not in INVALID_SEQUENCES
            # Check if column 5 contains an invalid sequence.
            is_val5_invalid = not val5 or val5.lower() in INVALID_SEQUENCES

            # Rule: Assign '1' ONLY if column 4 is valid AND column 5 is invalid.
            if is_val4_valid and is_val5_invalid:
                value_to_add = '1'
            # For ALL other fallback cases, assign '0'.
            else:
                value_to_add = '0'

        # --- FINAL LOGIC END ---
        
        # Write the original row plus the newly determined value.
        writer.writerow(row + [str(value_to_add)])

print(f'Fused file written to {output_path}')