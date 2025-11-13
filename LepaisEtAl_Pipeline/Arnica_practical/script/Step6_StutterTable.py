import csv
import re
import sys 
from pathlib import Path
from collections import defaultdict

# --- Parameter Part ---
try:
    SCORE_MINUS_1 = float(sys.argv[1])
except (IndexError, ValueError):
    SCORE_MINUS_1 = 0.5  # Default value
try:
    SCORE_OTHER_DIFF = float(sys.argv[2])
except (IndexError, ValueError):
    SCORE_OTHER_DIFF = 0.1  # Default value

print(f"✅ Using comparison scores:")
print(f"  - Score for diff (-1): {SCORE_MINUS_1}")
print(f"  - Score for diff (1, 2, -2): {SCORE_OTHER_DIFF}")
# --- End Parameter Part ---


def parse_repeat_pattern(s):
    return re.findall(r'([A-Za-z]+)\((\d+)\)', s)


def compare_patterns_asym(bench, comp, score_val_minus_1, score_val_other):
    if bench == comp:
        return 1.0
    if len(bench) == len(comp):
        diffs = []
        diff_idx = -1
        for i, ((m1, n1), (m2, n2)) in enumerate(zip(bench, comp)):
            if m1 != m2:
                return 0.0  # Motifs differ
            if n1 != n2:
                diffs.append((i, int(n1) - int(n2)))
                diff_idx = i
        
        if len(diffs) == 1:
            diff = diffs[0][1]
            if diff == -1:
                return score_val_minus_1
            elif diff in [1, 2, -2]:
                return score_val_other
                
    return 0.0


# Step 1: Read locus list
try:
    with open('../nSSR_LocusList.txt', 'r') as f:
        locus_list = [line.strip() for line in f if line.strip()]
except FileNotFoundError:
    print("Error: ../nSSR_LocusList.txt not found.")
    sys.exit(1)

# Step 2: Read dictionary CSV and collect patterns per locus
locus_to_patterns = defaultdict(list)
try:
    with open('../new_output/AlleleInfo_dictionary.csv', newline='') as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader)
        for row in reader:
            if len(row) < 4: 
                continue
            locus = row[-1]     
            
            # --- THIS IS THE FIX ---
            # Use row[1] (Name) to match Step4, not row[0] (ID)
            label = row[1]      
            # --- END OF FIX ---
            
            sequence = row[2]   # Sequence is 3rd column (index 2)
            if locus in locus_list:
                locus_to_patterns[locus].append((label, sequence))
except FileNotFoundError:
    print("Error: ../new_output/AlleleInfo_dictionary.csv not found.")
    sys.exit(1)


all_rows = []
all_header = ['Locus', 'Allele1', 'Allele2', 'Value']

# Step 3: For each locus, build and write the comparison matrix
for locus in locus_list:
    patterns = sorted(set(locus_to_patterns[locus]), key=lambda x: x[0])
    if not patterns:
        print(f'Warning: No patterns found for locus {locus} in dictionary. Skipping.')
        continue
        
    labels = [p[0] for p in patterns]
    sequences = [p[1] for p in patterns]
    parsed_patterns = [parse_repeat_pattern(seq) for seq in sequences]
    
    matrix = []
    for i, pat_bench in enumerate(parsed_patterns):
        row = []
        for j, pat_comp in enumerate(parsed_patterns):
            score = compare_patterns_asym(pat_bench, pat_comp, SCORE_MINUS_1, SCORE_OTHER_DIFF)
            row.append(score)
        matrix.append(row)
        
    output_file = f'../new_output/table/{locus}_SequenceComparison.csv'
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w', newline='') as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow(['Allele1', 'Allele2', 'Value'])
        for i, label1 in enumerate(labels):
            for j, label2 in enumerate(labels):
                writer.writerow([label1, label2, matrix[i][j]])
                all_rows.append([locus, label1, label2, matrix[i][j]])
    print(f'Wrote sequence comparison matrix for {locus} to {output_file}')

# Write combined output file
combined_output_file = '../new_output/table/All_SequenceComparisons.csv'
with open(combined_output_file, 'w', newline='') as allcsv:
    writer = csv.writer(allcsv)
    writer.writerow(all_header)
    writer.writerows(all_rows)
print(f'WWrote combined sequence comparison matrix to {combined_output_file}')