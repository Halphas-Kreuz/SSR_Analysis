import csv
import re
from pathlib import Path
from collections import defaultdict

def parse_repeat_pattern(s):
    # Returns a list of (motif, repeat_count) for all (motif)(n) in the string
    return re.findall(r'([A-Za-z]+)\((\d+)\)', s)

def compare_patterns_asym(bench, comp):
    # If patterns are identical, return 1
    if bench == comp:
        return 1.0
    # If all motifs are the same and only one repeat count differs, check direction
    if len(bench) == len(comp):
        diffs = []
        diff_idx = -1
        for i, ((m1, n1), (m2, n2)) in enumerate(zip(bench, comp)):
            if m1 != m2:
                return 0.0
            if n1 != n2:
                diffs.append((i, int(n1) - int(n2)))
                diff_idx = i
        if len(diffs) == 1:
            diff = diffs[0][1]
            if diff == 1:
                return 0.5
            elif diff == -1:
                return 0.1
    return 0.0

# Step 1: Read locus list
with open('../nSSR_LocusList.txt', 'r') as f:
    locus_list = [line.strip() for line in f if line.strip()]

# Step 2: Read dictionary CSV and collect patterns per locus
locus_to_patterns = defaultdict(list)
with open('../new_output/AlleleInfo_dictionary.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    header = next(reader)
    for row in reader:
        locus = row[-1]
        label = row[0]      # Use as row/column name
        sequence = row[1]   # Use for comparison
        if locus in locus_list:
            locus_to_patterns[locus].append((label, sequence))

all_rows = []
all_header = ['Locus', 'Allele1', 'Allele2', 'Value']

# Step 3: For each locus, build and write the comparison matrix
for locus in locus_list:
    # Get unique (label, sequence) pairs
    patterns = sorted(set(locus_to_patterns[locus]), key=lambda x: x[0])
    labels = [p[0] for p in patterns]
    sequences = [p[1] for p in patterns]
    parsed_patterns = [parse_repeat_pattern(seq) for seq in sequences]
    matrix = []
    for i, pat_bench in enumerate(parsed_patterns):
        row = []
        for j, pat_comp in enumerate(parsed_patterns):
            score = compare_patterns_asym(pat_bench, pat_comp)
            row.append(score)
        matrix.append(row)
    # Write CSV for this locus
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
print(f'Wrote combined sequence comparison matrix to {combined_output_file}')