import csv
import re
import sys 
from pathlib import Path
from collections import defaultdict

# --- Parameter Part ---
# Usage: python Step6_StutterTable.py [SCORE_MINUS_1] [SCORE_OTHER_DIFF] [MODE]

# Defaults
SCORE_MINUS_1 = 0.5
SCORE_OTHER_DIFF = 0.1
MODE_INPUT = 'human'

args = sys.argv[1:]

# 1. Parse Scores (Float detection)
if len(args) >= 1 and args[0].replace('.', '', 1).isdigit():
    SCORE_MINUS_1 = float(args[0])
    if len(args) >= 2 and args[1].replace('.', '', 1).isdigit():
        SCORE_OTHER_DIFF = float(args[1])

# 2. Parse Mode (String detection)
MODE_MAP = {
    'human': 1,      # Human_ID (Arm01_83)
    'genalex': 2,    # GenAlEx_ID (83)
    'original': 3    # Original_AlleleSeqCode (101)
}

for arg in args:
    if arg.lower() in MODE_MAP:
        MODE_INPUT = arg.lower()
        break

LABEL_COL_INDEX = MODE_MAP[MODE_INPUT]

print(f"✅ Configuration:")
print(f"  - Mode: {MODE_INPUT} (Using Column {LABEL_COL_INDEX} for names)")
print(f"  - Score (-1 step): {SCORE_MINUS_1}")
print(f"  - Score (+/- 1, 2 steps): {SCORE_OTHER_DIFF}")
# --- End Parameter Part ---


def parse_repeat_pattern(s):
    # Extracts pairs like ('CATA', '13') from "TGT...CATA(13)..."
    return re.findall(r'([A-Za-z]+)\((\d+)\)', s)


def compare_patterns_asym(bench_parsed, comp_parsed, score_val_minus_1, score_val_other):
    """
    Compares two parsed patterns.
    bench_parsed: list of (motif, count) tuples
    """
    # 1. Exact Match
    if bench_parsed == comp_parsed:
        return 1.0
        
    # 2. Structural Mismatch check
    # If they don't have the same number of motifs/regions, we can't compare repeats
    if len(bench_parsed) != len(comp_parsed):
        return 0.0

    diffs = []
    
    # Zip through the motifs in parallel
    for (motif1, count1), (motif2, count2) in zip(bench_parsed, comp_parsed):
        # If the motif sequence itself is different (e.g. CATA vs GATA), it's a mismatch
        if motif1 != motif2:
            return 0.0 
            
        # If counts differ, record the difference
        c1 = int(count1)
        c2 = int(count2)
        if c1 != c2:
            diff = c1 - c2
            diffs.append(diff)
    
    # 3. Analyze Differences
    # We only score "Stutter" if there is EXACTLY ONE difference in the repeat counts
    if len(diffs) == 1:
        d = diffs[0]
        # Logic: -1 means Comp has 1 more repeat than Bench (Common stutter)
        if d == -1:
            return score_val_minus_1
        # Other close steps
        elif d in [1, 2, -2]:
            return score_val_other
            
    return 0.0


# Step 1: Read locus list
try:
    with open('../nSSR_LocusList.txt', 'r') as f:
        locus_list = [line.strip() for line in f if line.strip()]
except FileNotFoundError:
    print("❌ Error: ../nSSR_LocusList.txt not found.")
    sys.exit(1)

# Step 2: Read dictionary CSV and collect patterns per locus
locus_to_data = defaultdict(list)
dict_path = '../new_output/Step3b_Patched_AlleleInfo_dictionary.csv'

try:
    with open(dict_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader) # Skip Header
        for row in reader:
            # Check row length (Needs at least 6 columns to reach Repeat_Pattern)
            if len(row) < 6: 
                continue
            
            # Header Ref: 
            # 0:Fused, 1:Human, 2:GenAlEx, 3:Orig, 4:Seq, 5:Repeat_Pattern, 6:Len, 7:Locus
            
            locus = row[7] if len(row) > 7 else row[-1] # Robust fetch
            
            label = row[LABEL_COL_INDEX] # The name (e.g. Arm01_83)
            pattern_str = row[5]         # The Repeat Pattern (e.g. ...CATA(13)...)
            
            if locus in locus_list:
                locus_to_data[locus].append((label, pattern_str))
                
except FileNotFoundError:
    print(f"❌ Error: {dict_path} not found.")
    print("👉 Please run Step 3 first.")
    sys.exit(1)


all_rows = []
# Step 6a usually expects just: Allele1, Allele2, Value
all_header = ['Allele1', 'Allele2', 'Value'] 

print("Processing loci patterns...")

# Step 3: Loop through loci
for locus in locus_list:
    # Get items for this locus: List of (Label, PatternString)
    items = locus_to_data[locus]
    if not items:
        continue
        
    # Sort items by Label to keep matrix deterministic
    items.sort(key=lambda x: x[0])
    
    labels = [x[0] for x in items]
    pattern_strings = [x[1] for x in items]
    
    # Pre-parse patterns to avoid re-parsing inside the loop (Optimization)
    parsed_patterns = [parse_repeat_pattern(p) for p in pattern_strings]
    
    matrix = []
    
    # Double loop to compare every allele against every allele
    for i in range(len(labels)):
        row_vals = []
        for j in range(len(labels)):
            score = compare_patterns_asym(
                parsed_patterns[i], 
                parsed_patterns[j], 
                SCORE_MINUS_1, 
                SCORE_OTHER_DIFF
            )
            row_vals.append(score)
            
            # Add to the big list
            all_rows.append([labels[i], labels[j], score])
            
        matrix.append(row_vals)
        
    # Optional: Write individual locus file
    output_file = f'../new_output/table/{locus}_SequenceComparison.csv'
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w', newline='') as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow(['Allele1', 'Allele2', 'Value'])
        for i, label1 in enumerate(labels):
            for j, label2 in enumerate(labels):
                writer.writerow([label1, label2, matrix[i][j]])

# Write combined output file
combined_output_file = '../new_output/table/Step6_All_SequenceComparisons.csv'
with open(combined_output_file, 'w', newline='') as allcsv:
    writer = csv.writer(allcsv)
    writer.writerow(all_header)
    writer.writerows(all_rows)

print(f"✅ Generated regex-based comparisons for {len(locus_list)} loci.")
print(f"✅ Saved combined table to: {combined_output_file}")