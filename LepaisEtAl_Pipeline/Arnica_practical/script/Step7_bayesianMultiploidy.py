#!/usr/bin/env python3
"""
Step 7: Bayesian Genotype Likelihood Calculator
- Reads the fused file (Alleles + Stutter Comparisons).
- Calculates probability of genotypes for a fixed Ploidy.
"""

import csv
import sys
import re
import math
from collections import Counter
from datetime import datetime
import os

# --- Parameter Part ---
# Usage: python Step7_bayesianMultiploidy.py [PLOIDY] [MODE]
# Example: python Step7_bayesianMultiploidy.py 2 human

# Defaults
PLOIDY = 2
MODE_INPUT = 'human'

# Argument Parsing
args = sys.argv[1:]

# 1. Parse Ploidy (First number found)
for arg in args:
    if arg.isdigit():
        PLOIDY = int(arg)
        break

# 2. Parse Mode (String matching known modes)
KNOWN_MODES = {'human', 'genalex', 'original'}
for arg in args:
    if arg.lower() in KNOWN_MODES:
        MODE_INPUT = arg.lower()
        break

print(f"✅ Configuration:")
print(f"  - Ploidy: {PLOIDY}")
print(f"  - Mode:   {MODE_INPUT}")

# --- File Paths ---
# INPUT: Output from Step 6a
INPUT_CSV = f"../new_output/Step6a_AlleleInfo_named_fused_{MODE_INPUT}.csv"

# OUTPUT: Likelihoods file
OUTPUT_CSV = f"../new_output/Step7_likelihoods_Ploidy{PLOIDY}_{MODE_INPUT}.csv"

# Pseudo-count scaling factor for multinomial likelihood
N_SCALE = 1000  

# --- Helper Functions ---

def _softmax_list(log_likelihoods):
    """Normalizes a list of log-likelihoods into probabilities."""
    if not log_likelihoods:
        return []
    max_log = max(log_likelihoods)
    if max_log == float("-inf"):
        return [1.0 / len(log_likelihoods)] * len(log_likelihoods)
    exps = [math.exp(logL - max_log) if logL > float("-inf") else 0.0 for logL in log_likelihoods]
    sum_exps = sum(exps)
    if sum_exps == 0.0:
        return [1.0 / len(log_likelihoods)] * len(log_likelihoods)
    return [e / sum_exps for e in exps]

def _log_term(count, prob):
    if prob <= 0.0:
        return 0.0 if count == 0 else float("-inf")
    return count * math.log(prob)

def get_genotype_hypotheses(ploidy):
    """Generates fixed list of hypotheses (e.g., AA, AB, BB for diploid)."""
    if ploidy == 1:
        hypotheses = [('A',)]
    elif ploidy == 2:
        hypotheses = [('A', 'A'), ('A', 'B')]
    elif ploidy == 3:
        hypotheses = [('A', 'A', 'A'), ('A', 'A', 'B'), ('A', 'B', 'C')]
    elif ploidy == 4:
        hypotheses = [
            ('A', 'A', 'A', 'A'), ('A', 'A', 'A', 'B'), 
            ('A', 'A', 'B', 'B'), ('A', 'A', 'B', 'C'), 
            ('A', 'B', 'C', 'D')
        ]
    else:
        print(f"❌ Error: Ploidy {ploidy} logic not implemented.")
        sys.exit(1)

    names = ["Prob_" + "".join(h) for h in hypotheses]
    return hypotheses, names

def calculate_likelihood(hypothesis, allele_map, observed_data, stutter_map, ploidy):
    # 1. Translate hypothesis placeholders to real names
    try:
        true_genotype = tuple(allele_map[p] for p in hypothesis)
    except KeyError:
        return float("-inf")
        
    genotype_counts = Counter(true_genotype)
    expected_freqs = {}
    all_observed = list(observed_data.keys())
    
    # 2. Calculate Expected Frequencies
    for read_allele in all_observed:
        expected_p = 0.0
        for true_allele, dose in genotype_counts.items():
            prob_true_given_H = dose / ploidy
            
            if read_allele == true_allele:
                # P(read A | true A) = 1 - sum(stutters)
                sum_stutters = 0.0
                for other in all_observed:
                    if other != true_allele:
                        sum_stutters += stutter_map.get((true_allele, other), 0.0)
                prob_read_given_true = max(0.0, 1.0 - sum_stutters)
            else:
                # P(read B | true A)
                prob_read_given_true = stutter_map.get((true_allele, read_allele), 0.0)
            
            expected_p += prob_read_given_true * prob_true_given_H
        expected_freqs[read_allele] = expected_p

    # 3. Multinomial Log-Likelihood
    log_likelihood = 0.0
    sum_expected = sum(expected_freqs.values())
    if sum_expected == 0.0: return float("-inf") 

    for allele, data in observed_data.items():
        # Normalize expected prob to ensure sum=1 within the observed set
        normalized_prob = expected_freqs[allele] / sum_expected
        log_likelihood += _log_term(data['count'], normalized_prob)

    return log_likelihood

# --- Main Execution ---

def main():
    hypotheses, hypothesis_names = get_genotype_hypotheses(PLOIDY)
    
    try:
        infile = open(INPUT_CSV, 'r', newline='')
    except FileNotFoundError:
        print(f"❌ Error: Input file not found: {INPUT_CSV}")
        print(f"👉 Solution: Check if Step 6a ran successfully with mode '{MODE_INPUT}'.")
        return

    reader = csv.reader(infile)
    try:
        header = next(reader)
    except StopIteration:
        return

    # Dynamically find columns
    allele_pattern = re.compile(r'^Allele(\d+)$')
    freq_pattern = re.compile(r'^Frequency(\d+)$')
    comp_pattern = re.compile(r'^Comp(\d+)_(\d+)$')
    
    allele_cols = {}
    freq_cols = {}
    comp_cols = {}
    
    # Check for required base columns
    try:
        # Some headers might be 'Locus' instead of 'Loci', let's be flexible
        col_map = {name: i for i, name in enumerate(header)}
        
        idx_serial = col_map.get('Serial', 0) # Default to 0 if missing, likely unused
        idx_locus = col_map.get('Locus') if 'Locus' in col_map else col_map.get('Loci')
        idx_sample = col_map.get('Sample') if 'Sample' in col_map else col_map.get('Sample_Name')

        if idx_locus is None: raise ValueError("Column 'Locus' or 'Loci' missing")
        if idx_sample is None: raise ValueError("Column 'Sample' missing")

    except ValueError as e:
        print(f"❌ Error in header: {e}")
        return

    for i, col_name in enumerate(header):
        if allele_pattern.match(col_name): allele_cols[col_name] = i
        elif freq_pattern.match(col_name): freq_cols[col_name] = i
        elif comp_pattern.match(col_name): comp_cols[col_name] = i
    
    sorted_allele_cols = sorted(allele_cols.items(), key=lambda x: int(allele_pattern.match(x[0]).group(1)))
    sorted_freq_cols = sorted(freq_cols.items(), key=lambda x: int(freq_pattern.match(x[0]).group(1)))
    
    # Verify we have frequencies
    if not sorted_freq_cols:
         print("❌ Error: No 'FrequencyX' columns found. Bayesian model requires allele frequencies.")
         sys.exit(1)

    placeholder_chars = [chr(65 + i) for i in range(len(sorted_allele_cols))]
    
    print(f"Processing data...")

    with open(OUTPUT_CSV, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        
        # Header
        allele_out_cols = [f"Allele{i+1}_Name" for i in range(PLOIDY)]
        new_header = ['Locus', 'Sample', 'Ploidy'] + allele_out_cols + hypothesis_names
        writer.writerow(new_header)
        
        count_processed = 0
        for row in reader:
            # Safe row access
            def get_val(idx): return row[idx] if idx < len(row) else None
            
            # 1. Parse all allele names
            all_row_names = {}
            for col, idx in sorted_allele_cols:
                val = get_val(idx)
                if val: val = val.strip()
                # Basic validity check
                if val and val.lower() not in ('n/a', 'no data', '', 'nan'):
                    all_row_names[col] = val
                else:
                    all_row_names[col] = None

            # 2. Prepare Output Names (The top N valid alleles found)
            # Or strictly Allele1..N columns? 
            # Let's use the columns corresponding to Ploidy for consistency with input
            out_names = []
            for i in range(PLOIDY):
                if i < len(sorted_allele_cols):
                    col_name = sorted_allele_cols[i][0]
                    name = all_row_names.get(col_name, 'N/A')
                    out_names.append(name if name else 'N/A')
                else:
                    out_names.append('N/A')

            base_data = [
                get_val(idx_locus),
                get_val(idx_sample),
                PLOIDY
            ] + out_names

            # 3. Build Contiguous Block of Valid Data
            observed_data = {}
            valid_placeholders = {}

            # Pair Allele col with Frequency col
            for i, ((a_col, a_idx), (f_col, f_idx)) in enumerate(zip(sorted_allele_cols, sorted_freq_cols)):
                name = all_row_names.get(a_col)
                
                # Check Frequency
                freq_val = 0.0
                freq_valid = False
                f_str = get_val(f_idx)
                if name and f_str:
                    try:
                        freq_val = float(f_str)
                        freq_valid = True
                    except ValueError:
                        pass
                
                # STOP condition: First invalid allele or invalid frequency stops the chain
                if not name or not freq_valid:
                    break
                
                observed_data[name] = {
                    'freq': freq_val,
                    'count': int(round(freq_val * N_SCALE))
                }
                valid_placeholders[placeholder_chars[i]] = name

            if not observed_data:
                # No valid data start
                writer.writerow(base_data + [0.0] * len(hypotheses))
                continue

            # 4. Build Stutter Map from Comp columns
            stutter_map = {}
            for col, idx in comp_cols.items():
                val = get_val(idx)
                if val is None: continue
                
                m = comp_pattern.match(col)
                n1, n2 = m.group(1), m.group(2)
                
                name1 = all_row_names.get(f'Allele{n1}')
                name2 = all_row_names.get(f'Allele{n2}')
                
                if name1 and name2:
                    try:
                        stutter_map[(name1, name2)] = float(val)
                    except ValueError:
                        stutter_map[(name1, name2)] = 0.0

            # 5. Run Hypotheses
            log_likelihoods = []
            for H in hypotheses:
                logL = calculate_likelihood(H, valid_placeholders, observed_data, stutter_map, PLOIDY)
                log_likelihoods.append(logL)

            # 6. Write Result
            probs = _softmax_list(log_likelihoods)
            
            # Handle all -inf case
            if log_likelihoods and all(L == float("-inf") for L in log_likelihoods):
                 writer.writerow(base_data + [0.0] * len(hypotheses))
            else:
                 writer.writerow(base_data + probs)
            
            count_processed += 1

    infile.close()
    print(f"✅ Finished! Processed {count_processed} rows.")
    print(f"✅ Output saved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()