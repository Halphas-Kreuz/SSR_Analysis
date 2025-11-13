#!/usr/bin/env python3
"""
Bayesian Genotype Likelihood Calculator (Fixed Ploidy Model) - v4_fixed

Modifications:
- (v4_fixed) Corrected a typo in the stutter_map 'except' block (namef1 -> name1).
- (v4) Implements "contiguous block" logic. The script now stops
       collecting observed data at the first invalid allele OR
       invalid frequency. Everything after this point is ignored.
- (v3) If ALL hypotheses fail, output all 0.0 probabilities.
- (v2) Outputs allele names (Allele1_Name, etc.) in the row.
- (v2) If Allele1 is invalid, outputs all 0.0 probabilities.
"""

import csv
import sys
import re
import math
from collections import Counter, defaultdict
from itertools import combinations_with_replacement
from datetime import datetime

# --- Parameters ---
try:
    # Read ploidy from command line (e.g., python script.py 3)
    PLOIDY = int(sys.argv[1])
    if PLOIDY < 1:
        raise ValueError()
except (IndexError, ValueError):
    print("No/invalid ploidy provided. Defaulting to diploid (Ploidy = 2).")
    PLOIDY = 2

print(f"✅ Running Bayesian model for fixed PLOIDY = {PLOIDY}")

# Pseudo-count scaling factor for multinomial likelihood
N_SCALE = 1000  
INPUT_CSV = "../new_output/AlleleInfo_named_fused.csv"
# OUTPUT_CSV = f"../new_output/likelihoods_Ploidy{PLOIDY}_{datetime.now():%Y%m%d_%H%M}.csv"
OUTPUT_CSV = f"../new_output/likelihoods.csv"

# --- Helper Functions ---

def _softmax_list(log_likelihoods: list[float]) -> list[float]:
    """
    Normalizes a list of log-likelihoods into probabilities.
    """
    if not log_likelihoods:
        return []
    
    max_log = max(log_likelihoods)
    if max_log == float("-inf"):
        # This case should be handled by the main loop,
        # but as a safety, return uniform.
        return [1.0 / len(log_likelihoods)] * len(log_likelihoods)
        
    # Calculate exponentials relative to the max (for numerical stability)
    exps = [math.exp(logL - max_log) if logL > float("-inf") else 0.0 for logL in log_likelihoods]
    sum_exps = sum(exps)
    
    if sum_exps == 0.0:
        # Fallback in an unlikely underflow case
        return [1.0 / len(log_likelihoods)] * len(log_likelihoods)
        
    return [e / sum_exps for e in exps]

def _log_term(count: int, prob: float) -> float:
    """Calculates a single term of the multinomial log-likelihood."""
    if prob <= 0.0:
        # If prob is 0, log-likelihood is -inf unless count is also 0
        return 0.0 if count == 0 else float("-inf")
    return count * math.log(prob)

def get_genotype_hypotheses(ploidy: int) -> tuple[list[tuple], list[str]]:
    """
    Generates the fixed list of hypotheses and their names for a given ploidy.
    'A' = Allele1, 'B' = Allele2, etc.
    """
    if ploidy == 1:
        alleles = ['A']
        hypotheses = [('A',)]
    elif ploidy == 2:
        alleles = ['A', 'B']
        hypotheses = [('A', 'A'), ('A', 'B')]
    elif ploidy == 3:
        alleles = ['A', 'B', 'C']
        hypotheses = [('A', 'A', 'A'), ('A', 'A', 'B'), ('A', 'B', 'C')]
    elif ploidy == 4:
        alleles = ['A', 'B', 'C', 'D']
        hypotheses = [
            ('A', 'A', 'A', 'A'), 
            ('A', 'A', 'A', 'B'), 
            ('A', 'A', 'B', 'B'), 
            ('A', 'A', 'B', 'C'), 
            ('A', 'B', 'C', 'D')
        ]
    else:
        raise ValueError(f"Ploidy {ploidy} is not supported by this script's fixed logic.")

    # Create clean names, e.g., "Prob_AAB"
    names = ["Prob_" + "".join(h) for h in hypotheses]
    return hypotheses, names

def calculate_likelihood(
    hypothesis: tuple,            # e.g., ('A', 'A', 'B')
    allele_map: dict,              # e.g., {'A': 'Allele1_Name', 'B': 'Allele2_Name', ...}
    observed_alleles_data: dict,   # e.g., {'MSCB_89': {'freq': 0.6, 'count': 600}, ...}
    stutter_map: dict,             # e.g., {('MSCB_89', 'MSCB_88'): 0.1, ...}
    ploidy: int
) -> float:
    """
    This is the new Bayesian core.
    Calculates the log-likelihood of one hypothesis given the observed data.
    """
    
    # 1. Translate hypothesis placeholders ('A', 'A', 'B') into real allele names
    #    e.g., ('MSCB_89', 'MSCB_89', 'MSCB_88')
    try:
        true_genotype = tuple(allele_map[placeholder] for placeholder in hypothesis)
    except KeyError:
        # This happens if, e.g., Ploidy=4 and H=ABCD, but Allele3 is "N/A"
        # and was not included in the 'allele_map'
        # This hypothesis is impossible.
        return float("-inf")
        
    # Count the "dose" of each true allele in the hypothesis
    # e.g., {'MSCB_89': 2, 'MSCB_88': 1}
    genotype_counts = Counter(true_genotype)

    # 2. Calculate the *Expected* frequency for *every* observed allele
    expected_freqs = {}
    
    # Get a list of all *valid* observed alleles (e.g., ['MSCB_89', 'MSCB_88', ...])
    all_observed_alleles = list(observed_alleles_data.keys())
    
    for read_allele in all_observed_alleles:
        expected_p = 0.0
        # This is the core math:
        # P(read_allele | H) = sum over true_alleles [ P(read_allele | true_allele) * P(true_allele | H) ]
        
        for true_allele, dose in genotype_counts.items():
            # P(true_allele | H) is the dose / ploidy
            prob_true_given_H = dose / ploidy
            
            # P(read_allele | true_allele) is the stutter rate
            if read_allele == true_allele:
                # This is P(read A | true A)
                # It's 1.0 minus the sum of all stutters *from* A
                sum_stutters_from_true = 0.0
                for other_allele in all_observed_alleles:
                    if other_allele != true_allele:
                        sum_stutters_from_true += stutter_map.get((true_allele, other_allele), 0.0)
                
                # Clip at 0 in case stutter sums > 1.0
                prob_read_given_true = max(0.0, 1.0 - sum_stutters_from_true)
            
            else:
                # This is P(read B | true A), a direct lookup
                prob_read_given_true = stutter_map.get((true_allele, read_allele), 0.0)
            
            expected_p += prob_read_given_true * prob_true_given_H
        
        expected_freqs[read_allele] = expected_p

    # 3. Calculate the final multinomial log-likelihood
    # We compare the list of *expected* frequencies with the *observed* counts
    log_likelihood = 0.0
    
    # Re-normalize expected frequencies just in case they don't sum to 1
    sum_expected = sum(expected_freqs.values())
    if sum_expected == 0.0:
        # This can happen if all probs are 0
        return float("-inf") 

    for allele_name, data in observed_alleles_data.items():
        observed_count = data['count']
        normalized_expected_prob = expected_freqs[allele_name] / sum_expected
        
        log_likelihood += _log_term(observed_count, normalized_expected_prob)

    return log_likelihood

# --- Main Execution ---

def main():
    
    # 1. Generate the fixed hypotheses for this ploidy run
    hypotheses, hypothesis_names = get_genotype_hypotheses(PLOIDY)
    print(f"Testing {len(hypotheses)} hypotheses: {hypothesis_names}")
    
    try:
        infile = open(INPUT_CSV, 'r', newline='')
    except FileNotFoundError:
        print(f"FATAL ERROR: Input file not found.\n{INPUT_CSV}", file=sys.stderr)
        return

    reader = csv.reader(infile)
    
    try:
        header = next(reader)
    except StopIteration:
        print(f"FATAL ERROR: Input file is empty.\n{INPUT_CSV}", file=sys.stderr)
        return

    # 2. Dynamically find all Allele, Frequency, and Comparison columns
    allele_pattern = re.compile(r'^Allele(\d+)$')
    freq_pattern = re.compile(r'^Frequency(\d+)$')
    comp_pattern = re.compile(r'^Comp(\d+)_(\d+)$')
    
    allele_cols = {} # 'Allele1': index
    freq_cols = {}   # 'Frequency1': index
    comp_cols = {}   # 'Comp1_2': index
    
    try:
        base_cols = {
            'Serial': header.index('Serial'),
            'Loci': header.index('Loci'),
            'Sample': header.index('Sample')
        }
    except ValueError as e:
        print(f"FATAL ERROR: Input file missing required column. {e}", file=sys.stderr)
        return

    for i, col_name in enumerate(header):
        if allele_pattern.match(col_name):
            allele_cols[col_name] = i
        elif freq_pattern.match(col_name):
            freq_cols[col_name] = i
        elif comp_pattern.match(col_name):
            comp_cols[col_name] = i
    
    sorted_allele_cols = sorted(allele_cols.items(), key=lambda item: int(allele_pattern.match(item[0]).group(1)))
    sorted_freq_cols = sorted(freq_cols.items(), key=lambda item: int(freq_pattern.match(item[0]).group(1)))
    
    # Placeholder map, e.g., {'A': 'Allele1', 'B': 'Allele2', ...}
    # We build this *once* for all alleles found
    placeholder_map = {chr(65 + i): name for i, (name, idx) in enumerate(sorted_allele_cols)}
    placeholder_chars = [chr(65 + i) for i in range(len(sorted_allele_cols))] # ['A', 'B', 'C', ...]
    
    print(f"Found {len(allele_cols)} Allele, {len(freq_cols)} Freq, and {len(comp_cols)} Comp columns.")

    # 3. Open output file and write the new "wide" header
    with open(OUTPUT_CSV, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        
        # Add AlleleX_Name columns to header
        allele_name_cols = []
        for i in range(PLOIDY):
            allele_name_cols.append(f"Allele{i+1}_Name")
        
        new_header = [
            'Serial', 'Loci', 'Sample', 'Ploidy'
        ] + allele_name_cols + hypothesis_names
        writer.writerow(new_header)
        
        # 4. Process each row from the input file
        for row in reader:
            
            # --- a. Parse all data for this row ---
            
            # First, get *all* allele names, even invalid ones.
            # This is for the stutter map, which needs all pairs.
            all_row_allele_names = {} # {'Allele1': 'MSCB_89', 'Allele2': 'N/A', ...}
            for col_name, idx in sorted_allele_cols:
                if idx >= len(row): # Handle rows that are too short
                    all_row_allele_names[col_name] = None
                    continue
                allele_name = row[idx].strip()
                if allele_name and allele_name.lower() not in ('n/a', 'no data', 'other sequences'):
                    all_row_allele_names[col_name] = allele_name
                else:
                    all_row_allele_names[col_name] = None # Mark as invalid

            # Get the allele names for the output row (up to PLOIDY)
            output_allele_names = []
            for i, (col_name, idx) in enumerate(sorted_allele_cols):
                if i >= PLOIDY:
                    break
                name = all_row_allele_names.get(col_name) # Get pre-parsed name
                name = name if name is not None else 'N/A'
                output_allele_names.append(name)
            
            # Pad with 'N/A' if row had fewer allele columns than PLOIDY
            while len(output_allele_names) < PLOIDY:
                output_allele_names.append('N/A')

            # Build the base of the output row
            base_data = [
                row[base_cols['Serial']],
                row[base_cols['Loci']],
                row[base_cols['Sample']],
                PLOIDY
            ] + output_allele_names

            # --- "Contiguous Block" Logic (v4) ---
            # Build the *contiguous* block of valid data for the model
            
            observed_alleles_data = {} # {'MSCB_89': {'freq': 0.5, 'count': 500}, ...}
            valid_allele_placeholders = {} # {'A': 'MSCB_89', 'B': 'MSCB_88'}

            for i, ((allele_col_name, allele_idx), (freq_col_name, freq_idx)) in enumerate(zip(sorted_allele_cols, sorted_freq_cols)):
                
                allele_name = all_row_allele_names.get(allele_col_name) # Get pre-parsed name
                placeholder = placeholder_chars[i]

                # Check for frequency validity
                freq_is_valid = False
                freq = 0.0
                if allele_name: # Only check freq if allele is valid
                    try:
                        if freq_idx >= len(row): # Row is too short
                            pass # freq_is_valid remains False
                        else:
                            freq = float(row[freq_idx])
                            freq_is_valid = True # Success!
                    except (ValueError, TypeError):
                        pass # freq_is_valid remains False

                # --- This is the "STOP" condition ---
                if not allele_name or not freq_is_valid:
                    # We found an invalid allele OR invalid frequency.
                    # Stop processing here. Do not look at Allele4, etc.
                    break 
                
                # If we are here, both allele and freq are valid.
                observed_alleles_data[allele_name] = {
                    'freq': freq,
                    'count': int(round(freq * N_SCALE))
                }
                valid_allele_placeholders[placeholder] = allele_name

            # --- (END "Contiguous Block" Logic) ---


            if not observed_alleles_data:
                # No valid data for this row (e.g., Allele1 was invalid)
                writer.writerow(base_data + [0.0] * len(hypotheses))
                continue

            # Build the full stutter map (uses all_row_allele_names)
            stutter_map = {}
            for col_name, idx in comp_cols.items():
                if idx >= len(row): continue # Row is too short
                
                match = comp_pattern.match(col_name)
                num1, num2 = match.group(1), match.group(2)
                
                allele1_col = f'Allele{num1}'
                allele2_col = f'Allele{num2}'

                name1 = all_row_allele_names.get(allele1_col)
                name2 = all_row_allele_names.get(allele2_col)
                
                if name1 and name2: # Only add to map if both alleles are valid
                    try:
                        stutter_map[(name1, name2)] = float(row[idx])
                    except (ValueError, TypeError):
                        # --- THIS IS THE FIX ---
                        stutter_map[(name1, name2)] = 0.0

            # --- b. Run all hypotheses ---
            log_likelihoods = []
            for H in hypotheses:
                # Pass the map of *contiguous* valid placeholders
                logL = calculate_likelihood(
                    hypothesis=H,
                    allele_map=valid_allele_placeholders,
                    observed_alleles_data=observed_alleles_data,
                    stutter_map=stutter_map,
                    ploidy=PLOIDY
                )
                log_likelihoods.append(logL)

            # --- c. Normalize and write results ---
            probabilities = _softmax_list(log_likelihoods)
            
            # If all hypotheses failed (all -inf), write all 0.0s
            if log_likelihoods and all(L == float("-inf") for L in log_likelihoods):
                 writer.writerow(base_data + [0.0] * len(hypotheses))
            else:
                 writer.writerow(base_data + probabilities)

    # 5. Done
    infile.close()
    print(f"\nSuccess! Wrote {reader.line_num - 1} rows to '{OUTPUT_CSV}'.")


if __name__ == "__main__":
    main()