#!/usr/bin/env python3
"""
Automatic Ploidy Pattern Likelihood Calculator (Fast Version)

This script implements the optimized logic:
1.  It automatically detects ploidy (n) from the 'Frequency<n>' columns.
2.  It generates all 'n' partitions (patterns) (e.g., [4], [3, 1], [2, 2]).
3.  It assumes the input data 'Frequency1', 'Frequency2', etc., are
    sorted in descending order of coverage.
4.  For each pattern (e.g., [3, 1]), it performs a *direct assignment*
    (3 parts to Allele1, 1 part to Allele2) and calculates a
    single likelihood score.
5.  It outputs the single best-fitting pattern for each row.
"""

from __future__ import annotations
import math
from datetime import datetime
import pandas as pd
import numpy as np
import sys
from collections import defaultdict
# We no longer need 'itertools.permutations'

# --- Parameters ---
N_SCALE = 1000  # Pseudo-count scaling factor
INPUT_CSV = "../new_output/AlleleInfo_named_fused.csv"
OUTPUT_CSV = f"../new_output/best_patterns_{datetime.now():%Y%m%d_%HM%S}.csv"

# Set a higher recursion depth for safely handling larger n values
sys.setrecursionlimit(2000)

# ##########################################################################
# LOGIC FROM: ploidy_calculator.py (Unchanged)
# ##########################################################################

def find_all_partitions(n):
    """
    Finds all integer partitions of n using a recursive (backtracking) approach.
    """
    all_partitions = []
    
    def get_partitions_recursive(target_sum, max_part, current_partition):
        if target_sum == 0:
            all_partitions.append(list(current_partition))
            return

        for i in range(min(max_part, target_sum), 0, -1):
            current_partition.append(i)
            get_partitions_recursive(target_sum - i, i, current_partition)
            current_partition.pop()

    get_partitions_recursive(n, n, [])
    return all_partitions

# ##########################################################################
# LIKELIHOOD FUNCTIONS (Unchanged)
# ##########################################################################

def _log_term(count: int, prob: float) -> float:
    """Return count*log(prob), with exact-zero handling."""
    if prob <= 0.0:
        return 0.0 if count == 0 else float("-inf")
    return count * math.log(prob)


def _log_likelihood_list(counts: list[int], cD: int, probs: list[float], pD: float) -> float:
    """
    Calculates the multinomial log-likelihood for a list of counts.
    """
    terms = [_log_term(c, p) for c, p in zip(counts, probs)]
    terms.append(_log_term(cD, pD))
    
    if any(t == float("-inf") for t in terms):
        if cD > 0 and pD == 0.0:
            return float("-inf")
        if any(c > 0 and p == 0.0 for c, p in zip(counts, probs)):
            return float("-inf")
            
    return sum(terms)


def get_emission_probs_list(genotype_counts: list[int], n: int, fD: float) -> tuple[list[float], float]:
    """
    Calculates the expected P(f1), P(f2)...P(fn) and P(D) for a
    given genotype, assuming e_cross = 0.
    """
    expected_freqs = [count / n for count in genotype_counts]
    probs = [f * (1.0 - fD) for f in expected_freqs]
    pD = fD
    
    return (probs, pD)

# ##########################################################################
# NEW COMBINED LOGIC (Simplified)
# ##########################################################################

def compute_best_pattern(frequencies: list[float], n: int, testable_patterns: list[str]) -> str:
    """
    For a given row (list of frequencies), finds the pattern with the
    maximum likelihood using the fast, direct-assignment method.
    
    frequencies: A list of [f1, f2, ... fn]
    testable_patterns: A list of all pattern strings, e.g., ["[4]", "[3, 1]", "[2, 2]", ...]
    """
    
    # --- 1. Handle Invalid Data ---
    if any(pd.isna(f) for f in frequencies):
        return "Invalid_Data"
        
    # --- 2. Get Counts ---
    f_sum = sum(frequencies)
    fD = max(0.0, 1.0 - f_sum)
    
    counts = [int(round(N_SCALE * f)) for f in frequencies]
    cD = int(round(N_SCALE * fD))
    
    # --- 3. Calculate LogL for all Patterns (Fast Method) ---
    pattern_logL = {pat: float("-inf") for pat in testable_patterns}

    # Iterate through each pattern (e.g., "[3, 1]")
    for pattern_key in pattern_logL.keys():
        
        # Convert string "[3, 1]" to list [3, 1]
        pattern = eval(pattern_key)
        
        # Create the single, direct-assignment genotype
        # e.g., [3, 1] for n=4 becomes [3, 1, 0, 0]
        genotype_counts = pattern.copy()
        genotype_counts.extend([0] * (n - len(genotype_counts)))
        
        # This is the ONLY genotype we test for this pattern
        # e.g., 3x Allele1, 1x Allele2, 0x Allele3, 0x Allele4
        
        # Get emission probabilities for this specific genotype
        (probs, pD) = get_emission_probs_list(genotype_counts, n, fD)
        
        # Calculate the log-likelihood
        logL = _log_likelihood_list(counts, cD, probs, pD)
            
        # Store this single score for the pattern
        pattern_logL[pattern_key] = logL

    # --- 4. Find the Best Pattern ---
    if all(logL == float("-inf") for logL in pattern_logL.values()):
        return "No_Fit"
        
    best_pattern = max(pattern_logL, key=pattern_logL.get)
    
    return best_pattern


def main():
    # --- 1. Load Data ---
    try:
        df = pd.read_csv(INPUT_CSV)
    except FileNotFoundError:
        print(f"Error: Input file not found at '{INPUT_CSV}'.", file=sys.stderr)
        return
    except Exception as e:
        print(f"An error occurred while reading the CSV file: {e}", file=sys.stderr)
        return

    # --- 2. Automatically Detect Ploidy (n) ---
    freq_cols = [col for col in df.columns if col.startswith('Frequency')]
    if not freq_cols:
        print(f"Error: No 'Frequency<n>' columns found in '{INPUT_CSV}'.", file=sys.stderr)
        return

    try:
        freq_numbers = [int(col.replace('Frequency', '')) for col in freq_cols]
        n_input = max(freq_numbers)
    except Exception as e:
        print(f"Error: Could not parse ploidy from column names: {freq_cols}", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
        return
        
    print(f"--- Automatically detected ploidy (n) = {n_input} ---")

    # --- 3. Generate All Patterns for this Ploidy ---
    all_patterns = find_all_partitions(n_input)
    all_pattern_keys = [str(p) for p in all_patterns]
    
    print(f"Generated {len(all_patterns)} total patterns to test.")
    print(f"Example patterns: {all_pattern_keys[:3]}...")

    # --- 4. Prepare Data for Calculation ---
    freq_col_names = [f'Frequency{i}' for i in range(1, n_input + 1)]
    
    if not all(col in df.columns for col in freq_col_names):
        print(f"Error: Input file is missing one or more required frequency columns.", file=sys.stderr)
        print(f"Expected: {freq_col_names}", file=sys.stderr)
        return
        
    print(f"Calculating best pattern for {len(df)} rows (fast mode)...")
    
    freq_num_col_names = []
    for i in range(1, n_input + 1):
        col_name = f'Frequency{i}'
        num_col_name = f'f_num_{i}'
        df[num_col_name] = pd.to_numeric(df[col_name], errors='coerce')
        freq_num_col_names.append(num_col_name)

    # --- 5. Apply Calculation ---
    df['Best_Pattern'] = df.apply(
        lambda row: compute_best_pattern(
            row[freq_num_col_names].tolist(),  # Pass the list of n frequencies
            n=n_input,
            testable_patterns=all_pattern_keys
        ),
        axis=1
    )
    
    # --- 6. Save Output ---
    output_cols_to_keep = []
    for col in ['Serial', 'Sample', 'Loci']:
        if col in df.columns:
            output_cols_to_keep.append(col)
            
    if not output_cols_to_keep:
        output_cols_to_keep = df.columns.tolist()[:3]
        
    output_cols_to_keep.append('Best_Pattern')
    output_df = df[output_cols_to_keep].copy()

    output_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSuccess! Wrote {len(output_df)} rows to '{OUTPUT_CSV}'.")


if __name__ == "__main__":
    main()