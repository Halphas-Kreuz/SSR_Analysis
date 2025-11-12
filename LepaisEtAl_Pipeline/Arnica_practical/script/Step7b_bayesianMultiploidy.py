#!/usr/bin/env python3
"""
Bayesian Genotype Likelihood Calculator (Diploid, Custom Errors) - v_final

This script calculates the likelihood of AA vs AB genotypes based on
observed read fractions and a locus-specific, asymmetric error map.

It implements the following logic for invalid data:
1.  If Allele1 or Frequency1 is invalid, the row is invalid.
    (like_AA = 0, like_AB = 0)
2.  If Allele1 is valid, but Allele2 or Frequency2 is invalid,
    this is a confident AA genotype.
    (like_AA = 1, like_AB = 0)
3.  If all data is valid, proceed with the Bayesian calculation.
"""

from __future__ import annotations
import math
from datetime import datetime
import pandas as pd
import numpy as np
import sys
from collections import defaultdict

# --- Parameters ---
N_SCALE = 1000  # Pseudo-count scaling factor
INPUT_CSV = "../new_output/AlleleInfo_ready.csv"
ERROR_CSV = "../new_output/table/All_SequenceComparisons.csv"
OUTPUT_CSV = f"../new_output/likelihoods_AA_AB_{datetime.now():%Y%m%d_%H%M%S}.csv"




def build_error_map(error_file_path: str) -> dict:
    """
    Loads the error CSV and builds a fast, nested lookup map.
    Structure: map[Locus][(True_Allele, Read_Allele)] = Probability
    """
    print(f"Building error map from '{error_file_path}'...")
    error_map = defaultdict(dict)
    try:
        df_err = pd.read_csv(error_file_path)
        
        required_cols = ['Locus', 'Allele1', 'Allele2', 'Value']
        if not all(col in df_err.columns for col in required_cols):
            print(f"Error: Error file must contain {required_cols}", file=sys.stderr)
            return None

        for row in df_err.itertuples(index=False):
            locus = row.Locus
            true_allele = row.Allele1
            read_allele = row.Allele2
            prob = row.Value
            error_map[locus][(true_allele, read_allele)] = prob
            
        print(f"Successfully built error map for {len(error_map)} loci.")
        return error_map

    except FileNotFoundError:
        print(f"Error: Error file not found at '{error_file_path}'", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error reading error file: {e}", file=sys.stderr)
        return None


def _validate_emissions_nonnegative(pA, pB, pD):
    if pA < -1e-12 or pB < -1e-12 or pD < -1e-12:
        raise ValueError(f"Negative emission probability encountered: {(pA, pB, pD)}")
    s = pA + pB + pD
    if not abs(s - 1.0) < 1e-9:
        raise ValueError(f"Emission probabilities do not sum to 1: {s}")


def _emissions(eA2B: float, eB2A: float, dA: float, dB: float):
    """
    Return emission triplets (pA, pB, pD) for AA, AB, BB in that order.
    """
    # True A (genotype AA)
    pA_AA = 1.0 - eA2B - dA
    pB_AA = eA2B
    pD_AA = dA

    # True B (genotype BB) - Needed for AB calculation
    pB_BB = 1.0 - eB2A - dB
    pA_BB = eB2A
    pD_BB = dB

    for trip in ((pA_AA, pB_AA, pD_AA), (pA_BB, pB_BB, pD_BB)):
        if any(v < -1e-12 for v in trip):
            raise ValueError(f"Invalid parameters: e + d > 1. eA2B={eA2B}, eB2A={eB2A}, dA={dA}, dB={dB}")

    def _clip(a, b, d):
        a = max(0.0, a); b = max(0.0, b); d = max(0.0, d)
        s = a + b + d
        if s == 0.0: return (1/3, 1/3, 1.3)
        return (a / s, b / s, d / s)

    pA_AA, pB_AA, pD_AA = _clip(pA_AA, pB_AA, pD_AA)
    pA_BB, pB_BB, pD_BB = _clip(pA_BB, pB_BB, pD_BB)

    # AB is a 50/50 mixture
    pA_AB = 0.5 * (pA_AA + pA_BB)
    pB_AB = 0.5 * (pB_AA + pB_BB)
    pD_AB = 0.5 * (pD_AA + pD_BB)

    for trip in ((pA_AA, pB_AA, pD_AA), (pA_AB, pB_AB, pD_AB), (pA_BB, pB_BB, pD_BB)):
        _validate_emissions_nonnegative(*trip)

    return (pA_AA, pB_AA, pD_AA), (pA_AB, pB_AB, pD_AB)


def _log_term(count: int, prob: float) -> float:
    if prob <= 0.0:
        return 0.0 if count == 0 else float("-inf")
    return count * math.log(prob)


def _log_likelihood(cA: int, cB: int, cD: int, pA: float, pB: float, pD: float) -> float:
    terms = (
        _log_term(cA, pA),
        _log_term(cB, pB),
        _log_term(cD, pD),
    )
    if any(t == float("-inf") for t in terms):
        if (cA > 0 and pA == 0.0) or (cB > 0 and pB == 0.0) or (cD > 0 and pD == 0.0):
            return float("-inf")
    return sum(terms)


def _softmax2(loga: float, logb: float):
    m = max(loga, logb)
    if m == float("-inf"):
        return (0.5, 0.5)
    ea = math.exp(loga - m) if loga != float("-inf") else 0.0
    eb = math.exp(logb - m) if logb != float("-inf") else 0.0
    s = ea + eb
    if s == 0.0:
        return (0.5, 0.5)
    return ea / s, eb / s


def compute_row(row: pd.Series, error_map: dict, N: int = N_SCALE) -> pd.Series:
    """
    Computes the AA and AB likelihoods using the new invalid data rules.
    """
    # 1. Get data from the row
    locus = row.get('Loci')
    seq_i = str(row.get('Allele1')) # 'A'
    seq_j = str(row.get('Allele2')) # 'B'
    fA = row.get('fA_num')
    fB = row.get('fB_num')

    # 2. Implement the new Invalid Data Logic
    
    # Rule 1: Fatal Invalidity (Allele1 or Freq1 is invalid)
    if pd.isna(fA) or seq_i in ("No data", "N/A","Other sequences"):
        return pd.Series({"like_AA": 0, "like_AB": 0})

    # Rule 2: Pure AA (Allele2 or Freq2 is invalid, but A is valid)
    if pd.isna(fB) or seq_j in ("No data", "N/A","Other sequences"):
        return pd.Series({"like_AA": 1, "like_AB": 0})

    # Rule 3: Valid Data. Proceed with Bayesian calculation.
    
    # Check if locus is in our error map
    if locus not in error_map:
        eA2B = 0.0
        eB2A = 0.0
    else:
        locus_errors = error_map[locus]
        eA2B = locus_errors.get((seq_i, seq_j), 0.0) # P(read B | true A)
        eB2A = locus_errors.get((seq_j, seq_i), 0.0) # P(read A | true B)

    # 4. Calculate counts and discard
    fD = max(0.0, 1.0 - fA - fB)
    cA = int(round(N * fA))
    cB = int(round(N * fB))
    cD = int(round(N * fD))

    # 5. Get Emission Models
    try:
        (pA_AA, pB_AA, pD_AA), (pA_AB, pB_AB, pD_AB) = _emissions(
            eA2B=eA2B, 
            eB2A=eB2A, 
            dA=fD, 
            dB=fD
        )
    except ValueError:
        # Error (e.g., e + d > 1) means params are impossible
        return pd.Series({"like_AA": 0.5, "like_AB": 0.5}) # Fallback uncertainty

    # 6. Log-likelihoods (Only for AA and AB)
    logL_AA = _log_likelihood(cA, cB, cD, pA_AA, pB_AA, pD_AA)
    logL_AB = _log_likelihood(cA, cB, cD, pA_AB, pB_AB, pD_AB)

    # 7. Normalize (Softmax on 2 variables)
    like_AA, like_AB = _softmax2(logL_AA, logL_AB)

    return pd.Series({"like_AA": like_AA, "like_AB": like_AB})


def main():
    # 1. Build the error map
    error_map = build_error_map(ERROR_CSV)
    if error_map is None:
        print("Aborting due to error map failure.")
        return

    # 2. Load the input data
    try:
        df = pd.read_csv(INPUT_CSV)
    except FileNotFoundError:
        print(f"Error: Input file not found at '{INPUT_CSV}'.", file=sys.stderr)
        return
    except Exception as e:
        print(f"An error occurred while reading the input CSV: {e}", file=sys.stderr)
        return

    # 3. Check for required columns
    required_cols = ['Loci', 'Sample', 'Allele1', 'Allele2', 'Frequency1', 'Frequency2']
    if not all(col in df.columns for col in required_cols):
        print(f"Error: The input file '{INPUT_CSV}' must contain columns: {required_cols}", file=sys.stderr)
        return
    
    print(f"Calculating likelihoods for {len(df)} rows...")
    
    # 4. Coerce frequencies to numeric type
    df['fA_num'] = pd.to_numeric(df['Frequency1'], errors='coerce')
    df['fB_num'] = pd.to_numeric(df['Frequency2'], errors='coerce')

    # 5. Apply the compute_row function to each row
    results_df = df.apply(
        lambda row: compute_row(
            row,
            error_map=error_map,
            N=N_SCALE
        ), 
        axis=1
    )
    
  # 6. Prepare final output DataFrame
    # Added 'Serial' to this list
    output_cols = ['Serial', 'Loci', 'Sample', 'Allele1', 'Allele2']
    final_df = df[output_cols].copy()
    
    # Add the new calculated likelihoods
    final_df['like_AA'] = results_df['like_AA']
    final_df['like_AB'] = results_df['like_AB']
    
    # 7. Convert 0.0/1.0 to 0/1 for cleaner output
    final_df['like_AA'] = final_df['like_AA'].astype(object)
    final_df['like_AB'] = final_df['like_AB'].astype(object)
    final_df.loc[final_df['like_AA'] == 1.0, 'like_AA'] = 1
    final_df.loc[final_df['like_AA'] == 0.0, 'like_AA'] = 0
    final_df.loc[final_df['like_AB'] == 1.0, 'like_AB'] = 1
    final_df.loc[final_df['like_AB'] == 0.0, 'like_AB'] = 0

    # 8. Save the final DataFrame to a new CSV file
    final_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSuccess! Wrote {len(final_df)} rows to '{OUTPUT_CSV}'.")

if __name__ == "__main__":
    main()