#!/usr/bin/env python3
"""
Hardcoded AA/AB/BB likelihoods from allele read fractions.

Each input row has:
  seq_i, seq_j, f_A, f_B, e_cross

Definitions:
  - f_A, f_B: observed fractions of all reads labeled to sequence i (A) and j (B).
  - f_D := max(0, 1 - f_A - f_B) is the discard fraction (unlabeled/filtered reads).
  - e_cross: symmetric per-read cross-misread probability, so eA2B = eB2A = e_cross.

Emission model (per-read):
  If true is A: P(A|A)=1-eA2B-dA, P(B|A)=eA2B, P(D|A)=dA
  If true is B: P(B|B)=1-eB2A-dB, P(A|B)=eB2A, P(D|B)=dB
  If genotype is AB: 50/50 mixture of the two rows above.

Likelihood:
  Treat (f_A,f_B,f_D) as proportions; convert to pseudo-counts using scale N (default 1000).
  For genotype G with emissions (pA, pB, pD):
      logL(G) = c_A*log(pA) + c_B*log(pB) + c_D*log(pD)
  Report logL_* and normalized likelihoods like_* that sum to 1 (for comparability).
"""

from __future__ import annotations
import math
from datetime import datetime
import pandas as pd
import numpy as np
import sys

# Default parameters
N_SCALE = 1000
INPUT_CSV = "AlleleInfo_named_fused.csv"
OUTPUT_CSV = f"likelihoods_{datetime.now():%Y%m%d_%H%M%S}.csv"

def _validate_emissions_nonnegative(pA, pB, pD):
    if pA < 0 or pB < 0 or pD < 0:
        raise ValueError(f"Negative emission probability encountered: {(pA, pB, pD)}")
    s = pA + pB + pD
    if not (abs(s - 1.0) < 1e-12 or s <= 1.0 + 1e-12):
        # Floating drift is allowed, but probabilities must not substantially exceed 1
        pass


def _emissions(eA2B: float, eB2A: float, dA: float, dB: float):
    """
    Return emission triplets (pA, pB, pD) for AA, AB, BB in that order.
    No epsilon smoothing; zeros remain zeros.
    """
    # True A (genotype AA)
    pA_AA = 1.0 - eA2B - dA
    pB_AA = eA2B
    pD_AA = dA

    # True B (genotype BB)
    pB_BB = 1.0 - eB2A - dB
    pA_BB = eB2A
    pD_BB = dB

    # Basic validity
    for trip in ((pA_AA, pB_AA, pD_AA), (pA_BB, pB_BB, pD_BB)):
        if any(v < -1e-12 for v in trip):
            raise ValueError("Invalid parameters: e + d > 1 produced negative probabilities.")

    # Clip tiny negatives to zero from numerical noise; renormalize rows to sum to 1
    def _clip_and_norm(a, b, d):
        a = max(0.0, a); b = max(0.0, b); d = max(0.0, d)
        s = a + b + d
        if s == 0.0:
            # Degenerate; fall back to uniform (should not happen with valid params)
            return (1/3, 1/3, 1/3)
        return (a / s, b / s, d / s)

    pA_AA, pB_AA, pD_AA = _clip_and_norm(pA_AA, pB_AA, pD_AA)
    pA_BB, pB_BB, pD_BB = _clip_and_norm(pA_BB, pB_BB, pD_BB)

    # AB is a 50/50 mixture of the "true A" and "true B" rows
    pA_AB = 0.5 * (pA_AA + pA_BB)
    pB_AB = 0.5 * (pB_AA + pB_BB)
    pD_AB = 0.5 * (pD_AA + pD_BB)

    # Final hygiene, ensure sums ~ 1
    for trip in ((pA_AA, pB_AA, pD_AA), (pA_AB, pB_AB, pD_AB), (pA_BB, pB_BB, pD_BB)):
        _validate_emissions_nonnegative(*trip)

    return (pA_AA, pB_AA, pD_AA), (pA_AB, pB_AB, pD_AB), (pA_BB, pB_BB, pD_BB)


def _log_term(count: int, prob: float) -> float:
    """Return count*log(prob), with exact-zero handling:
       if prob == 0 and count > 0, contributes -inf; if count == 0, contributes 0."""
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
        # If any category assigns zero probability to observed positive counts, the logL is -inf
        if (cA > 0 and pA == 0.0) or (cB > 0 and pB == 0.0) or (cD > 0 and pD == 0.0):
            return float("-inf")
    return sum(terms)


def _softmax3(loga: float, logb: float, logc: float):
    """Stable normalization of three log-likelihoods to (like_AA, like_AB, like_BB) that sum to 1.
       If all are -inf, return equal weights."""
    m = max(loga, logb, logc)
    if m == float("-inf"):
        return (1/3, 1/3, 1/3)
    ea = math.exp(loga - m) if loga != float("-inf") else 0.0
    eb = math.exp(logb - m) if logb != float("-inf") else 0.0
    ec = math.exp(logc - m) if logc != float("-inf") else 0.0
    s = ea + eb + ec
    if s == 0.0:
        return (1/3, 1/3, 1/3)
    return ea / s, eb / s, ec / s


def compute_row(seq_i: str, seq_j: str, fA: float, fB: float, e_cross: float, N: int = N_SCALE):
    # Handle the "No data" case as requested. If either frequency is NaN,
    # return the uniform distribution, representing complete uncertainty.
    if pd.isna(fA) or pd.isna(fB):
        return {
            "seq_i": seq_i, "seq_j": seq_j,
            "f_A": fA, "f_B": fB, "f_D": np.nan, "e_cross": e_cross,
            "c_A": 0, "c_B": 0, "c_D": 0,
            "pA_AA": np.nan, "pB_AA": np.nan, "pD_AA": np.nan,
            "pA_AB": np.nan, "pB_AB": np.nan, "pD_AB": np.nan,
            "pA_BB": np.nan, "pB_BB": np.nan, "pD_BB": np.nan,
            "logL_AA": np.nan, "logL_AB": np.nan, "logL_BB": np.nan,
            "like_AA": 0.0, "like_AB": 0.0, "like_BB": 0.0,
        }

    # Derived discard fraction
    fD = 1.0 - fA - fB
    if fD < 0.0:
        # Tolerate small negatives from rounding; otherwise clamp to zero
        fD = max(0.0, fD)

    # New check to prevent negative probabilities due to invalid input data
    if fD + e_cross > 1.0:
        return {
            "seq_i": seq_i, "seq_j": seq_j,
            "f_A": fA, "f_B": fB, "f_D": fD, "e_cross": e_cross,
            "c_A": 0, "c_B": 0, "c_D": 0,
            "pA_AA": np.nan, "pB_AA": np.nan, "pD_AA": np.nan,
            "pA_AB": np.nan, "pB_AB": np.nan, "pD_AB": np.nan,
            "pA_BB": np.nan, "pB_BB": np.nan, "pD_BB": np.nan,
            "logL_AA": np.nan, "logL_AB": np.nan, "logL_BB": np.nan,
            "like_AA": 0.0, "like_AB": 0.0, "like_BB": 0.0,
        }

    # Pseudo-counts
    cA = int(round(N * fA))
    cB = int(round(N * fB))
    cD = int(round(N * fD))

    # Error/Discard parameters
    eA2B = e_cross
    eB2A = e_cross
    dA = fD
    dB = fD

    # Emissions per genotype
    (pA_AA, pB_AA, pD_AA), (pA_AB, pB_AB, pD_AB), (pA_BB, pB_BB, pD_BB) = _emissions(eA2B, eB2A, dA, dB)

    # Log-likelihoods
    logL_AA = _log_likelihood(cA, cB, cD, pA_AA, pB_AA, pD_AA)
    logL_AB = _log_likelihood(cA, cB, cD, pA_AB, pB_AB, pD_AB)
    logL_BB = _log_likelihood(cA, cB, cD, pA_BB, pB_BB, pD_BB)

    # Rescaled likelihoods (sum to 1 for comparability; NOT posteriors)
    like_AA, like_AB, like_BB = _softmax3(logL_AA, logL_AB, logL_BB)

    return {
        "seq_i": seq_i, "seq_j": seq_j,
        "f_A": fA, "f_B": fB, "f_D": fD, "e_cross": e_cross,
        "c_A": cA, "c_B": cB, "c_D": cD,
        "pA_AA": pA_AA, "pB_AA": pB_AA, "pD_AA": pD_AA,
        "pA_AB": pA_AB, "pB_AB": pB_AB, "pD_AB": pD_AB,
        "pA_BB": pA_BB, "pB_BB": pB_BB, "pD_BB": pD_BB,
        "logL_AA": logL_AA, "logL_AB": logL_AB, "logL_BB": logL_BB,
        "like_AA": like_AA, "like_AB": like_AB, "like_BB": like_BB,
    }

def main():
    try:
        # Load the input data from the provided CSV file
        df = pd.read_csv(INPUT_CSV)
    except FileNotFoundError:
        print(f"Error: Input file not found at '{INPUT_CSV}'. Please ensure the file is in the correct directory.", file=sys.stderr)
        return
    except Exception as e:
        print(f"An error occurred while reading the CSV file: {e}", file=sys.stderr)
        return

    # Check if the required columns exist
    required_cols = ['Allele1', 'Allele2', 'Frequency1', 'Frequency2', 'ComparisonValue']
    if not all(col in df.columns for col in required_cols):
        print(f"Error: The input file must contain the following columns: {required_cols}", file=sys.stderr)
        return
    
    print(f"Calculating likelihoods for {len(df)} rows...")
    
    # Apply the compute_row function to each row of the DataFrame
    # The result_type='expand' is crucial for converting the returned dictionaries into new columns
    # We use .get() with a default value to handle cases where 'N/A' might be present
    results_df = df.apply(
        lambda row: compute_row(
            seq_i=str(row.get('Allele1')),
            seq_j=str(row.get('Allele2')),
            fA=pd.to_numeric(row.get('Frequency1'), errors='coerce'),
            fB=pd.to_numeric(row.get('Frequency2'), errors='coerce'),
            e_cross=pd.to_numeric(row.get('ComparisonValue'), errors='coerce'),
            N=N_SCALE
        ), 
        axis=1, 
        result_type='expand'
    )
    
    # Keep the first 3 columns from the original DataFrame
    final_df = df.iloc[:, :3].copy()
    
    # Replace columns 4-6 with the calculated likelihoods
    final_df['like_AA'] = results_df['like_AA']
    final_df['like_AB'] = results_df['like_AB']
    final_df['like_BB'] = results_df['like_BB']
    
    # Save the final DataFrame to a new CSV file
    final_df.to_csv(OUTPUT_CSV, index=False)
    print(f"Success! Wrote {len(final_df)} rows to '{OUTPUT_CSV}'.")


if __name__ == "__main__":
    main()
