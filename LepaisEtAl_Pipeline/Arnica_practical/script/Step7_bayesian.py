#!/usr/bin/env python3
"""
bayes_genotype_posterior.py

Posterior formula (given "different" observation D):
  P(homo | D) = (e * h) / (e * h + (1 - e) * (1 - h))

Optionally, you can output P(hetero | D) instead, which is simply 1 - P(homo | D).


Angela wrote this to be robust to a label column in the first position.
"""

import argparse
import pandas as pd
import numpy as np
import os


# --- Hardcoded settings ---
INPUT_DIR = "../new_output/table/"
OUTPUT_DIR = "../new_output/comparison/"
ERROR = 0.1
MIN_PROB = 0.01
MAX_PROB = 0.99
MODE = "homo"  # or "hetero"
NO_INDEX_HEURISTIC = False
# --------------------------


def parse_args():
    p = argparse.ArgumentParser(description="Bayesian posterior for genotype calls from prior homozygosity.")
    p.add_argument("--input", "-i", required=True, help="Input CSV: cells are prior P(homozygous).")
    p.add_argument("--output", "-o", required=True, help="Output CSV for posterior probabilities.")
    p.add_argument("--error", "-e", type=float, default=0.1,
                   help="Symmetric error rate e (default: 0.1).")
    p.add_argument("--min-prob", type=float, default=0.01,
                   help="Softening for exact zeros: replace 0 with this (default: 0.01).")
    p.add_argument("--max-prob", type=float, default=0.99,
                   help="Softening for exact ones: replace 1 with this (default: 0.99).")
    p.add_argument("--mode", choices=["homo", "hetero"], default="homo",
                   help="Which posterior to output: homo or hetero (default: homo).")
    p.add_argument("--no-index-heuristic", action="store_true",
                   help="Disable first-column-as-index heuristic.")
    return p.parse_args()


def maybe_set_index(df: pd.DataFrame, disable: bool) -> pd.DataFrame:
    if disable:
        return df
    if df.shape[1] == 0:
        return df
    first_col = df.columns[0]
    # Heuristic: if the first column is non-numeric OR has non-unique values, treat as labels.
    if (not pd.api.types.is_numeric_dtype(df[first_col])) or (not df[first_col].is_unique):
        return df.set_index(first_col)
    return df


def process_file(input_path, output_path):
    df = pd.read_csv(input_path)
    df = maybe_set_index(df, disable=NO_INDEX_HEURISTIC)
    h = df.apply(pd.to_numeric, errors="coerce").clip(lower=0.0, upper=1.0)
    h = h.replace(0.0, MIN_PROB).replace(1.0, MAX_PROB)
    e = ERROR
    P_D_given_H = 1.0 - e
    P_D_given_notH = e
    den = (P_D_given_notH * h) + (P_D_given_H * (1.0 - h))
    den = den.replace(0, np.nan)
    posterior_homo = (P_D_given_notH * h) / den
    if MODE == "homo":
        out = posterior_homo
    else:
        out = 1.0 - posterior_homo

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    out.to_csv(output_path, index=True)


def main():
    for root, dirs, files in os.walk(INPUT_DIR):
        for file in files:
            if file.endswith(".csv"):
                input_path = os.path.join(root, file)
                # Get relative path from INPUT_DIR
                rel_path = os.path.relpath(input_path, INPUT_DIR)
                # Remove .csv extension and add _posterior.csv
                rel_base = os.path.splitext(rel_path)[0] + "_posterior.csv"
                output_path = os.path.join(OUTPUT_DIR, rel_base)
                print(f"Processing {input_path} -> {output_path}")
                process_file(input_path, output_path)


if __name__ == "__main__":
    main()
