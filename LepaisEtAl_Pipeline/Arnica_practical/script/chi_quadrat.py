import csv
from scipy.stats import chisquare
import numpy as np

def analyze_single_row_patterns(data_row):
    """
    Performs Chi-Squared Goodness-of-Fit tests for a single row of data
    [total_sum_count, v1_count, v2_count] and returns the p-values
    for V1V1, V1V2, and V2V2 patterns.
    """
    # Unpack the data row. total_sum_count is present but not used in
    # the chi-square calculations for V1/V2 proportions.
    # The sum v1_count + v2_count is used as the sample size for the test.
    _total_sum_count, v1_count, v2_count = data_row

    observed_counts = np.array([v1_count, v2_count])
    observed_total_for_test = v1_count + v2_count

    p_value_v1v1 = None
    p_value_v1v2 = None
    p_value_v2v2 = None

    if observed_total_for_test == 0:
        # If no V1 or V2 counts, cannot perform any test on these variants
        return [p_value_v1v1, p_value_v1v2, p_value_v2v2]

    # --- Scenario 1: Test for V1V1 Pattern (Expected: 100% V1, 0% V2) ---
    # If any V2 is observed, it fundamentally contradicts a pure V1V1 pattern.
    if v2_count > 0:
        p_value_v1v1 = 0.0  # Immediate rejection
    elif v1_count == 0: # If there's no V1, it can't be V1V1. Only happens if observed_total_for_test is non-zero (i.e. v2_count must be positive)
        p_value_v1v1 = 0.0
    else:
        # If v2_count is 0, then observed_counts is [v1_count, 0]
        # Expected counts are [observed_total_for_test, 0] which is [v1_count, 0]
        # This is a perfect match, chi2 = 0, p-value = 1.0
        expected_counts_v1v1 = np.array([float(observed_total_for_test), 0.0]) # Ensure floats for expected
        try:
            _, p_value_v1v1 = chisquare(f_obs=observed_counts, f_exp=expected_counts_v1v1)
        except ValueError:
            # This catch is mostly for safety, should ideally not be hit with the prior logic
            p_value_v1v1 = 0.0 # Default to rejection if an unexpected error occurs

    # --- Scenario 2: Test for V2V2 Pattern (Expected: 0% V1, 100% V2) ---
    # If any V1 is observed, it fundamentally contradicts a pure V2V2 pattern.
    if v1_count > 0:
        p_value_v2v2 = 0.0  # Immediate rejection
    elif v2_count == 0: # If there's no V2, it can't be V2V2. Only happens if observed_total_for_test is non-zero (i.e. v1_count must be positive)
        p_value_v2v2 = 0.0
    else:
        # If v1_count is 0, then observed_counts is [0, v2_count]
        # Expected counts are [0, observed_total_for_test] which is [0, v2_count]
        # This is a perfect match, chi2 = 0, p-value = 1.0
        expected_counts_v2v2 = np.array([0.0, float(observed_total_for_test)]) # Ensure floats
        try:
            _, p_value_v2v2 = chisquare(f_obs=observed_counts, f_exp=expected_counts_v2v2)
        except ValueError:
            p_value_v2v2 = 0.0

    # --- Scenario 3: Test for V1V2 Pattern (assuming 50/50 distribution) ---
    # This test requires both expected counts to be non-zero for a valid chisquare
    # application that includes both categories.
    expected_counts_v1v2 = np.array([observed_total_for_test / 2, observed_total_for_test / 2])

    # If observed_total_for_test is 0, this is handled at the start.
    # If it's 1, then expected counts are [0.5, 0.5], which are valid but small.
    # chisquare handles small expected values, but a warning about accuracy might be issued by scipy internally.
    try:
        _, p_value_v1v2 = chisquare(f_obs=observed_counts, f_exp=expected_counts_v1v2)
    except ValueError:
        # This might happen if observed_total_for_test is 0 (already handled),
        # or if expected_counts become NaN for some unexpected reason.
        p_value_v1v2 = None # Indicate test couldn't be performed

    return [p_value_v1v1, p_value_v1v2, p_value_v2v2]

# Your provided file path
csv_file = "../new_output/AlleleInfo.csv"

# Example data for testing without a file (if you want to run standalone)
# data_for_testing = [
#     [121, 85, 1],
#     [4599, 4171, 1],
#     [3635, 1513, 1283],
#     [10, 0, 0], # total 10, v1=0, v2=0 (meaning other variants contributed to total_sum_count)
#     [5, 5, 0], # Pure V1
#     [7, 0, 7], # Pure V2
#     [20, 10, 10], # Pure V1V2
#     [100, 99, 1], # V1V1 very likely, but with one V2
#     [100, 1, 99], # V2V2 very likely, but with one V1
#     [2, 1, 0], # Small N, one V1, no V2
#     [2, 0, 1], # Small N, no V1, one V2
#     [2, 1, 1], # Small N, V1V2
#     [0, 0, 0], # No counts at all
# ]
# for i, row in enumerate(data_for_testing):
#     results = analyze_single_row_patterns(row)
#     formatted_results = [f"{p:.4f}" if p is not None else "N/A" for p in results]
#     print(f"Row {i+1} ({row}): [V1V1 P-value: {formatted_results[0]}, V1V2 P-value: {formatted_results[1]}, V2V2 P-value: {formatted_results[2]}]")


# Actual CSV file processing
try:
    with open(csv_file, newline='') as file:
        reader = csv.reader(file)
        header = next(reader)  # Skip header if present
        print("Analysis Results from CSV:")
        print("Row (Original Data): [V1V1 P-value, V1V2 P-value, V2V2 P-value]")
        print("-" * 70)
        for i, row in enumerate(reader, 1):
            try:
                # Extract the 6th to 8th elements (index 5, 6, 7)
                elements_str = row[5:8]
                # Replace 'N/A' with '1' and convert to int
                # IMPORTANT: Replacing 'N/A' with 1 might not be statistically sound
                # if 'N/A' truly means missing data or 0.
                # Consider replacing with 0 or handling as missing.
                # For this cross-check, I'll stick to your logic, but flag it.
                elements = [int(x) if x.strip().lower() != 'n/a' else 1 for x in elements_str]

                results = analyze_single_row_patterns(elements)
                formatted_results = [f"{p:.4f}" if p is not None else "N/A" for p in results]
                print(f"Row {i} ({elements}): [V1V1 P-value: {formatted_results[0]}, V1V2 P-value: {formatted_results[1]}, V2V2 P-value: {formatted_results[2]}]")
            except IndexError:
                print(f"Warning: Row {i} has fewer than 8 columns. Skipping: {row}")
            except ValueError as e:
                print(f"Error processing row {i} ({row}): {e}. Check data format (expected integers).")
except FileNotFoundError:
    print(f"Error: CSV file not found at '{csv_file}'. Please check the path.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")