import csv
from scipy.stats import chisquare
import numpy as np

def analyze_single_row_patterns(data_row):
    """
    Performs Chi-Squared Goodness-of-Fit tests for a single row of data
    [total_sum_count, v1_count, v2_count] and returns the p-values
    for V1V1, V1V2, and V2V2 patterns.
    """
    _total_sum_count, v1_count, v2_count = data_row

    observed_counts = np.array([v1_count, v2_count])
    observed_total_for_test = v1_count + v2_count

    p_value_v1v1 = None
    p_value_v1v2 = None
    p_value_v2v2 = None

    if observed_total_for_test == 0:
        return [p_value_v1v1, p_value_v1v2, p_value_v2v2]

    # --- Scenario 1: Test for V1V1 Pattern (Expected: 100% V1, 0% V2) ---
    if v2_count > 0:
        p_value_v1v1 = 0.0 # Observed V2, so immediately reject pure V1V1
    else: # v2_count is 0, meaning all observed counts are V1
        # This is a perfect fit for V1V1
        p_value_v1v1 = 1.0

    # --- Scenario 2: Test for V2V2 Pattern (Expected: 0% V1, 100% V2) ---
    if v1_count > 0:
        p_value_v2v2 = 0.0 # Observed V1, so immediately reject pure V2V2
    else: # v1_count is 0, meaning all observed counts are V2
        # This is a perfect fit for V2V2
        p_value_v2v2 = 1.0

    # --- Scenario 3: Test for V1V2 Pattern (assuming 50/50 distribution) ---
    expected_counts_v1v2 = np.array([observed_total_for_test / 2, observed_total_for_test / 2])
    
    # Ensure expected counts are not zero, which would cause issues for chisquare
    # This scenario implies observed_total_for_test > 0, which is handled at function start.
    # If observed_total_for_test is odd (e.g. 1), expected will be 0.5, which chisquare handles.
    try:
        _, p_value_v1v2 = chisquare(f_obs=observed_counts, f_exp=expected_counts_v1v2)
    except ValueError:
        p_value_v1v2 = None # Indicate test couldn't be performed


    return [p_value_v1v1, p_value_v1v2, p_value_v2v2]

# Your provided file path
csv_file = "../new_output/AlleleInfo.csv"

print("Analysis Results from CSV:")
print("Row (Data: Total, V1, V2) -> [V1V1 P-value, V1V2 P-value, V2V2 P-value] -> Best Fit")
print("-" * 100)

try:
    with open(csv_file, newline='') as file:
        reader = csv.reader(file)
        header = next(reader)  # Skip header
        for i, row in enumerate(reader, 1):
            try:
                elements_str = row[5:8]
                # Convert to int, replacing 'N/A' with 0 (more statistically appropriate than 1)
                elements = [int(x) if x.strip().lower() != 'n/a' else 0 for x in elements_str]

                results = analyze_single_row_patterns(elements)
                
                formatted_results = [f"{p:.4f}" if p is not None else "N/A" for p in results]

                # Determine the 'best fit' among the patterns
                p_values = [p for p in results if p is not None] # Filter out None values
                
                best_fit_pattern = "None of the above (all rejected or N/A)"
                if p_values: # If there are any valid p-values
                    max_p_value = -1.0 # Initialize with a value lower than any possible p-value (p-values are [0, 1])
                    best_pattern_index = -1
                    
                    for idx, p in enumerate(results):
                        if p is not None and p > max_p_value:
                            max_p_value = p
                            best_pattern_index = idx

                    if best_pattern_index == 0:
                        best_fit_pattern = "V1V1"
                    elif best_pattern_index == 1:
                        best_fit_pattern = "V1V2"
                    elif best_pattern_index == 2:
                        best_fit_pattern = "V2V2"
                    
                    # If the best p-value is still very low, none is a 'good' fit
                    # Consider a threshold for what constitutes a 'good' fit (e.g., p > 0.05)
                    if max_p_value < 0.05: # Using a common alpha level for "poor fit"
                        best_fit_pattern = f"{best_fit_pattern} (poor fit, p={max_p_value:.4f})"
                    else: # If max_p_value is >= 0.05
                         best_fit_pattern = f"{best_fit_pattern} (good fit, p={max_p_value:.4f})"


                print(f"Row {i} ({elements}): [{formatted_results[0]}, {formatted_results[1]}, {formatted_results[2]}] -> {best_fit_pattern}")

            except IndexError:
                print(f"Warning: Row {i} has fewer than 8 columns. Skipping: {row}")
            except ValueError as e:
                print(f"Error processing row {i} ({row}): {e}. Check data format (expected integers).")
except FileNotFoundError:
    print(f"Error: CSV file not found at '{csv_file}'. Please check the path.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")