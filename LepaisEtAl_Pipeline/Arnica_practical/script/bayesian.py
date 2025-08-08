import csv
from scipy.stats import chisquare
import numpy as np

# (The calculate_pattern_relative_possibility_from_pvalues function remains the same as before)
def calculate_pattern_relative_possibility_from_pvalues(data_row):
    _total_sum_count, v1_count, v2_count = data_row
    observed_counts = np.array([v1_count, v2_count])
    observed_total_for_test = v1_count + v2_count
    p_value_v1v1 = None
    p_value_v1v2 = None
    p_value_v2v2 = None

    if observed_total_for_test == 0:
        return [None, None, None], [None, None, None] # Return p-values and possibilities

    if v2_count > 0:
        p_value_v1v1 = 0.0
    else:
        p_value_v1v1 = 1.0

    if v1_count > 0:
        p_value_v2v2 = 0.0
    else:
        p_value_v2v2 = 1.0

    expected_counts_v1v2 = np.array([observed_total_for_test / 2, observed_total_for_test / 2])
    try:
        _, p_value_v1v2 = chisquare(f_obs=observed_counts, f_exp=expected_counts_v1v2)
    except ValueError:
        p_value_v1v2 = None 

    raw_p_values = [p_value_v1v1, p_value_v1v2, p_value_v2v2]
    
    possibility_scores = []
    for p_val in raw_p_values:
        if p_val is None:
            possibility_scores.append(0.0) 
        elif p_val == 0.0:
            possibility_scores.append(1e-300) 
        else:
            possibility_scores.append(p_val)

    total_score = sum(possibility_scores)
    if total_score == 0:
        normalized_possibilities = [0.0, 0.0, 0.0]
    else:
        normalized_possibilities = [score / total_score for score in possibility_scores]

    return normalized_possibilities, raw_p_values # Return both normalized scores and raw p-values


# Your provided file path
csv_file = "../new_output/AlleleInfo.csv"

print("Analysis Results from CSV (Relative Possibilities based on p-values):")
print("Row (Data: Total, V1, V2) -> [P(V1V1), P(V1V2), P(V2V2)] -> Best Fit (Absolute Fit Status, max_P=Normalized_Score)")
print("-" * 120)

try:
    with open(csv_file, newline='') as file:
        reader = csv.reader(file)
        header = next(reader)  # Skip header
        for i, row in enumerate(reader, 1):
            try:
                elements_str = row[5:8]
                elements = [int(x) if x.strip().lower() != 'n/a' else 0 for x in elements_str]

                # Get both normalized possibilities and raw p-values
                possibilities, raw_p_values = calculate_pattern_relative_possibility_from_pvalues(elements)
                
                formatted_possibilities = [f"{p:.4f}" if p is not None else "N/A" for p in possibilities]

                best_fit_pattern_name = "N/A" # Name of the pattern (V1V1, V1V2, V2V2)
                best_fit_status = "" # Good fit, poor fit, all impossible etc.
                best_normalized_score = 0.0 # Max of the normalized possibilities
                best_raw_p_value = -1.0 # P-value of the best fitting model (for absolute assessment)

                # Determine the best fit pattern based on normalized possibilities
                if possibilities and all(p is not None for p in possibilities):
                    max_possibility = -1.0
                    best_pattern_index = -1
                    
                    for idx, p in enumerate(possibilities):
                        if p > max_possibility:
                            max_possibility = p
                            best_pattern_index = idx

                    best_normalized_score = max_possibility

                    if best_pattern_index == 0:
                        best_fit_pattern_name = "V1V1"
                    elif best_pattern_index == 1:
                        best_fit_pattern_name = "V1V2"
                    elif best_pattern_index == 2:
                        best_fit_pattern_name = "V2V2"
                    
                    best_raw_p_value = raw_p_values[best_pattern_index]

                    # Now, assess the *absolute* quality of the best fit using its raw p-value
                    if best_raw_p_value is None:
                        best_fit_status = "Test not applicable"
                    elif best_raw_p_value >= 0.05: # Conventional threshold for significance
                        best_fit_status = "Good fit"
                    else: # best_raw_p_value < 0.05
                        best_fit_status = "Poor fit" # Statistically significant deviation

                    # Edge case: if all normalized probabilities are 0 (e.g. from all original p-values being None)
                    if np.isclose(best_normalized_score, 0.0) and best_raw_p_value is None:
                         best_fit_pattern_name = "N/A"
                         best_fit_status = "No data for test"
                    elif np.isclose(best_normalized_score, 0.0) and best_raw_p_value == 0.0:
                         best_fit_pattern_name = "All" # All patterns
                         best_fit_status = "Extremely poor fit (all rejected)"


                print(f"Row {i} ({elements}): [{formatted_possibilities[0]}, {formatted_possibilities[1]}, {formatted_possibilities[2]}] -> {best_fit_pattern_name} ({best_fit_status}, max_P={best_normalized_score:.4f}, raw_p={best_raw_p_value:.4f})")

            except IndexError:
                print(f"Warning: Row {i} has fewer than 8 columns. Skipping: {row}")
            except ValueError as e:
                print(f"Error processing row {i} ({row}): {e}. Check data format (expected integers).")
except FileNotFoundError:
    print(f"Error: CSV file not found at '{csv_file}'. Please check the path.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")