import csv
import pandas as pd
from collections import Counter

# Read the first file into a matrix
file_path = '../new_output/AlleleInfo.csv'
df = pd.read_csv(file_path)

alleles = df[['Allele1', 'Allele2']]

# Filter out rows with NA values
alleles = alleles.dropna()

# Combine the 4th and 5th elements into a single long array
result_array = alleles.values.flatten()

# Count the frequency of each element
frequency_count = Counter(result_array)

# Sort the frequency count from highest to lowest
sorted_frequency = sorted(frequency_count.items(), key=lambda x: x[1], reverse=True)

# Save the sorted frequency as a CSV file
output_file = '../new_output/sorted_frequency.csv'
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Element', 'Frequency'])  # Write header
    writer.writerows(sorted_frequency)        # Write sorted data

print(f"Sorted frequency saved to {output_file}")
