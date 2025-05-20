import csv
import pandas as pd
from collections import Counter

# Read the first file (Allele Info , also Table 1 ) into a matrix
file_path = '../new_output/AlleleInfo.csv'
df = pd.read_csv(file_path)

alleles = df[['Allele1', 'Allele2']]

#Read the second file (AlleleInformationFile_nSSR_FullLength_ParameterSet2_sa70_sb10_m10_n20.csv) into a matrix
file_path2 = '../AlleleInformationFile_nSSR_FullLength_ParameterSet2_sa70_sb10_m10_n20.csv'
df2 = pd.read_csv(file_path2, sep=" ")  # Specify space as the delimiter
print(df2)  

# Combine the 4th and 5th elements into a single long array
result_array = alleles.values.flatten()
 
# Drop NA elements from the long array
result_array = result_array[~pd.isna(result_array)]

# Drop all elements called "No data"
result_array = result_array[result_array != "No data"] #this results array should be used already

# Count the frequency of each element
frequency_count = Counter(result_array)

# Sort the frequency count from highest to lowest
sorted_frequency = sorted(frequency_count.items(), key=lambda x: x[1], reverse=True)

# Create a dictionary to store the first occurrence of each element
first_occurrence = {}
first_occurrence_Sample = {}

# matching_rows = []
# for element in result_array:
#     # Check if the element exists in the 5th column of df2
#     if element in seq_for_compare:
#         # Add the row(s) from df2 where the element exists
#         rows = df2[df2.iloc[:, 4] == element]
#         matching_rows.append(rows)

# print(matching_rows)

# Iterate through the sorted frequency list
for element, _ in sorted_frequency:
    # Find the first row where the element appears in either 'Allele1' or 'Allele2'
    first_row = df[(df['Allele1'] == element) | (df['Allele2'] == element)].iloc[0]
    
    # Store the first element (identification number) of the row
    first_occurrence[element] = first_row.iloc[0]  # Assuming the first column is the ID
    first_occurrence_Sample[element] = first_row.iloc[1]  # Assuming the second column is the sample name



# Save the first occurrence data along with frequency to a new CSV file
output_file_with_first_occurrence = '../new_output/sorted_frequency_with_first_occurrence.csv'
with open(output_file_with_first_occurrence, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['First Occurrence Serial', 'First Occurence Sample','Element', 'Frequency', ])  # Write header
    for element, freq in sorted_frequency:
        writer.writerow([first_occurrence[element], first_occurrence_Sample[element],element, freq ])  # Write data

print(f"Sorted frequency with first occurrence saved to {output_file_with_first_occurrence}")

# Save the sorted frequency as a CSV file
output_file = '../new_output/sorted_frequency.csv'
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Element', 'Frequency'])  # Write header
    writer.writerows(sorted_frequency)        # Write sorted data

print(f"Sorted frequency saved to {output_file}")
