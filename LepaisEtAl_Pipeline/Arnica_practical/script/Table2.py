import csv
import pandas as pd

# Read the first file into a matrix

file_path = '../new_output/AlleleInfo.csv'
df = pd.read_csv(file_path)


alleles = df[['Allele1', 'Allele2']]

# Filter out rows with NA values
alleles = alleles.dropna()

# Combine the 4th and 5th elements into a single long array
result_array = alleles.values.flatten()

# Print the result
print(result_array)

# def process_array_with_serial(array):
#     element_registry = {}  # Dictionary to store element serial numbers and counts
#     result = []  # List to store the processed results

#     for element in array:
#         if element not in element_registry:
#             # Register the element with a new serial number and set appear_number to 1
#             serial_number = len(element_registry) + 1
#             element_registry[element] = {"serial_number": serial_number, "appear_number": 1}
#         else:
#             # Increment the appear_number for the existing element
#             element_registry[element]["appear_number"] += 1

#         # Append the result for the current element
#         result.append({
#             "element": element,
#             "serial_number": element_registry[element]["serial_number"],
#             "appear_number": element_registry[element]["appear_number"]
#         })

#     return result

# # Example usage
# array = ["A", "B", "A", "C", "B", "A", "D", "C", "C"]
# processed_array = process_array_with_serial(array)

# # Print the results
# for entry in processed_array:
#     print(f"Element: {entry['element']}, Serial Number: {entry['serial_number']}, Appear Number: {entry['appear_number']}")
# for i in range(len(matrix1)):
#     seq = matrix1[i][3]
