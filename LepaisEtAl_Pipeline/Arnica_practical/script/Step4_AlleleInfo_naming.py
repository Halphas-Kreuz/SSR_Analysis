import csv

def load_dictionary(dict_path):
    dictionary = []
    with open(dict_path, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            if len(row) >= 3:
                dictionary.append(row)
    return dictionary

def search_in_dictionary(dictionary, search_string):
    for row in dictionary:
        if row[3] == search_string:
            return row[1]
    return search_string  # If not found, keep original

alleleinfo_path = '../new_output/AlleleInfo.csv'
dict_path = '../new_output/AlleleInfo_dictionary.csv'
output_path = '../new_output/AlleleInfo_named.csv'

# Load dictionary
dictionary = load_dictionary(dict_path)

# Read, convert, and write to a new file
with open(alleleinfo_path, 'r') as infile, open(output_path, 'w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)
    header = next(reader)
    writer.writerow(header)
    for row in reader:
        # Replace 4th and 5th columns (index 3 and 4) with names from dictionary
        row[3] = search_in_dictionary(dictionary, row[3])
        row[4] = search_in_dictionary(dictionary, row[4])
        writer.writerow(row)

print(f"AlleleInfo with names saved to {output_path}")