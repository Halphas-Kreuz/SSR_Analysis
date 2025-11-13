import csv
import collections

def extract_allele_name(locus, raw_allele_string):
    """
    Cleans the raw allele string, removing the locus prefix
    or converting non-data to 'NA'.
    """
    # Handle missing data or non-allele entries
    if raw_allele_string in ["No data", "Other sequences", "", None]:
        return "NA"
    
    # Check if the string starts with the locus name and an underscore
    prefix = locus + "_"
    if raw_allele_string.startswith(prefix):
        # Extract the part after the prefix
        return raw_allele_string[len(prefix):]
    
    # Fallback: if the string doesn't match the pattern (e.g., a raw sequence),
    # we return it as-is, though "NA" might be safer.
    # Let's return "NA" for any un-patterned data to keep the output clean.
    # If you want to keep 'TTATACACAT...', change this to: return raw_allele_string
    if not raw_allele_string.replace('_', '').isalnum():
        return "NA" # Catches full sequences

    # Default return if it's some other format
    return raw_allele_string

def convert_to_wide_format(input_filename, output_filename):
    """
    Reads the long-format likelihoods file and converts it to
    the wide, two-column genotype format.
    """
    
    # This nested dictionary will hold our restructured data:
    # {sample_id: {locus_name: (allele1, allele2)}}
    genotype_data = collections.defaultdict(dict)
    
    # We use these to maintain a consistent order in the output file
    all_loci = set()
    # Use OrderedDict.fromkeys to get unique samples in their original order
    all_samples = collections.OrderedDict() 

    print(f"Reading and processing '{input_filename}'...")

    try:
        with open(input_filename, mode='r', encoding='utf-8') as infile:
            reader = csv.DictReader(infile)
            
            for row in reader:
                sample = row['Sample']
                locus = row['Loci']
                
                # Add to our ordered sets
                all_loci.add(locus)
                all_samples[sample] = None # Add sample to the ordered dict

                allele_1_raw = row['Allele1']
                allele_2_raw = row['Allele2']
                is_homozygous = (row['like_AA'] == '1')
                
                # Determine the genotype based on the user's logic
                allele_name_1 = extract_allele_name(locus, allele_1_raw)
                
                if is_homozygous:
                    # Genotype is AA (Allele1, Allele1)
                    allele_name_2 = allele_name_1
                else:
                    # Genotype is AB (Allele1, Allele2)
                    allele_name_2 = extract_allele_name(locus, allele_2_raw)
                
                # Store the genotype
                genotype_data[sample][locus] = (allele_name_1, allele_name_2)

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    # Now, write the processed data to the new CSV file
    
    # Get a sorted list of loci for the header
    sorted_loci = sorted(list(all_loci))
    
    print(f"Writing data to '{output_filename}'...")
    
    try:
        with open(output_filename, mode='w', newline='', encoding='utf-8') as outfile:
            writer = csv.writer(outfile)
            
            # 1. Write the Header Row
            header = ['Indiv']
            for locus in sorted_loci:
                header.append(locus)  # Locus name
                header.append('')     # Blank column
            writer.writerow(header)
            
            # 2. Write the Data Rows
            for sample in all_samples.keys():
                data_row = [sample]
                for locus in sorted_loci:
                    # Get the genotype, default to ('NA', 'NA') if missing
                    genotype = genotype_data[sample].get(locus, ('NA', 'NA'))
                    data_row.append(genotype[0])
                    data_row.append(genotype[1])
                
                writer.writerow(data_row)
                
    except Exception as e:
        print(f"An error occurred while writing the file: {e}")
        return

    print(f"\nSuccess! Conversion complete. Output file is '{output_filename}'")

# --- Main execution ---
if __name__ == "__main__":
    INPUT_FILE = "../new_output/likelihoods_AA_AB_20251112_203730.csv"
    OUTPUT_FILE = "../new_output/genotype_table_formatted.csv"

    convert_to_wide_format(INPUT_FILE, OUTPUT_FILE)