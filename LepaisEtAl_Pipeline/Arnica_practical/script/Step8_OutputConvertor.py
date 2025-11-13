import csv
import collections
import argparse
import sys

def clean_allele_name(raw_name):
    """
    Cleans a raw allele string, converting various 'no data'
    formats to a single 'NA' standard.
    """
    if raw_name in ["N/A", "No data", "Other sequences", "", None]:
        return "N/A"
    
    # Check for raw sequence data (non-alphanumeric after simple cleaning)
    # and convert to NA, as it's not a simple allele name.
    # This checks if the string, after removing _ and -, is alphanumeric.
    if not raw_name.replace('_', '').replace('-', '').isalnum():
        return "N/A"
    
    # Otherwise, return the full, original allele name
    return raw_name

def process_diploid_file(input_filename, output_filename):
    """
    Processes a diploid (ploidy 2) likelihoods file.
    Assumes 'like_AA', 'Allele1', 'Allele2' columns.
    Outputs 2 allele columns per locus.
    """
    print(f"Processing '{input_filename}' as DIPLOID (ploidy 2)...")
    # defaultdict simplifies adding new loci to a sample
    # OrderedDict preserves the order of samples as they appear in the file
    genotype_data = collections.defaultdict(dict)
    all_loci = set()
    all_samples = collections.OrderedDict()

    try:
        with open(input_filename, mode='r', encoding='utf-8') as infile:
            reader = csv.DictReader(infile)
            
            for row in reader:
                sample = row['Sample']
                locus = row['Loci']
                
                all_loci.add(locus)
                all_samples[sample] = None # Use OrderedDict as an ordered set

                allele_1_raw = row.get('Allele1', 'No data')
                allele_2_raw = row.get('Allele2', 'No data')
                # Use .get() for safety, default to '0' if column is missing
                is_homozygous = (row.get('like_AA', '0') == '1')
                
                allele_name_1 = clean_allele_name(allele_1_raw)
                
                if is_homozygous:
                    # For AA, both alleles are the same
                    allele_name_2 = allele_name_1
                else:
                    # For AB, use the second allele name
                    allele_name_2 = clean_allele_name(allele_2_raw)
                
                genotype_data[sample][locus] = (allele_name_1, allele_name_2)

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred reading the file: {e}")
        sys.exit(1)

    # Sort loci alphabetically for consistent column order
    sorted_loci = sorted(list(all_loci))
    
    print(f"Writing diploid data to '{output_filename}'...")
    try:
        with open(output_filename, mode='w', newline='', encoding='utf-8') as outfile:
            # Use tab as the delimiter
            writer = csv.writer(outfile, delimiter='\t')
            
            # Header: Indiv, Locus1, , Locus2, , ...
            header = ['Indiv']
            for locus in sorted_loci:
                header.append(locus)
                header.append('') # Add one blank column for diploid
            writer.writerow(header)
            
            # Data Rows
            for sample in all_samples.keys():
                data_row = [sample]
                for locus in sorted_loci:
                    # Get the genotype, or default to ('NA', 'NA') if missing
                    genotype = genotype_data[sample].get(locus, ('NA', 'NA'))
                    data_row.extend(genotype)
                writer.writerow(data_row)
                
    except Exception as e:
        print(f"An error occurred writing the file: {e}")
        sys.exit(1)

def process_triploid_file(input_filename, output_filename):
    """
    Processes a triploid (ploidy 3) likelihoods file.
    Assumes probability columns (Prob_AAA, Prob_AAB, Prob_ABC)
    and 'Allele1_Name', 'Allele2_Name', etc.
    Outputs 3 allele columns per locus.
    """
    print(f"Processing '{input_filename}' as TRIPLOID (ploidy 3)...")
    genotype_data = collections.defaultdict(dict)
    all_loci = set()
    all_samples = collections.OrderedDict()
    
    # Define the probability columns and their corresponding genotype structures
    PROB_COLS_3N = {
        'Prob_AAA': ['A', 'A', 'A'],
        'Prob_AAB': ['A', 'A', 'B'],
        'Prob_ABC': ['A', 'B', 'C']
    }
    
    try:
        with open(input_filename, mode='r', encoding='utf-8') as infile:
            reader = csv.DictReader(infile)
            
            for row in reader:
                sample = row['Sample']
                locus = row['Loci']
                
                all_loci.add(locus)
                all_samples[sample] = None 

                # Map allele placeholders (A,B,C) to actual allele names
                allele_map = {
                    'A': clean_allele_name(row.get('Allele1_Name', 'N/A')),
                    'B': clean_allele_name(row.get('Allele2_Name', 'N/A')),
                    'C': clean_allele_name(row.get('Allele3_Name', 'N/A')),
                }
                
                max_prob = -1.0
                winner_genotype_code = None
                
                # Find the highest probability
                for prob_col_name in PROB_COLS_3N.keys():
                    try:
                        # Get probability, default to 0.0 if missing
                        prob = float(row.get(prob_col_name, 0.0))
                        if prob > max_prob:
                            max_prob = prob
                            winner_genotype_code = PROB_COLS_3N[prob_col_name]
                    except (ValueError, TypeError):
                        continue # Skip if prob is not a valid number
                
                # Build the final 3-allele genotype
                final_genotype_list = ['NA'] * 3 # Default to all NA
                if max_prob > 0 and winner_genotype_code:
                    # Map 'A', 'B', 'C' codes to the actual allele names
                    final_genotype_list = [allele_map[code] for code in winner_genotype_code]
                
                genotype_data[sample][locus] = final_genotype_list

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred reading the file: {e}")
        sys.exit(1)

    sorted_loci = sorted(list(all_loci))
    
    print(f"Writing triploid data to '{output_filename}'...")
    try:
        with open(output_filename, mode='w', newline='', encoding='utf-8') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            
            # Header: Indiv, Locus1, , , Locus2, , , ...
            header = ['Indiv']
            for locus in sorted_loci:
                header.append(locus)
                header.extend([''] * 2) # 1 locus name + 2 blank columns
            writer.writerow(header)
            
            # Data Rows
            for sample in all_samples.keys():
                data_row = [sample]
                for locus in sorted_loci:
                    genotype = genotype_data[sample].get(locus, ['NA'] * 3)
                    data_row.extend(genotype)
                writer.writerow(data_row)
                
    except Exception as e:
        print(f"An error occurred writing the file: {e}")
        sys.exit(1)

def process_tetraploid_file(input_filename, output_filename):
    """
    Processes a tetraploid (ploidy 4) likelihoods file.
    Assumes probability columns (Prob_AAAA, Prob_AABC, etc.)
    and 'Allele1_Name', 'Allele2_Name', etc.
    Outputs 4 allele columns per locus.
    """
    print(f"Processing '{input_filename}' as TETRAPLOID (ploidy 4)...")
    genotype_data = collections.defaultdict(dict)
    all_loci = set()
    all_samples = collections.OrderedDict()
    
    # Define the probability columns and their corresponding genotype structures
    PROB_COLS_4N = {
        'Prob_AAAA': ['A', 'A', 'A', 'A'],
        'Prob_AAAB': ['A', 'A', 'A', 'B'],
        'Prob_AABB': ['A', 'A', 'B', 'B'],
        'Prob_AABC': ['A', 'A', 'B', 'C'],
        'Prob_ABCD': ['A', 'B', 'C', 'D']
    }
    
    try:
        with open(input_filename, mode='r', encoding='utf-8') as infile:
            reader = csv.DictReader(infile)
            
            for row in reader:
                sample = row['Sample']
                locus = row['Loci']
                
                all_loci.add(locus)
                all_samples[sample] = None 

                # Map allele placeholders (A,B,C,D) to actual allele names
                allele_map = {
                    'A': clean_allele_name(row.get('Allele1_Name', 'N/A')),
                    'B': clean_allele_name(row.get('Allele2_Name', 'N/A')),
                    'C': clean_allele_name(row.get('Allele3_Name', 'N/A')),
                    'D': clean_allele_name(row.get('Allele4_Name', 'N/A'))
                }
                
                max_prob = -1.0
                winner_genotype_code = None
                
                # Find the highest probability
                for prob_col_name in PROB_COLS_4N.keys():
                    try:
                        prob = float(row.get(prob_col_name, 0.0))
                        if prob > max_prob:
                            max_prob = prob
                            winner_genotype_code = PROB_COLS_4N[prob_col_name]
                    except (ValueError, TypeError):
                        continue # Skip if prob is not a valid number
                
                # Build the final 4-allele genotype
                final_genotype_list = ['NA'] * 4 # Default to all NA
                if max_prob > 0 and winner_genotype_code:
                    # Map 'A', 'B', 'C', 'D' codes to the actual allele names
                    final_genotype_list = [allele_map[code] for code in winner_genotype_code]
                
                genotype_data[sample][locus] = final_genotype_list

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred reading the file: {e}")
        sys.exit(1)

    sorted_loci = sorted(list(all_loci))
    
    print(f"Writing tetraploid data to '{output_filename}'...")
    try:
        with open(output_filename, mode='w', newline='', encoding='utf-8') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            
            # Header: Indiv, Locus1, , , , Locus2, , , , ...
            header = ['Indiv']
            for locus in sorted_loci:
                header.append(locus)
                header.extend([''] * 3) # 1 locus name + 3 blank columns
            writer.writerow(header)
            
            # Data Rows
            for sample in all_samples.keys():
                data_row = [sample]
                for locus in sorted_loci:
                    genotype = genotype_data[sample].get(locus, ['NA'] * 4)
                    data_row.extend(genotype)
                writer.writerow(data_row)
                
    except Exception as e:
        print(f"An error occurred writing the file: {e}")
        sys.exit(1)

def main():
    # --- HARDCODED FILE PATHS ---
    # !! Change these values to match your file paths !!
    INPUT_FILE = "../new_output/likelihoods.csv"
    OUTPUT_FILE = "../new_output/final_genotypes.csv"
    # ------------------------------

    parser = argparse.ArgumentParser(
        description="Converts a 'likelihoods' CSV file from long format to a wide, tab-delimited genotype format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""

"""
    )
    
    # Optional positional argument: Ploidy
    parser.add_argument(
        'ploidy', 
        type=int, 
        nargs='?',      # Makes the argument optional
        default=2,      # Default value if not provided
        choices=[2, 3, 4], # Restrict allowed values
        help='Optional: Expected ploidy (2, 3, or 4). Default: 2'
    )
    
    args = parser.parse_args()
    
    # Run the correct processing function based on the ploidy argument
    if args.ploidy == 2:
        process_diploid_file(INPUT_FILE, OUTPUT_FILE)
    elif args.ploidy == 3:
        process_triploid_file(INPUT_FILE, OUTPUT_FILE)
    elif args.ploidy == 4:
        process_tetraploid_file(INPUT_FILE, OUTPUT_FILE)
    else:
        # This check is technically redundant due to 'choices' but is good practice
        print(f"Error: Ploidy {args.ploidy} is not supported. Please use 2, 3, or 4.")
        sys.exit(1)
        
    print(f"\nSuccess! Conversion complete. Output file is '{OUTPUT_FILE}'")

if __name__ == "__main__":
    main()