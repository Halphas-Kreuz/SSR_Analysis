import csv
import collections
import argparse
import sys
import os

# --- Step 8: Final Output Converter ---
# Converts the probabilistic likelihoods file (from Step 7) 
# into a standard "Wide Format" Genotype Table for downstream analysis.

def clean_allele_name(raw_name):
    """
    Cleans a raw allele string.
    """
    if raw_name in ["N/A", "No data", "Other sequences", "", "nan", None]:
        return "N/A"
    
    # Check if alphanumeric (allowing - and _)
    if not raw_name.replace('_', '').replace('-', '').replace('.', '').isalnum():
        return "N/A"
    
    return raw_name

def process_file_generic(input_filename, output_filename, ploidy):
    """
    Generic processor for ANY ploidy (2, 3, or 4).
    Finds the genotype hypothesis with the maximum probability.
    """
    print(f"⚙️  Processing Ploidy {ploidy} data...")
    print(f"   Input:  {input_filename}")
    
    # 1. Define Hypotheses Maps based on Step 7 output headers
    # Maps Column Name -> List of Allele Placeholders
    if ploidy == 2:
        PROB_COLS = {
            'Prob_AA': ['A', 'A'],
            'Prob_AB': ['A', 'B']
        }
    elif ploidy == 3:
        PROB_COLS = {
            'Prob_AAA': ['A', 'A', 'A'],
            'Prob_AAB': ['A', 'A', 'B'],
            'Prob_ABC': ['A', 'B', 'C']
        }
    elif ploidy == 4:
        PROB_COLS = {
            'Prob_AAAA': ['A', 'A', 'A', 'A'],
            'Prob_AAAB': ['A', 'A', 'A', 'B'],
            'Prob_AABB': ['A', 'A', 'B', 'B'],
            'Prob_AABC': ['A', 'A', 'B', 'C'],
            'Prob_ABCD': ['A', 'B', 'C', 'D']
        }
    else:
        print(f"❌ Error: Ploidy {ploidy} not supported.")
        sys.exit(1)

    genotype_data = collections.defaultdict(dict)
    all_loci = set()
    all_samples = collections.OrderedDict()

    try:
        with open(input_filename, mode='r', encoding='utf-8') as infile:
            reader = csv.DictReader(infile)
            
            # Verify critical columns exist
            fieldnames = reader.fieldnames
            if not fieldnames or 'Allele1_Name' not in fieldnames:
                 print(f"❌ Error: Input file format mismatch. Could not find 'Allele1_Name'.")
                 print(f"   Headers found: {fieldnames[:5]}...")
                 sys.exit(1)

            for row in reader:
                sample = row.get('Sample') or row.get('Sample_Name')
                locus = row.get('Locus') or row.get('Loci')
                
                if not sample or not locus: continue

                all_loci.add(locus)
                all_samples[sample] = None 

                # 2. Map Allele Placeholders (A, B, C...) to Actual Names
                # Step 7 outputs: Allele1_Name, Allele2_Name, ...
                allele_map = {}
                chars = ['A', 'B', 'C', 'D', 'E']
                for i in range(ploidy):
                    col_name = f"Allele{i+1}_Name"
                    val = row.get(col_name, 'N/A')
                    allele_map[chars[i]] = clean_allele_name(val)

                # 3. Find Winner (Max Probability)
                max_prob = -1.0
                winner_structure = None
                
                # Iterate through expected probability columns (e.g. Prob_AA, Prob_AB)
                for col_name, structure in PROB_COLS.items():
                    try:
                        val_str = row.get(col_name)
                        if val_str:
                            prob = float(val_str)
                            if prob > max_prob:
                                max_prob = prob
                                winner_structure = structure
                    except ValueError:
                        continue

                # 4. Construct Final Genotype List
                final_genotype = ['N/A'] * ploidy # Default
                
                # Threshold: You can set a minimum confidence here (e.g. > 0.0)
                if max_prob > 0.0 and winner_structure:
                    # Map ['A', 'B'] -> ['Arm01_83', 'Arm01_92']
                    final_genotype = [allele_map[code] for code in winner_structure]
                
                genotype_data[sample][locus] = final_genotype

    except FileNotFoundError:
        print(f"❌ Error: Input file not found: {input_filename}")
        print("   Did you run Step 7 with the same Ploidy and Mode?")
        sys.exit(1)

    # Write Output
    sorted_loci = sorted(list(all_loci))
    print(f"💾 Writing {len(all_samples)} samples and {len(sorted_loci)} loci to '{output_filename}'...")
    
    try:
        with open(output_filename, mode='w', newline='', encoding='utf-8') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            
            # Header Construction: Indiv, Loc1, , Loc2, , ...
            # We add (Ploidy - 1) blank columns after each Locus name
            header = ['Indiv']
            for locus in sorted_loci:
                header.append(locus)
                header.extend([''] * (ploidy - 1))
            writer.writerow(header)
            
            # Data Rows
            for sample in all_samples.keys():
                data_row = [sample]
                for locus in sorted_loci:
                    # Get genotype or default NA
                    gt = genotype_data[sample].get(locus, ['N/A'] * ploidy)
                    data_row.extend(gt)
                writer.writerow(data_row)
        
        print(f"✅ Success! Final file ready: {output_filename}")

    except Exception as e:
        print(f"❌ Error writing file: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Step 8: Final Output Converter")
    
    parser.add_argument('ploidy', type=int, nargs='?', default=2, 
                        help='Ploidy level (2, 3, or 4). Default: 2')
    
    parser.add_argument('mode', type=str, nargs='?', default='human',
                        help="Mode used in previous steps ('human', 'genalex', 'original'). Default: human")

    args = parser.parse_args()
    
    # Construct Filenames Dynamically
    # Matches Step 7 output: ../new_output/Step7_likelihoods_Ploidy2_human.csv
    input_file = f"../new_output/Step7_likelihoods_Ploidy{args.ploidy}_{args.mode}.csv"
    output_file = f"../new_output/Step8_final_genotypes_Ploidy{args.ploidy}_{args.mode}.txt"

    process_file_generic(input_file, output_file, args.ploidy)

if __name__ == "__main__":
    main()