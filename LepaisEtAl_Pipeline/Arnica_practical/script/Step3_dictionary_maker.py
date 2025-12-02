import csv
import os
import sys
import glob

# --- Configuration ---
# Match any AlleleInformationFile txt
DEFAULT_INPUT_PATTERN = '../AlleleInformationFile_*.txt'
OUTPUT_FILE = '../new_output/Step3_AlleleInfo_dictionary.csv'

def txt_to_csv(txt_path, csv_path):
    """Converts a tab-delimited text file to a CSV file."""
    if not os.path.exists(csv_path):
        try:
            with open(txt_path, 'r') as infile, open(csv_path, 'w', newline='') as outfile:
                reader = csv.reader(infile, delimiter='\t')
                writer = csv.writer(outfile)
                for row in reader:
                    writer.writerow(row)
            print(f"✅ Converted {txt_path} to {csv_path}")
        except Exception as e:
            print(f"❌ Error converting TXT to CSV: {e}")
            sys.exit(1)
    else:
        print(f"ℹ️  CSV file already exists: {csv_path}")

def load_and_transform_data(csv_path):
    """
    Loads data and creates THREE naming conventions, plus preserves the Annotated Sequence.
    Input Columns based on user sample:
    0: Locus
    1: AlleleSequenceAnnotated (e.g., ...CATA(13)...) -> WE NEED THIS
    2: AlleleSeqCode (Original ID)
    3: Occurances
    4: AlleleSequence (Full)
    5: AlleleLength
    """
    try:
        with open(csv_path, 'r') as file:
            reader = csv.reader(file)
            table = [row for row in reader]
            
        if not table:
            return []

        # Skip header if it exists
        data_rows = table[1:] 

        matrix = []
        # Dictionary to track ID occurrences: Key = "LocusName_Length"
        id_counts = {}

        for row in data_rows:
            # Ensure row has required columns
            if len(row) >= 5:
                locus_name = row[0]
                annotated_seq = row[1] # <--- NEW: Extract Annotated Sequence
                original_id = row[2] 
                sequence = row[4]
                
                # Calculate or get length
                # Using python len() is safer than trusting the file column if it varies
                seq_length = len(sequence) 
                
                # --- The Logic for 3 IDs ---
                
                # Base Key for counting duplicates (e.g. "Arm11_110")
                base_key = f"{locus_name}_{seq_length}"

                if base_key in id_counts:
                    id_counts[base_key] += 1
                    count = id_counts[base_key]
                    
                    # 1. Human ID: Arm11_110_2
                    human_id = f"{locus_name}_{seq_length}_{count}"
                    
                    # 2. GenAlEx ID: 110.2
                    genalex_id = f"{seq_length}.{count}"
                    
                else:
                    id_counts[base_key] = 1
                    count = 1
                    
                    # 1. Human ID: Arm11_110
                    human_id = f"{locus_name}_{seq_length}"
                    
                    # 2. GenAlEx ID: 110
                    genalex_id = f"{seq_length}"

                # Key used to match Step 2 data
                fused_index = f"{locus_name}_{original_id}"

                # Assemble Output Row
                output_row = [
                    fused_index,   # 0
                    human_id,      # 1
                    genalex_id,    # 2
                    original_id,   # 3
                    sequence,      # 4
                    annotated_seq, # 5 <--- NEW: Inserted here
                    seq_length,    # 6
                    locus_name     # 7
                ]
                matrix.append(output_row)
        
        return matrix

    except Exception as e:
        print(f"❌ An error occurred: {e}")
        return []

if __name__ == "__main__":
    # 1. Find Input
    found_files = glob.glob(DEFAULT_INPUT_PATTERN)
    if not found_files:
        print(f"❌ No input file found matching: {DEFAULT_INPUT_PATTERN}")
        sys.exit(1)
    
    txt_path = found_files[0]
    csv_path = txt_path.replace('.txt', '.csv')

    # 2. Convert
    # Note: If the CSV already exists from a previous run but is 'bad', 
    # we might want to force regenerate it. But strictly strictly speaking, 
    # the python script reads the CSV. If you updated the TXT, delete the old CSV first.
    if os.path.exists(csv_path):
        os.remove(csv_path) # Force remove old CSV to ensure fresh conversion
        
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    txt_to_csv(txt_path, csv_path)

    # 3. Process
    print(f"⚙️  Processing dictionary from {txt_path}...")
    matrix = load_and_transform_data(csv_path)

    # 4. Write Output
    if matrix:
        with open(OUTPUT_FILE, 'w', newline='') as outfile:
            writer = csv.writer(outfile)
            # Updated Header
            writer.writerow([
                'Fused_Index', 
                'Human_ID', 
                'GenAlEx_ID', 
                'Original_AlleleSeqCode', 
                'Full_Sequence', 
                'Repeat_Pattern', # <--- NEW Header
                'Length', 
                'Locus_Name'
            ])
            writer.writerows(matrix)
        print(f"✅ Dictionary created: {OUTPUT_FILE}")
        print("   -> Added column: Repeat_Pattern")
    else:
        print("⚠️  No data processed.")