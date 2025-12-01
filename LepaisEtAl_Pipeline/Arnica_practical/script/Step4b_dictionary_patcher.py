import pandas as pd
import sys
import os
import collections

# --- Configuration ---
DICT_PATH = '../new_output/AlleleInfo_dictionary.csv'
NAMED_FILE_PATH = '../new_output/AlleleInfo_named_human.csv' 

def patch_dictionary():
    print(f"⚙️  Running Dictionary Patcher (Step 3b) - Smart Increment Mode...")
    
    # 1. Load Dictionary
    if not os.path.exists(DICT_PATH):
        print(f"❌ Error: Dictionary not found at {DICT_PATH}")
        sys.exit(1)
        
    try:
        # Read as string to preserve formatting, but we'll convert codes to int for math
        df_dict = pd.read_csv(DICT_PATH, dtype=str)
        print(f"   Loaded dictionary: {len(df_dict)} entries.")
    except Exception as e:
        print(f"❌ Error reading dictionary: {e}")
        sys.exit(1)

    # 2. Analyze Existing Codes & Counts (Per Locus)
    # We need two trackers:
    # A. Max Original Code per Locus (to increment 105 -> 106)
    # B. Max Human ID Count per Locus_Length (to increment .2 -> .3)
    
    locus_max_code = collections.defaultdict(int)
    locus_len_count = collections.defaultdict(int) # Key: (Locus, Length)

    for idx, row in df_dict.iterrows():
        locus = row['Locus_Name']
        length = str(row['Length'])
        human_id = str(row['Human_ID'])
        orig_code_str = str(row['Original_AlleleSeqCode'])
        
        # A. Track Max Original Code
        if pd.notna(locus) and orig_code_str.isdigit():
            val = int(orig_code_str)
            if val > locus_max_code[locus]:
                locus_max_code[locus] = val

        # B. Track Human ID Count (e.g. Arm01_110_2 -> count 2)
        # Format: Locus_Length_Count
        parts = human_id.split('_')
        current_c = 1
        # Heuristic: If last part is a digit and != length, treat as count
        if len(parts) >= 3 and parts[-1].isdigit() and parts[-1] != length:
             current_c = int(parts[-1])
        
        key = (locus, length)
        if current_c > locus_len_count[key]:
            locus_len_count[key] = current_c

    # 3. Load the Named File to find Orphans
    if not os.path.exists(NAMED_FILE_PATH):
        print(f"❌ Error: Input file not found at {NAMED_FILE_PATH}")
        sys.exit(1)
        
    df_named = pd.read_csv(NAMED_FILE_PATH, dtype=str)
    
    allele_cols = [c for c in df_named.columns if c.startswith('Allele')]
    existing_sequences = set(df_dict['Full_Sequence'].dropna().unique())
    
    new_entries = []
    seen_orphans = set() 

    print(f"   Scanning {len(df_named)} rows for orphan sequences...")

    for idx, row in df_named.iterrows():
        locus = row.get('Loci') 
        if pd.isna(locus): continue
        
        for col in allele_cols:
            val = row[col]
            if pd.isna(val): continue
            val = val.strip()
            
            # --- DETECTION LOGIC (Same as before) ---
            # Length > 20, ACGTN only, no underscores
            if len(val) > 20 and all(c.upper() in 'ACGTN' for c in val) and '_' not in val:
                
                if val in existing_sequences: continue
                if (locus, val) in seen_orphans: continue
                
                seen_orphans.add((locus, val))
                
                # --- REGISTER NEW ORPHAN ---
                seq_len = str(len(val))
                
                # 1. Determine Count for Human/GenAlEx ID
                key_len = (locus, seq_len)
                new_count = locus_len_count[key_len] + 1
                locus_len_count[key_len] = new_count # Update tracker
                
                human_id = f"{locus}_{seq_len}_{new_count}"
                genalex_id = f"{seq_len}.{new_count}"
                
                # 2. Determine Original Code (Incrementing per Locus)
                # If Locus has no codes yet, default to a safe start like 5000 or 1
                # But usually it has some. We take max + 1.
                current_max_code = locus_max_code[locus]
                if current_max_code == 0:
                    current_max_code = 5000 # Fallback if locus is totally new
                
                new_code_val = current_max_code + 1
                locus_max_code[locus] = new_code_val # Update tracker