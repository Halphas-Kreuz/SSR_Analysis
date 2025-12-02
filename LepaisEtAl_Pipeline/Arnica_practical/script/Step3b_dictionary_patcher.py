import pandas as pd
import sys
import os
import collections

# --- Configuration ---
# Usage: python Step3b_DictionaryPatcher.py [MODE]
# Default Mode: human

DICT_PATH = '../new_output/Step3_AlleleInfo_dictionary.csv'
PATCHED_DICT_PATH = '../new_output/Step3b_Patched_AlleleInfo_dictionary.csv'

def patch_dictionary():
    # 1. Parse Mode from Command Line
    mode = 'human'
    if len(sys.argv) > 1:
        mode = sys.argv[1].lower()
    
    # Dynamic Filename based on Mode
    named_file_path = f'../new_output/Step4_AlleleInfo_named_{mode}.csv'

    print(f"⚙️  Running Dictionary Patcher v2 (Priority Mode)")
    print(f"   Mode: {mode}")
    print(f"   Target Input: {named_file_path}")
    
    # 2. Load Dictionary
    if not os.path.exists(DICT_PATH):
        print(f"❌ Error: Dictionary not found at {DICT_PATH}")
        sys.exit(1)
        
    try:
        df_dict = pd.read_csv(DICT_PATH, dtype=str)
        print(f"   Loaded dictionary: {len(df_dict)} entries.")
    except Exception as e:
        print(f"❌ Error reading dictionary: {e}")
        sys.exit(1)

    # 3. Analyze Existing Codes (for auto-increment)
    locus_max_code = collections.defaultdict(int)
    locus_len_count = collections.defaultdict(int) 

    for idx, row in df_dict.iterrows():
        locus = str(row['Locus_Name'])
        human_id = str(row['Human_ID'])
        length = str(row['Length'])
        
        # Track Max Code
        orig_code = str(row['Original_AlleleSeqCode'])
        if orig_code.isdigit():
            val = int(orig_code)
            if val > locus_max_code[locus]:
                locus_max_code[locus] = val

        # Track Count (Arm01_110_2)
        parts = human_id.split('_')
        current_c = 1
        if len(parts) >= 3 and parts[-1].isdigit() and parts[-1] != length:
             current_c = int(parts[-1])
        
        key = (locus, length)
        if current_c > locus_len_count[key]:
            locus_len_count[key] = current_c

    # 4. Load the Named File
    if not os.path.exists(named_file_path):
        print(f"❌ Error: Input file not found at {named_file_path}")
        print(f"   (Did you run Step 4 with mode '{mode}' first?)")
        sys.exit(1)
        
    df_named = pd.read_csv(named_file_path, dtype=str)
    
    allele_cols = [c for c in df_named.columns if c.startswith('Allele')]
    
    # Cache existing sequences (normalized)
    existing_sequences = set(df_dict['Full_Sequence'].dropna().str.strip().str.upper())
    
    new_entries = []
    seen_orphans = set() 

    # Debug counters
    count_scanned = 0
    count_ignored_short = 0
    count_ignored_chars = 0
    count_ignored_exists = 0

    print(f"   Scanning {len(df_named)} rows for orphan sequences...")

    for idx, row in df_named.iterrows():
        locus = row.get('Loci') 
        if pd.isna(locus): continue
        locus = str(locus).strip()
        
        for col in allele_cols:
            val = row[col]
            if pd.isna(val): continue
            
            # CLEANUP
            val = str(val).strip().replace('"', '').replace("'", "")
            val_upper = val.upper()
            
            count_scanned += 1
            
            # --- DETECTION LOGIC ---
            if '_' in val: continue # Looks like a name already
            if len(val) < 20: 
                count_ignored_short += 1
                continue
            if not all(c in 'ACGTN' for c in val_upper):
                count_ignored_chars += 1
                continue
            if val_upper in existing_sequences:
                count_ignored_exists += 1
                continue
            if (locus, val_upper) in seen_orphans: 
                continue
            
            seen_orphans.add((locus, val_upper))
            
            # --- REGISTER NEW ORPHAN ---
            seq_len = str(len(val))
            
            # Count Logic
            key_len = (locus, seq_len)
            new_count = locus_len_count[key_len] + 1
            locus_len_count[key_len] = new_count
            
            # Code Logic (Increment)
            current_max = locus_max_code[locus]
            if current_max == 0: current_max = 5000
            new_code = current_max + 1
            locus_max_code[locus] = new_code
            
            human_id = f"{locus}_{seq_len}_{new_count}"
            genalex_id = f"{seq_len}.{new_count}"
            original_code = str(new_code)
            fused_index = f"{locus}_{original_code}"
            
            new_row = {
                'Fused_Index': fused_index,
                'Human_ID': human_id,
                'GenAlEx_ID': genalex_id,
                'Original_AlleleSeqCode': original_code,
                'Full_Sequence': val,
                'Repeat_Pattern': val, 
                'Length': seq_len,
                'Locus_Name': locus
            }
            new_entries.append(new_row)

    # 5. Save Updates (PREPENDING for Priority)
    if new_entries:
        print(f"⚡ Found {len(new_entries)} new orphan sequences!")
        ex = new_entries[0]
        print(f"   Example: {ex['Locus_Name']} -> {ex['Human_ID']} (Code {ex['Original_AlleleSeqCode']})")
        
        df_new = pd.DataFrame(new_entries)
        
        # New entries go FIRST
        df_combined = pd.concat([df_new, df_dict], ignore_index=True)
        
        df_combined.to_csv(PATCHED_DICT_PATH, index=False)
        print(f"✅ Dictionary updated! Saved to {PATCHED_DICT_PATH}")
    else:
        print("✅ No orphan sequences found. Dictionary is up to date.")

if __name__ == "__main__":
    patch_dictionary()