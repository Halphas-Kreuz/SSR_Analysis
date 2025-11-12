import os
import sys

# --- CONFIGURATION ---
# IMPORTANT: Change this to the path of your data folder!
# Use '.' to run in the same directory as the script.
# Example (Mac/Linux): '/Users/yourname/Desktop/data'
# Example (Windows): 'C:\\Users\\yourname\\Desktop\\data'
DIRECTORY_PATH = '../samples'
# ---------------------


def clean_and_rename_files(directory):
    """
    Deletes I1/I2 files and renames R1/R2 files in a specified directory.
    """
    print(f"--- 1. Deletion Pass in: {directory} ---")
    deleted_count = 0
    
    # We must get a list of files first
    try:
        all_files = os.listdir(directory)
    except FileNotFoundError:
        print(f"ERROR: Directory not found: {directory}")
        print("Please update the 'DIRECTORY_PATH' variable in the script.")
        return
    except Exception as e:
        print(f"An error occurred listing files: {e}")
        return

    # --- Pass 1: Deletion (NOW DELETING I1 and I2) ---
    for filename in all_files:
        if not filename.endswith('.fastq'):
            continue  # Skip non-fastq files

        if '_I1_' in filename or '_I2_' in filename:
            try:
                file_path = os.path.join(directory, filename)
                os.remove(file_path)
                print(f"DELETED: {filename}")
                deleted_count += 1
            except Exception as e:
                print(f"Error deleting {filename}: {e}")
    
    if deleted_count == 0:
        print("No I1 or I2 files found to delete.")
    else:
        print(f"Total files deleted: {deleted_count}")

    # --- Pass 2: Renaming (NOW RENAMING R1 and R2) ---
    print(f"\n--- 2. Renaming Pass in: {directory} ---")
    renamed_count = 0
    
    # We must re-list the directory to get the remaining files
    try:
        remaining_files = os.listdir(directory)
    except Exception as e:
        print(f"An error occurred re-listing files: {e}")
        return

    for filename in remaining_files:
        # Skip files that don't match the expected pattern
        if not filename.endswith('.fastq') or '_L001' not in filename:
            continue

        try:
            # Find the base name by splitting at _L001
            base_name = filename.split('_L001')[0]
            new_name = None

            # Case A: Rename _R1_ files to _A
            if '_R1_' in filename:
                new_name = f"{base_name}_A.fastq"
            
            # Case B: Rename _R2_ files to _B
            elif '_R2_' in filename:
                new_name = f"{base_name}_B.fastq"

            # If a new name was determined, proceed with renaming
            if new_name:
                old_path = os.path.join(directory, filename)
                new_path = os.path.join(directory, new_name)
                
                # Safety check to avoid overwriting files
                if os.path.exists(new_path):
                    print(f"SKIPPED (target already exists): {filename} -> {new_name}")
                    continue
                
                os.rename(old_path, new_path)
                print(f"RENAMED: {filename} -> {new_name}")
                renamed_count += 1
                
        except Exception as e:
            print(f"Error renaming {filename}: {e}")

    if renamed_count == 0:
        print("No R1 or R2 files found to rename.")
    else:
        print(f"Total files renamed: {renamed_count}")

    print("\n--- Process Finished ---")


# --- Main execution ---
if __name__ == "__main__":
    # Safety check if the user is running in the current directory
    if DIRECTORY_PATH == '.':
        print("---------------------------------------------------------")
        print("WARNING: Script is set to run in the current directory ('.')")
        print("This will DELETE and RENAME files in the same folder")
        print("as this script. This action is IRREVERSIBLE.")
        print("---------------------------------------------------------")
        
        try:
            # Ask for confirmation before proceeding
            confirm = input("Are you sure you want to continue? (y/n): ").strip().lower()
            if confirm != 'y':
                print("Operation cancelled. Please edit the 'DIRECTORY_PATH' variable.")
                sys.exit() # Exit the script
        except EOFError:
             print("\nOperation cancelled.")
             sys.exit()

    clean_and_rename_files(DIRECTORY_PATH)