import os

def process_file(input_path, output_path):
    if not os.path.exists(input_path):
        print(f"Error: Input file not found at '{input_path}'")
        return

    try:
        with open(input_path, 'r', encoding='utf-8') as infile, \
             open(output_path, 'w', encoding='utf-8') as outfile:
            
            for line in infile:
                columns = line.split()
                if columns:
                    new_line = ' '.join(columns)
                    outfile.write(new_line + '\n')
        
        print(f"Success: Processed '{input_path}' and saved to '{output_path}'")

    except Exception as e:
        print(f"An error occurred: {e}")


input_folder_path = "../"
output_folder_path = "../new_output"
input_filename = "AlleleInformationFile_nSSR_FullLength_ParameterSet1_sa50_sb10_m15_n20.txt"
output_filename = "AlleleInformationFile_nSSR_FullLength_ParameterSet1_sa50_sb10_m15_n20.csv"


full_input_path = os.path.join(input_folder_path, input_filename)
full_output_path = os.path.join(output_folder_path, output_filename)

process_file(full_input_path, full_output_path)
