import subprocess
import sys
import os

# --- Helper Function ---
def run_command(command, step_name):
    """
    Runs a shell command and checks for errors.
    If an error occurs, it prints a message and exits the script.
    """
    print(f"\n--- Starting: {step_name} ---")
    print(f"Running command: {command}")
    
    try:
        # Run the command.
        # check=True: This is CRITICAL. It makes the script automatically
        #              fail if the command returns a non-zero exit code.
        # shell=True: Allows us to pass the command as a single string.
        # text=True:  Shows output as text (not bytes).
        subprocess.run(
            command, 
            shell=True, 
            check=True, 
            text=True, 
            stdout=sys.stdout, 
            stderr=sys.stderr
        )
        
        print(f"--- Finished: {step_name} ---")
        
    except subprocess.CalledProcessError as e:
        print(f"\n**************************************************")
        print(f"ERROR: Step '{step_name}' FAILED!")
        print(f"Failed Command: {e.cmd}")
        print(f"Return Code: {e.returncode}")
        print("**************************************************")
        
        # Exit the entire helper script. We don't want to continue
        # if a step in the pipeline failed.
        sys.exit(1)

# ===================================================================
# --- MAIN PIPELINE SCRIPT ---
# ===================================================================

def main():
    print("=======================================")
    print("  Starting Analysis Pipeline Helper  ")
    print("=======================================")

    # --- Check Data Naming ---
    # We ask the user if the data is already clean.
    data_is_clean = False
    while True:
        user_input = input("Are your data files already named correctly? (y/n): ").strip().lower()
        if user_input == 'y':
            data_is_clean = True
            break
        elif user_input == 'n':
            data_is_clean = False
            break
        else:
            print("Invalid input. Please enter 'y' or 'n'.")

    # --- Get Ploidy Number ---
    # We ask the user for the ploidy number and make sure it's a valid integer.
    ploidy_number = None
    while True:
        user_input = input("Please enter the ploidy number (e.g., 2): ").strip()
        try:
            # Convert the input to an integer
            ploidy_number = int(user_input)
            if ploidy_number > 0:
                print(f"Ploidy set to: {ploidy_number}")
                break
            else:
                print("Ploidy must be a positive number.")
        except ValueError:
            print("Invalid input. Please enter a whole number.")

    # --- Define and Run the Pipeline ---
    
    # This is the list of commands to run, in order.
    pipeline_steps = []

    # --- Step 0: Optional Name Cleaner ---
    if not data_is_clean:
        pipeline_steps.append({
            "name": "Step 0: Clean Data Names",
            "command": "python3 Step0_NameCleaning.py"
        })
    else:
        print("\nSkipping Step 0 (Name Cleaning) as requested.")

    # --- Sequential Pipeline Steps ---
    #
    # !!! IMPORTANT !!!
    # You MUST edit the paths and script names below to match
    # your actual file structure.
    #
    
    # --- Step 1
    pipeline_steps.append({
        "name": "Step 1: filter the necessary line",
        "command": f"bash Step1_filter_2_line.sh {ploidy_number}"
    })
    
    # --- Step 2
    pipeline_steps.append({
        "name": "Step 2: database generation",
        "command": f"python3 Step2_improved_generator.py {ploidy_number}"
    })

    # --- Step 3:
    pipeline_steps.append({
        "name": "Step 3: write the dictionary",
        "command": "python3 Step3_dictionary_maker.py"
    })

    # --- Step 4
    pipeline_steps.append({
        "name": "Step 4: Alleleinfo naming ",
        "command": f"python3 Step4_AlleleInfo_naming.py {ploidy_number}"
    })
    
    # --- Step 5: 
    pipeline_steps.append({
        "name": "Step 5: Old Genetype table translation ",
        "command": "python3 Step5_GeneTable_unique.py"
    })

    # --- Step 6: 
    pipeline_steps.append({
        "name": "Step 6: Stutter Table construction",
        "command": "python3 Step6_StutterTable.py"
    })

    # --- Step 7:
    pipeline_steps.append({
        "name": "Step 6a: preparation for the final bayesian calculation",
        "command": "python3 Step6a_BayesianPredecessor.py"
    })

    # --- Step 8:
    pipeline_steps.append({
        "name": "Step 7: bayesian",
        "command": f"python3 Step7_bayesianMultiploidy.py {ploidy_number}"
    })

    # --- Step 9
    pipeline_steps.append({
        "name": "Step 8: output formatting",
        "command": f"python3 Step8_OutputConvertor.py {ploidy_number}"
    })


    # --- Execute the Pipeline ---
    # Now, we loop through our defined steps and run them one by one.
    for step in pipeline_steps:
        run_command(step["command"], step["name"])

    # --- Finish ---
    print("\n=======================================")
    print("  Pipeline Completed Successfully!  ")
    print("=======================================")


# This makes sure the `main()` function is called when
# you run: python pipeline_helper.py
if __name__ == "__main__":
    main()