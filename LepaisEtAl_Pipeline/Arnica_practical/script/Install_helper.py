import subprocess
import sys
import time
import os

def run_pipeline():
    print("=======================================")
    print("  Starting Analysis Pipeline Helper  ")
    print("=======================================")

    # --- 1. Check Data Naming ---
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

    # --- 2. Get Ploidy Number ---
    ploidy_number = None
    while True:
        user_input = input("Please enter the ploidy number (e.g., 2): ").strip()
        try:
            ploidy_number = int(user_input)
            if ploidy_number > 0:
                print(f"Ploidy set to: {ploidy_number}")
                break
            else:
                print("Ploidy must be a positive number.")
        except ValueError:
            print("Invalid input. Please enter a whole number.")

    # --- 3. Get Mode (Default to human) ---
    mode = 'human'
    while True:
        user_input = input("Select naming mode [human/genalex/original] (default: human): ").strip().lower()
        if user_input == '':
            mode = 'human'
            break
        elif user_input in ['human', 'genalex', 'original']:
            mode = user_input
            break
        else:
            print("Invalid mode. Please choose 'human', 'genalex', or 'original'.")
    print(f"Mode set to: {mode}")

    print("\n---------------------------------------")
    print(" Initializing Pipeline...")
    print("---------------------------------------")

    # --- Define Pipeline Steps ---
    
    # This is the list of commands to run, in order.
    pipeline_steps = []

    # --- Step 0: Optional Name Cleaner ---
    if not data_is_clean:
        pipeline_steps.append({
            "name": "Step 0: Clean Data Names",
            "command": "python3 Step0_NameCleaning.py"
        })
    else:
        print(" Skipping Step 0 (Name Cleaning) as requested.")

    # --- Step 1: Filter lines (Bash script)
    # Ensure Step1_filter_2_line.sh exists and is executable
    pipeline_steps.append({
        "name": "Step 1: Filter necessary lines (Bash)",
        "command": f"bash Step1_filter_2_line.sh {ploidy_number}"
    })
    
    # --- Step 2: Database generation
    pipeline_steps.append({
        "name": "Step 2: Database Generation",
        "command": f"python3 Step2_improved_generator.py {ploidy_number}"
    })

    # --- Step 3: Dictionary Creation
    pipeline_steps.append({
        "name": "Step 3: Dictionary Creation",
        "command": "python3 Step3_dictionary_maker.py"
    })

    # --- Step 4 (Pass 1): Initial Naming
    pipeline_steps.append({
        "name": "Step 4 (Pass 1): Initial Naming",
        "command": f"python3 Step4_AlleleInfo_naming.py {ploidy_number} {mode}"
    })

    # --- Step 3b: Dictionary Patcher
    pipeline_steps.append({
        "name": "Step 3b: Dictionary Patcher",
        "command": f"python3 Step3b_dictionary_patcher.py {mode}"
    })

    # --- Step 4 (Pass 2): Final Naming
    pipeline_steps.append({
        "name": "Step 4 (Pass 2): Final Naming",
        "command": f"python3 Step4_AlleleInfo_naming.py {ploidy_number} {mode}"
    })

    # --- Step 5: GeneTable Conversion
    pipeline_steps.append({
        "name": "Step 5: Unique GeneTable",
        "command": f"python3 Step5_GeneTable_unique.py {mode}"
    })

    # --- Step 6: Stutter Table Calculation
    pipeline_steps.append({
        "name": "Step 6: Stutter Table",
        "command": f"python3 Step6_StutterTable.py 0.5 0.1 {mode}"
    })

    # --- Step 6a: Bayesian Predecessor
    pipeline_steps.append({
        "name": "Step 6a: Bayesian Predecessor",
        "command": f"python3 Step6a_BayesianPredecessor.py {mode}"
    })

    # --- Step 7: Bayesian Likelihoods
    pipeline_steps.append({
        "name": "Step 7: Bayesian Genotype Likelihoods",
        "command": f"python3 Step7_bayesianMultiploidy.py {ploidy_number} {mode}"
    })

    # --- Step 8: Final Output Conversion
    pipeline_steps.append({
        "name": "Step 8: Final Output Converter",
        "command": f"python3 Step8_OutputConvertor.py {ploidy_number} {mode}"
    })

    # ==========================================
    # 🏃 EXECUTION LOOP
    # ==========================================
    
    total_start_time = time.time()

    for step in pipeline_steps:
        step_name = step["name"]
        cmd_str = step["command"]
        
        print(f"\n🔹 [Running] {step_name}...")
        # print(f"   Command: {cmd_str}") # Uncomment for debugging
        
        step_start = time.time()
        try:
            # shell=True allows running simple string commands like "bash script.sh arg"
            subprocess.run(cmd_str, shell=True, check=True)
            
            elapsed = time.time() - step_start
            print(f"✅ [Success] {step_name} ({elapsed:.2f}s)")
            
        except subprocess.CalledProcessError:
            print(f"\n❌ [FAILED] Pipeline stopped at: {step_name}")
            print("   Please check the error message above.")
            sys.exit(1)
        except Exception as e:
            print(f"\n❌ [ERROR] Unexpected error at {step_name}: {e}")
            sys.exit(1)

    total_elapsed = time.time() - total_start_time
    print(f"\n==================================================")
    print(f"   PIPELINE COMPLETE SUCCESSFULLY!")
    print(f"   Total Time: {total_elapsed:.2f}s")
    print(f"   Final Output: ../new_output/final_genotypes_Ploidy{ploidy_number}_{mode}.txt")
    print(f"==================================================")

if __name__ == "__main__":
    run_pipeline()