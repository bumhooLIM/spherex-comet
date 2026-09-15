import sys
import time
import subprocess
import pandas as pd
from pathlib import Path

# Assuming your directory configuration module is accessible
import directory 

def main():
    print("========== BATCH PHOTOMETRY ORCHESTRATOR ==========")
    
    DB_DIR = directory.DB_DIR
    OUT_DIR = directory.APPHOT_DIR
    fpath_db = DB_DIR / "db_filtered.parq"
    
    # 1. Extract Unique Targets Efficiently
    print(f"Reading database to extract unique targets...")
    # By specifying columns=['objdesig'], Pandas only loads that single column into RAM, 
    # making this read operation nearly instantaneous.
    df_targets = pd.read_parquet(fpath_db, columns=['objdesig'])
    unique_targets = df_targets['objdesig'].unique()
    
    total_targets = len(unique_targets)
    print(f"Found {total_targets} unique Solar System Objects to process.\n")
    
    start_time = time.time()
    
    # 2. Iterate and Execute
    for i, target in enumerate(unique_targets, start=1):
        print(f"\n" + "="*60)
        print(f"[{i}/{total_targets}] Initiating Pipeline for Target: {target}")
        print("="*60)
        
        # FAILSAFE: Check if this target is already finished. 
        # If the target's CSV file exists, we skip it.
        expected_output = OUT_DIR / f"{target}.csv"
        if expected_output.exists():
            print(f"  -> [SKIPPING] Output file '{target}.csv' already exists.")
            continue
            
        # 3. Build and Run the Command
        # sys.executable guarantees it uses the exact same Python environment you are currently in.
        cmd = [sys.executable, "apphot.py", "--objdesig", str(target)]
        
        try:
            # subprocess.run will execute apphot.py and stream its prints directly to your console
            subprocess.run(cmd, check=True)
            
        except subprocess.CalledProcessError as e:
            # If apphot.py crashes for this specific comet (e.g., missing FITS file), 
            # we catch the error, log it, and gracefully move on to the next comet.
            print(f"\n[ERROR] Pipeline failed for {target}!")
            
            # Log the failure so you can investigate later
            with open("failed_targets.log", "a") as f:
                f.write(f"{target}\n")
            print(f"  -> Target '{target}' logged to failed_targets.log. Moving to next target...")

    # 4. Final Summary
    elapsed = time.time() - start_time
    hours, rem = divmod(elapsed, 3600)
    minutes, seconds = divmod(rem, 60)
    print("\n" + "="*60)
    print(f"BATCH PROCESSING COMPLETE.")
    print(f"Total orchestration time: {int(hours)}h {int(minutes)}m {int(seconds)}s.")
    print("="*60)

if __name__ == "__main__":
    main()