#!/bin/bash

# Configuration
IMAGE="go2432/charmm:latest"
CURRENT_DIR=$(pwd)

# --- Check for Docker Image ---
if [[ "$(docker images -q $IMAGE 2> /dev/null)" == "" ]]; then
  echo "Image $IMAGE not found locally. Please build it first."
  exit 1
fi

# Helper function
run_in_docker() {
    echo "----------------------------------------------------"
    echo "Running: $1"
    echo "----------------------------------------------------"
    docker run --rm -v "$CURRENT_DIR:/data" -w /data "$IMAGE" /bin/bash -c "$1"
}

# --- Steps 1-4: Pre-processing (Same as before) ---
run_in_docker "convpdb.pl -sel 1:5000 -out charmm22 -renumber 1 -segnames -cleanaux myc.pdb > myc_fixed.pdb"
run_in_docker "convpdb.pl -sel 1:5000 -out charmm22 -renumber 1 -segnames -cleanaux max.pdb > max_fixed.pdb"
run_in_docker "charmm < myc_step1_add_h.inp > log_myc"
run_in_docker "charmm < max_step1_add_h.inp > log_max"
run_in_docker "bash ./myc_step3_fix_charmm_pdb_for_amber.sh"
run_in_docker "bash ./max_step3_fix_charmm_pdb_for_amber.sh"
run_in_docker "python3 myc_max_alignment_preproc.py"

# --- Step 5: Protonation & Hydrogen Optimization ---
# Output: myc_max_ph75.pdb + myc_max_ph75.log
echo "Running pdb4amber with --reduce optimization..."
run_in_docker "pdb4amber -i myc_max_amber_rotated.pdb -o myc_max_ph75.pdb --reduce -l myc_max_ph75.log"

# --- Step 6: STAGE 1 - Solvate ---
# Uses stage1_solvate.in (must be in current dir)
echo "Stage 1: Solvating to count waters..."
run_in_docker "tleap -f stage1_solvate.in"

# --- Step 7: STAGE 2 - Calculate Ions & Gen Script ---
# Uses calc_ions_and_gen_stage2.py (must be in current dir)
echo "Calculating ions and generating Stage 2 script..."
run_in_docker "python3 calc_ions_and_gen_stage2.py"

# --- Step 8: Run Final TLeap ---
# Runs the script generated in Step 7
echo "Stage 2: Adding ions and saving final parameters..."
run_in_docker "tleap -f stage2_ionize.in"

echo "All steps complete."
