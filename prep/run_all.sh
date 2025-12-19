#!/bin/bash

# Configuration
IMAGE="go2432/charmm:latest"
CURRENT_DIR=$(pwd)

if [[ "$(docker images -q $IMAGE 2> /dev/null)" == "" ]]; then
  echo "Image $IMAGE not found."
  exit 1
fi

run_in_docker() {
    echo "----------------------------------------------------"
    echo "Running: $1"
    echo "----------------------------------------------------"
    docker run --rm -v "$CURRENT_DIR:/data" -w /data "$IMAGE" /bin/bash -c "$1"
}

# --- Steps 1-4: Pre-processing ---
run_in_docker "convpdb.pl -sel 1:5000 -out charmm22 -renumber 1 -segnames -cleanaux myc.pdb > myc_fixed.pdb"
run_in_docker "convpdb.pl -sel 1:5000 -out charmm22 -renumber 1 -segnames -cleanaux max.pdb > max_fixed.pdb"
run_in_docker "charmm < myc_step1_add_h.inp > log_myc"
run_in_docker "charmm < max_step1_add_h.inp > log_max"
run_in_docker "bash ./myc_step3_fix_charmm_pdb_for_amber.sh"
run_in_docker "bash ./max_step3_fix_charmm_pdb_for_amber.sh"
run_in_docker "python3 myc_max_alignment_preproc.py"

# --- Step 5: Manual Fix (Deleting bad Hydrogens) ---
echo "Removing conflicting HD1 atoms manually..."
run_in_docker "sed -i '/HD1/d' myc_max_amber_rotated.pdb"

# --- Step 6: STAGE 1 - Solvate ---
echo "Stage 1: Solvating..."
run_in_docker "tleap -f stage1_solvate.in"

# --- Step 7: STAGE 2 - Calculate Ions & Gen Script ---
echo "Calculating ions (assuming charge +10.23)..."
run_in_docker "python3 calc_ions_hardcode_charge_and_gen_stage2.py"

# --- Step 8: Run Final TLeap ---
echo "Stage 2: Adding ions and saving final parameters..."
run_in_docker "tleap -f stage2_ionize.in"

echo "All steps complete."
