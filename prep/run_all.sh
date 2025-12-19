#!/bin/bash

# Configuration
IMAGE="go2432/charmm:latest"
CURRENT_DIR=$(pwd)

# --- Check for Docker Image ---
if [[ "$(docker images -q $IMAGE 2> /dev/null)" == "" ]]; then
  echo "Image $IMAGE not found locally."
else
  echo "Image $IMAGE found locally. Proceeding..."
fi

# Helper function
run_in_docker() {
    echo "Running: $1"
    docker run --rm \
        -v "$CURRENT_DIR:/data" \
        -w /data \
        "$IMAGE" \
        /bin/bash -c "$1"
}

# --- Step 1: Pre-processing ---
run_in_docker "convpdb.pl -sel 1:5000 -out charmm22 -renumber 1 -segnames -cleanaux myc.pdb > myc_fixed.pdb"
run_in_docker "convpdb.pl -sel 1:5000 -out charmm22 -renumber 1 -segnames -cleanaux max.pdb > max_fixed.pdb"

# --- Step 2: Run CHARMM ---
run_in_docker "charmm < myc_step1_add_h.inp > log_myc"
run_in_docker "charmm < max_step1_add_h.inp > log_max"

# --- Step 3: Run Fix Scripts (Renames Residues, but leaves wrong H atoms) ---
run_in_docker "bash ./myc_step3_fix_charmm_pdb_for_amber.sh"
run_in_docker "bash ./max_step3_fix_charmm_pdb_for_amber.sh"

# --- Step 4: Combine PDBs (Aligns them) ---
run_in_docker "python3 myc_max_alignment_preproc.py"

# --- Step 4.5: FIX HISTIDINE CONFLICTS (Crucial Step) ---
# The previous scripts renamed the residues to HIE/HID but left the CHARMM 'HD1' atoms.
# We delete 'HD1' so TLeap doesn't crash. TLeap will rebuild the correct hydrogen.
run_in_docker "sed -i '/HD1/d' myc_max_amber_rotated.pdb"

# --- Step 5: Run TLeap ---
run_in_docker "tleap -f prepareSys_myc_max_ff19sb_opc.in"

echo "All steps complete."
