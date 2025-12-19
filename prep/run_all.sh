#!/bin/bash

# Configuration
IMAGE="go2432/charmm:latest"
CURRENT_DIR=$(pwd)

# --- Check for Docker Image ---
if [[ "$(docker images -q $IMAGE 2> /dev/null)" == "" ]]; then
  echo "Image $IMAGE not found locally. Pulling now..."
  docker pull $IMAGE
else
  echo "Image $IMAGE found locally. Proceeding..."
fi

# Helper function to run commands inside the Docker container
run_in_docker() {
    echo "Running: $1"
    docker run --rm \
        -v "$CURRENT_DIR:/data" \
        -w /data \
        "$IMAGE" \
        /bin/bash -c "$1"
}

# --- Step 1: Pre-processing with convpdb.pl ---
run_in_docker "convpdb.pl -sel 1:5000 -out charmm22 -renumber 1 -segnames -cleanaux myc.pdb > myc_fixed.pdb"
run_in_docker "convpdb.pl -sel 1:5000 -out charmm22 -renumber 1 -segnames -cleanaux max.pdb > max_fixed.pdb"

# --- Step 2: Run CHARMM ---
run_in_docker "charmm < myc_step1_add_h.inp > log_myc"
run_in_docker "charmm < max_step1_add_h.inp > log_max"

# --- Step 3: Run Fix Scripts ---
run_in_docker "bash ./myc_step3_fix_charmm_pdb_for_amber.sh"
run_in_docker "bash ./max_step3_fix_charmm_pdb_for_amber.sh"

# --- Step 4: Combine PDBs, align ---
run_in_docker "python myc_max_alignment_preproc.py"

echo "All steps complete."
