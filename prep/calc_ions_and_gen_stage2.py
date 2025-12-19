import re
import sys

# --- Configuration ---
TLEAP_LOG = 'stage1.log'           # Reading TLeap log for accurate charge
SOLVATED_PDB = 'temp_solvated.pdb'
OUTPUT_TLEAP_SCRIPT = 'stage2_ionize.in'
MOLARITY = 0.15
WATER_MOLARITY = 55.5

# 1. Get Protein Charge from TLeap Log
charge = 0.0
found_charge = False

try:
    with open(TLEAP_LOG, 'r') as f:
        # We look for: "Total unperturbed charge:   -4.000000"
        for line in f:
            if "Total unperturbed charge" in line or "Total perturbed charge" in line:
                match = re.search(r'charge:\s+(-?\d+\.\d+)', line)
                if match:
                    charge = float(match.group(1))
                    found_charge = True

    if found_charge:
        print(f"  > Found Protein Charge (from TLeap): {charge}")
    else:
        print(f"  > ERROR: Could not find 'Total unperturbed charge' in {TLEAP_LOG}")
        sys.exit(1)

except FileNotFoundError:
    print(f"  > ERROR: {TLEAP_LOG} not found. Did Stage 1 run?")
    sys.exit(1)

# 2. Count Waters from Stage 1 PDB
num_waters = 0
try:
    with open(SOLVATED_PDB, 'r') as f:
        for line in f:
            # Count residues named WAT/OPC/TIP3
            if line.startswith('ATOM') and ('O   WAT' in line or 'O   OPC' in line): 
                num_waters += 1
    print(f"  > Detected Water Molecules: {num_waters}")
except FileNotFoundError:
    print(f"  > ERROR: {SOLVATED_PDB} not found.")
    sys.exit(1)

# 3. Calculate 0.15M Salt Ions
salt_pairs = int(round((MOLARITY / WATER_MOLARITY) * num_waters))
print(f"  > Salt Ions required for {MOLARITY}M: {salt_pairs} pairs")

# 4. Calculate Neutralization Ions
na_neutral = 0
cl_neutral = 0
rounded_charge = int(round(charge))

if rounded_charge < 0:
    na_neutral = abs(rounded_charge)
elif rounded_charge > 0:
    cl_neutral = rounded_charge

# 5. Sum Total Ions
total_na = salt_pairs + na_neutral
total_cl = salt_pairs + cl_neutral

print(f"  > Total Na+ to add: {total_na} ({salt_pairs} salt + {na_neutral} neut)")
print(f"  > Total Cl- to add: {total_cl} ({salt_pairs} salt + {cl_neutral} neut)")

# 6. Write Final TLeap Script (Stage 2)
tleap_content = f"""# stage2_ionize.in
source leaprc.protein.ff19SB
source leaprc.water.opc
set default PBradii mbondi

loadamberparams frcmod.phmd
loadoff phmd.lib

# Load the ALREADY SOLVATED structure
model = loadPdb {SOLVATED_PDB}

# Add Ions (Replaces random waters)
addIonsRand model Na+ {total_na} Cl- {total_cl}

charge model

savePdb model myc_max_cphmd.pdb
saveAmberParm model myc_max_cphmd.prmtop myc_max_cphmd.rst7
quit
"""

with open(OUTPUT_TLEAP_SCRIPT, 'w') as f:
    f.write(tleap_content)

print(f"  > Successfully generated {OUTPUT_TLEAP_SCRIPT}")
