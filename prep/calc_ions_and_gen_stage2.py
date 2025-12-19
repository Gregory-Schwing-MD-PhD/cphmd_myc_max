import re
import sys

# --- Configuration ---
LOG_FILE = 'myc_max_ph75.log'
SOLVATED_PDB = 'temp_solvated.pdb'
OUTPUT_TLEAP_SCRIPT = 'stage2_ionize.in'
MOLARITY = 0.15
WATER_MOLARITY = 55.5

# 1. Get Protein Charge from Log
charge = 0.0
try:
    with open(LOG_FILE, 'r') as f:
        content = f.read()
        match = re.search(r'Total charge:\s+(-?\d+\.?\d*)', content)
        if match:
            charge = float(match.group(1))
            print(f"  > Found Protein Charge: {charge}")
        else:
            print("  > WARNING: Could not find 'Total charge' in log. Assuming 0.0")
except FileNotFoundError:
    print(f"  > ERROR: {LOG_FILE} not found. Did pdb4amber run?")
    sys.exit(1)

# 2. Count Waters from Stage 1 PDB
num_waters = 0
try:
    with open(SOLVATED_PDB, 'r') as f:
        for line in f:
            # We count Oxygen atoms belonging to WAT (or OPC/TIP3)
            if line.startswith('ATOM') and ('O   WAT' in line or 'O   OPC' in line): 
                num_waters += 1
    print(f"  > Detected Water Molecules: {num_waters}")
except FileNotFoundError:
    print(f"  > ERROR: {SOLVATED_PDB} not found. Did Stage 1 TLeap run?")
    sys.exit(1)

# 3. Calculate 0.15M Salt Ions
# Formula: N_pairs = (M_salt / M_water) * N_water
salt_pairs = int(round((MOLARITY / WATER_MOLARITY) * num_waters))
print(f"  > Salt Ions required for {MOLARITY}M: {salt_pairs} pairs")

# 4. Calculate Neutralization Ions
na_neutral = 0
cl_neutral = 0
rounded_charge = int(round(charge))

if rounded_charge < 0:
    na_neutral = abs(rounded_charge) # Add Na+ to neutralize neg charge
elif rounded_charge > 0:
    cl_neutral = rounded_charge      # Add Cl- to neutralize pos charge

# 5. Sum Total Ions
total_na = salt_pairs + na_neutral
total_cl = salt_pairs + cl_neutral

print(f"  > Total Na+ to add: {total_na} ({salt_pairs} salt + {na_neutral} neut)")
print(f"  > Total Cl- to add: {total_cl} ({salt_pairs} salt + {cl_neutral} neut)")

# 6. Write Final TLeap Script (Stage 2)
tleap_content = f"""# stage2_ionize.in - Generated automatically
source leaprc.protein.ff19SB
source leaprc.water.opc
set default PBradii mbondi

loadamberparams frcmod.phmd
loadoff phmd.lib

# Load the ALREADY SOLVATED structure from Stage 1
model = loadPdb {SOLVATED_PDB}

# Add Ions (0.15M + Neutralization)
# We use addIonsRand to replace random water molecules with ions
addIonsRand model Na+ {total_na} Cl- {total_cl}

# Check final charge
charge model

# Save Final Outputs
savePdb model myc_max_cphmd.pdb
saveAmberParm model myc_max_cphmd.prmtop myc_max_cphmd.rst7
quit
"""

with open(OUTPUT_TLEAP_SCRIPT, 'w') as f:
    f.write(tleap_content)

print(f"  > Successfully generated {OUTPUT_TLEAP_SCRIPT}")
