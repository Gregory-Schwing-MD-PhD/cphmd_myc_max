import sys

# --- Configuration ---
SOLVATED_PDB = 'temp_solvated.pdb'
OUTPUT_TLEAP_SCRIPT = 'stage2_ionize.in'
MOLARITY = 0.15
WATER_MOLARITY = 55.5
PROTEIN_CHARGE = 10.23  # Hardcoded assumption

# 1. Count Waters from Stage 1 PDB
num_waters = 0
try:
    with open(SOLVATED_PDB, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and ('O   WAT' in line or 'O   OPC' in line): 
                num_waters += 1
    print(f"  > Detected Water Molecules: {num_waters}")
except FileNotFoundError:
    print(f"  > ERROR: {SOLVATED_PDB} not found.")
    sys.exit(1)

# 2. Calculate 0.15M Salt Ions
# Formula: N_pairs = (M_salt / M_water) * N_water
salt_pairs = int(round((MOLARITY / WATER_MOLARITY) * num_waters))
print(f"  > Salt Ions required for {MOLARITY}M: {salt_pairs} pairs")

# 3. Calculate Neutralization Ions
# Charge is +10.23 -> rounds to +10
# We need 10 Cl- ions to neutralize.
na_neutral = 0
cl_neutral = 0
rounded_charge = int(round(PROTEIN_CHARGE))

if rounded_charge < 0:
    na_neutral = abs(rounded_charge) # Add Na+ if negative
elif rounded_charge > 0:
    cl_neutral = rounded_charge      # Add Cl- if positive (this case)

# 4. Sum Total Ions
total_na = salt_pairs + na_neutral
total_cl = salt_pairs + cl_neutral

print(f"  > Protein Charge: {PROTEIN_CHARGE} (Rounds to {rounded_charge})")
print(f"  > Total Na+ to add: {total_na} ({salt_pairs} salt + {na_neutral} neut)")
print(f"  > Total Cl- to add: {total_cl} ({salt_pairs} salt + {cl_neutral} neut)")

# 5. Write Final TLeap Script (Stage 2)
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
