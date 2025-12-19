#!/usr/bin/env python
# coding: utf-8

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.coordinates.PDB import PDBWriter
import numpy as np
import mdtraj as md
from scipy.spatial.transform import Rotation as R
import warnings

# Suppress warnings to keep output clean
warnings.filterwarnings('ignore')

# ---------------------------------------------------------
# Part 1: Load, Renumber, Rename Chains/Segments, and Merge
# ---------------------------------------------------------

myc_rcsb = mda.Universe("myc_amber.pdb")
max_rcsb = mda.Universe("max_amber.pdb")

# 1. RENUMBER RESIDUES
# ---------------------------------------------------------
# Get the last residue number from Myc to calculate the offset
offset = myc_rcsb.residues.resids[-1]

# Apply offset to Max RESIDUES so numbering continues after Myc
max_rcsb.residues.resids += offset


# 2. ASSIGN SEGMENT IDs (MYC / MAX)
# ---------------------------------------------------------
# We set this on the .segments group to satisfy MDAnalysis hierarchy
myc_rcsb.atoms.segments.segids = ['MYC'] * len(myc_rcsb.atoms.segments)
max_rcsb.atoms.segments.segids = ['MAX'] * len(max_rcsb.atoms.segments)


# 3. ASSIGN CHAIN IDs (A / B)
# ---------------------------------------------------------
# We use add_TopologyAttr to safely force assignment to all atoms
myc_rcsb.add_TopologyAttr('chainID', ['A'] * len(myc_rcsb.atoms))
max_rcsb.add_TopologyAttr('chainID', ['B'] * len(max_rcsb.atoms))


# 4. MERGE AND WRITE
# ---------------------------------------------------------
merged = mda.Merge(myc_rcsb.atoms, max_rcsb.atoms)

output_pdb_file = "myc_max_amber.pdb"

with PDBWriter(output_pdb_file) as writer:
    merged.trajectory[0]  # Access the first frame
    writer.write(merged.atoms)

print(f"Merged frame written to {output_pdb_file}")
print(f"Myc: Chain A, Seg MYC (Residues 1-{offset})")
print(f"Max: Chain B, Seg MAX (Residues {offset+1}-{max_rcsb.residues.resids[-1]})")


# ---------------------------------------------------------
# Part 2: PCA Rotation (Unchanged)
# ---------------------------------------------------------

def load_pdb(file_path):
    """Load PDB file and extract atomic coordinates."""
    traj = md.load(file_path)
    return traj.xyz[0]  # Returns coordinates of the first frame

def compute_principal_axes(coords):
    """Compute principal axes using PCA."""
    mean = np.mean(coords, axis=0)
    centered_coords = coords - mean
    cov_matrix = np.cov(centered_coords.T)
    eigvals, eigvecs = np.linalg.eigh(cov_matrix)
    idx = np.argsort(eigvals)[::-1]
    eigvecs = eigvecs[:, idx]
    return eigvecs

def align_axis_to_z(coords, axis):
    """Align a given axis to the z-axis."""
    z_axis = np.array([0, 0, 1])
    rotation = R.align_vectors([z_axis], [axis])[0]
    rotated_coords = rotation.apply(coords)
    return rotated_coords

def rotate_by_90_degrees_y(coords):
    """Rotate coordinates by 90 degrees around the y-axis."""
    rotation = R.from_euler('x', -90, degrees=True)
    rotated_coords = rotation.apply(coords)
    return rotated_coords

def write_pdb(coords, original_pdb, output_pdb):
    """Write coordinates to a new PDB file using MDTraj."""
    traj = md.load(original_pdb)
    traj.xyz[0] = coords
    traj.save(output_pdb)

# Perform the rotation on the newly merged file
coords = load_pdb(output_pdb_file)
principal_axes = compute_principal_axes(coords)

# Use only the first principal axis
first_axis = principal_axes[:, 0]

# Align the first principal axis to the z-axis
aligned_coords = align_axis_to_z(coords, first_axis)

# Rotate by 90 degrees around the y-axis
final_coords = rotate_by_90_degrees_y(aligned_coords)

# Write to output PDB file
output_pdb = "myc_max_amber_rotated.pdb"
write_pdb(final_coords, output_pdb_file, output_pdb)
print(f"Written aligned and rotated structure: {output_pdb}")
