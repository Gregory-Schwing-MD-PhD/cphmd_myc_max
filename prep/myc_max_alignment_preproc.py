#!/usr/bin/env python
# coding: utf-8

# In[1]:


import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.tests.datafiles import CRD, PSF, DCD, DCD2
import nglview as nv


# In[4]:


myc_rcsb = mda.Universe("myc_amber.pdb")
nv.show_mdanalysis(myc_rcsb)


# In[5]:


max_rcsb = mda.Universe("max_amber.pdb")
nv.show_mdanalysis(max_rcsb)


# In[6]:


merged = mda.Merge(myc_rcsb.atoms, max_rcsb.atoms)
nv.show_mdanalysis(merged)

from MDAnalysis.coordinates.PDB import PDBWriter
# Define the output PDB file
output_pdb_file = "myc_max_amber.pdb"

# Create a PDBWriter object to write the PDB file
with PDBWriter(output_pdb_file) as writer:
    # Write only the first frame
    merged.trajectory[0]  # Access the first frame
    writer.write(merged.atoms)  # Write atoms to the PDB file

print(f"Frame written to {output_pdb_file}")


# In[28]:


import numpy as np
import mdtraj as md
from scipy.spatial.transform import Rotation as R

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


# In[29]:


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


# In[30]:


myc_max_amber_rotated = mda.Universe(output_pdb)
nv.show_mdanalysis(myc_max_amber_rotated)


# In[ ]:





# In[ ]:




