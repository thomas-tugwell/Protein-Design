# A script to dock the pore structures of two proteins
# Generated additional "flipped" orientations in case space_dock.py fails
# Created by Thomas Tugwell, 3/8/2025

import copy
import numpy as np
from Bio import PDB

############################################
# Utility Functions
############################################

def compute_plane_center_normal(coords):
    """Compute the best-fit plane and normal using SVD."""
    center = np.mean(coords, axis=0)
    X = coords - center
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    normal = Vt[-1] / np.linalg.norm(Vt[-1])
    return center, normal

def rotation_axis_angle(axis, angle):
    """Rodrigues rotation matrix for a given axis-angle."""
    ax = axis / np.linalg.norm(axis)
    x, y, z = ax
    c, s = np.cos(angle), np.sin(angle)
    C = 1 - c
    return np.array([
        [c + x*x*C, x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s, c + y*y*C, y*z*C - x*s],
        [z*x*C - y*s, z*y*C + x*s, c + z*z*C]
    ])

def rotation_from_vecs(v1, v2):
    """Compute a rotation matrix to align v1 to v2."""
    v1, v2 = v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2)
    cross = np.cross(v1, v2)
    dot = np.dot(v1, v2)
    if np.linalg.norm(cross) < 1e-8:
        return np.eye(3) if dot > 0 else rotation_axis_angle(np.cross(v1, [1,0,0]), np.pi)
    return rotation_axis_angle(cross / np.linalg.norm(cross), np.arccos(dot))

def apply_rotation_translation(structure, R, t):
    """Apply rigid transformation to all atoms in the structure."""
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.coord = R @ atom.coord + t

def flip_structure_across_plane(structure, center, normal):
    """Flips a structure across a given plane (defined by center and normal)."""
    flipped_structure = copy.deepcopy(structure)
    for model in flipped_structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coord = atom.coord
                    # Compute reflection
                    d = np.dot(coord - center, normal) * 2
                    reflected_coord = coord - d * normal
                    atom.coord = reflected_coord
    return flipped_structure

############################################
# Main Script
############################################

pdb_mobile  = "3_9_chainsfrom1_noend.pdb"
pdb_station = "pent_pore_chains_noend.pdb"

parser = PDB.PDBParser(QUIET=True)
structure1_orig = parser.get_structure("mobile", pdb_mobile)
structure2_orig = parser.get_structure("station", pdb_station)

# Define pentameric ring residues
chains = ["A", "B", "C", "D", "E"]
res1, res2 = 75, 109  # Mobile (res 75), Stationary (res 109)

# Extract ring Cα coordinates
coords1 = np.array([structure1_orig[0][ch][res1]["CA"].coord for ch in chains])
coords2 = np.array([structure2_orig[0][ch][res2]["CA"].coord for ch in chains])

# Compute plane centers & normals
c1_center, c1_normal = compute_plane_center_normal(coords1)
c2_center, c2_normal = compute_plane_center_normal(coords2)

# Align structure1’s ring plane to structure2’s ring plane
R_plane = rotation_from_vecs(c1_normal, c2_normal)
t_plane = c2_center - R_plane @ c1_center

# Apply transformation to align the mobile protein
structure1_aligned = copy.deepcopy(structure1_orig)
apply_rotation_translation(structure1_aligned, R_plane, t_plane)

# Create the flipped version of the aligned structure
structure1_flipped = flip_structure_across_plane(structure1_aligned, c2_center, c2_normal)

# Generate 10 translated PDBs for both orientations
n_steps = 10
offset_values = np.linspace(20.0, -5.0, n_steps)

from Bio import PDB
for i, offset in enumerate(offset_values):
    # Generate translations for normal-aligned structure
    struct_mobile_copy = copy.deepcopy(structure1_aligned)
    apply_rotation_translation(struct_mobile_copy, np.eye(3), offset * c2_normal)

    # Generate translations for flipped structure
    struct_mobile_flipped_copy = copy.deepcopy(structure1_flipped)
    apply_rotation_translation(struct_mobile_flipped_copy, np.eye(3), offset * c2_normal)

    # Combine structures for normal orientation
    combined_struct = PDB.Structure.Structure("Combined_Normal")
    combined_model = PDB.Model.Model(0)
    combined_struct.add(combined_model)

    # Add translated structure1 (mobile, chains A–E)
    for chain in struct_mobile_copy[0]:
        combined_model.add(chain)

    # Rename structure2's chains A-E → F-J
    structure2_copy = copy.deepcopy(structure2_orig)
    start_char = ord("F")
    for j, chain in enumerate(structure2_copy[0]):
        chain.id = chr(start_char + j)
        combined_model.add(chain)

    filename_normal = f"ringplane_translation_normal_{i}.pdb"
    io = PDB.PDBIO()
    io.set_structure(combined_struct)
    io.save(filename_normal)

    # Combine structures for flipped orientation
    combined_struct_flipped = PDB.Structure.Structure("Combined_Flipped")
    combined_model_flipped = PDB.Model.Model(0)
    combined_struct_flipped.add(combined_model_flipped)

    # Add translated flipped structure
    for chain in struct_mobile_flipped_copy[0]:
        combined_model_flipped.add(chain)

    # Add renamed stationary structure again
    structure2_copy_flipped = copy.deepcopy(structure2_orig)
    start_char = ord("F")
    for j, chain in enumerate(structure2_copy_flipped[0]):
        chain.id = chr(start_char + j)
        combined_model_flipped.add(chain)

    filename_flipped = f"ringplane_translation_flipped_{i}.pdb"
    io = PDB.PDBIO()
    io.set_structure(combined_struct_flipped)
    io.save(filename_flipped)

    print(f"Saved {filename_normal} and {filename_flipped} with offset={offset:.1f} Å along ring plane normal.")

print("\n✅ Done! Created 20 PDBs with:")
print("   - **10 PDBs for the normal orientation**")
print("   - **10 PDBs for the flipped orientation**")
print("   - **Both translated along the ring plane normal**")
