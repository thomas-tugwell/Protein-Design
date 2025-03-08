# A script to dock the pore structures of two proteins
# Created by Thomas Tugwell, 3/8/2025

import copy
import numpy as np
from Bio import PDB

############################################
# Utility Functions
############################################

def compute_plane_center_normal(coords):
    """
    Given Nx3 coords for ring Cα atoms, returns:
      center (3-vector): mean of coords
      normal (3-vector): plane's normal via SVD.
    """
    center = np.mean(coords, axis=0)
    X = coords - center
    # SVD to find best-fit plane
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    normal = Vt[-1]
    normal /= np.linalg.norm(normal)
    return center, normal

def rotation_axis_angle(axis, angle):
    """
    Rodrigues rotation matrix for rotating by 'angle' radians about 'axis' (3D).
    """
    ax = axis / np.linalg.norm(axis)
    x, y, z = ax
    c = np.cos(angle)
    s = np.sin(angle)
    C = 1 - c
    R = np.array([
        [c + x*x*C,   x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s, c + y*y*C,   y*z*C - x*s],
        [z*x*C - y*s, z*y*C + x*s, c + z*z*C]
    ])
    return R

def rotation_from_vecs(v1, v2):
    """
    Returns a 3x3 rotation matrix that rotates v1 onto v2.
    If they're parallel or zero, handle special cases.
    """
    v1u = v1 / np.linalg.norm(v1)
    v2u = v2 / np.linalg.norm(v2)

    cross = np.cross(v1u, v2u)
    dot = np.dot(v1u, v2u)
    if np.linalg.norm(cross) < 1e-8:
        # parallel or opposite
        if dot < 0:
            # opposite
            perp = np.cross(v1u, [1,0,0])
            if np.linalg.norm(perp) < 1e-8:
                perp = np.cross(v1u, [0,1,0])
            perp /= np.linalg.norm(perp)
            return rotation_axis_angle(perp, np.pi)
        else:
            # same direction
            return np.eye(3)
    else:
        axis = cross / np.linalg.norm(cross)
        angle = np.arccos(dot)
        return rotation_axis_angle(axis, angle)

def apply_rotation_translation(structure, R, t):
    """
    x' = R*x + t for all atoms in 'structure'.
    """
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coord = atom.coord
                    atom.coord = R @ coord + t

############################################
# Main Script
############################################

# We do NOT change these filenames:
pdb_mobile  = "3_9_chainsfrom1_noend.pdb"
pdb_station = "pent_pore_chains_noend.pdb"

parser = PDB.PDBParser(QUIET=True)
structure1_orig = parser.get_structure("mobile",   pdb_mobile)
structure2_orig = parser.get_structure("station", pdb_station)

# Both are pentamers with ring-lining residues:
chains = ["A","B","C","D","E"]
res1 = 75   # ring residue in structure1
res2 = 109  # ring residue in structure2

# 1) Extract ring coords (Cα from each chain's residue)
coords1 = []
coords2 = []
for ch in chains:
    coords1.append(structure1_orig[0][ch][res1]["CA"].coord)
    coords2.append(structure2_orig[0][ch][res2]["CA"].coord)
coords1 = np.array(coords1)
coords2 = np.array(coords2)

# 2) Compute plane center & normal for each ring
c1_center, c1_normal = compute_plane_center_normal(coords1)
c2_center, c2_normal = compute_plane_center_normal(coords2)

# 3) Align structure1's ring plane to structure2's:
#    - rotate c1_normal -> c2_normal
#    - translate c1_center -> c2_center

R_plane = rotation_from_vecs(c1_normal, c2_normal)
t_plane = c2_center - R_plane @ c1_center

import copy
structure1_aligned = copy.deepcopy(structure1_orig)
apply_rotation_translation(structure1_aligned, R_plane, t_plane)

# 4) Generate 10 PDBs, each translating structure1_aligned along c2_normal
#    from +20 A to -5 A

n_steps = 10
start_offset = 20.0
end_offset   = -5.0
offset_values = np.linspace(start_offset, end_offset, n_steps)

from Bio import PDB
for i, offset in enumerate(offset_values):
    struct_mobile_copy = copy.deepcopy(structure1_aligned)
    shift_vec = offset * c2_normal
    apply_rotation_translation(struct_mobile_copy, np.eye(3), shift_vec)

    # Combine with structure2
    combined_struct = PDB.Structure.Structure("Combined")
    combined_model = PDB.Model.Model(0)
    combined_struct.add(combined_model)

    # add the shifted structure1
    for chain in struct_mobile_copy[0]:
        combined_model.add(chain)

    # rename structure2's chains A-E -> F-J
    structure2_copy = copy.deepcopy(structure2_orig)
    start_char = ord("F")
    for j, chain in enumerate(structure2_copy[0]):
        chain.id = chr(start_char + j)
        combined_model.add(chain)

    filename = f"ringplane_translation_{i}.pdb"
    io = PDB.PDBIO()
    io.set_structure(combined_struct)
    io.save(filename)
    print(f"Saved {filename} with offset={offset:.1f} along ring plane normal.")

print("\nDone! Created 10 PDBs translating structure1 along the ring-plane normal.")
