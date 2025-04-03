import MDAnalysis as mda
import numpy as np

# Conversion factor
angstrom_to_bohr = 1.0 / 0.529177

# Load the topology and trajectory
u = mda.Universe("mdout.gro", "traj_centered.xtc")

# Open file for writing
with open("bohr_coordinates.txt", "w") as f:
    # Write header
    f.write("resid resname atomname atom timestep x y z\n")

    for ts in u.trajectory:
        timestep = ts.frame
        coords_bohr = ts.positions * angstrom_to_bohr

        for atom, pos in zip(u.atoms, coords_bohr):
            x, y, z = pos
            f.write(f"{atom.resid:4d} {atom.resname:>6s} {atom.name:>6s} {atom.index+1:5d} {timestep:5d} {x:.6f} {y:.6f} {z:.6f}\n")

print(" Saved with atom identity and Bohr coordinates in 'bohrIndxTimestep.txt'")

