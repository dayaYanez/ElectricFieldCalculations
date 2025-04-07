import numpy as np
import pandas as pd
from pathlib import Path
from MDAnalysis.analysis import distances

DISTANCE_CUTOFF = 37.794519772  # 20Å in Bohr
O_CHARGE = -1.04
HW_CHARGE = 0.52


def load_coordinates(file_path):
    cols = ["resid", "resname", "atomname", "atom", "timestep", "x", "y", "z"]
    df = pd.read_csv(file_path, delim_whitespace=True, names=cols, skiprows=1)
    return df


def unit_vector(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else np.zeros_like(v)


def compute_electric_field(ni_pos, water_molecules):
    ei_vector = np.zeros(3)
    for ow, mw, hw1, hw2 in water_molecules:
        dist_ow = np.linalg.norm(ow - ni_pos)
        if dist_ow < DISTANCE_CUTOFF:
            for pos, charge in zip([mw, hw1, hw2], [O_CHARGE, HW_CHARGE, HW_CHARGE]):
                r_vec = pos - ni_pos
                r = np.linalg.norm(r_vec)
                if r != 0:
                    r_hat = r_vec / r
                    ei_vector += (charge / r**2) * r_hat
    return ei_vector


def compute_bisector(n_center, n1, n2):
    v1 = unit_vector(n1 - n_center)
    v2 = unit_vector(n2 - n_center)
    return unit_vector(v1 + v2)


def main():
    df = load_coordinates("bohr_coordinates.txt")

    timesteps = df["timestep"].unique()
    results = []

    for ts in timesteps:
        frame = df[df["timestep"] == ts]

        water_molecules = []
        for resid, group in frame.groupby("resid"):
            try:
                ow = group[group["atomname"] == "OW"][['x', 'y', 'z']].values[0]
                mw = group[group["atomname"] == "MW"][['x', 'y', 'z']].values[0]
                hw1 = group[group["atomname"] == "HW1"][['x', 'y', 'z']].values[0]
                hw2 = group[group["atomname"] == "HW2"][['x', 'y', 'z']].values[0]
                water_molecules.append((ow, mw, hw1, hw2))
            except IndexError:
                continue

        try:
            n2 = frame[frame["atomname"] == "N2"][['x', 'y', 'z']].values[0]
            n3 = frame[frame["atomname"] == "N3"][['x', 'y', 'z']].values[0]
            n4 = frame[frame["atomname"] == "N4"][['x', 'y', 'z']].values[0]
        except IndexError:
            continue

        ei_n2 = compute_electric_field(n2, water_molecules)
        ei_n3 = compute_electric_field(n3, water_molecules)
        ei_n4 = compute_electric_field(n4, water_molecules)

        bisector_n2 = compute_bisector(n2, n3, n4)
        bisector_n3 = compute_bisector(n3, n2, n4)
        bisector_n4 = compute_bisector(n4, n2, n3)

        proj_n2 = np.dot(bisector_n2, ei_n2)
        proj_n3 = np.dot(bisector_n3, ei_n3)
        proj_n4 = np.dot(bisector_n4, ei_n4)

        results.append({
            "timestep": ts,
            "N2": proj_n2,
            "N3": proj_n3,
            "N4": proj_n4
        })

    df_out = pd.DataFrame(results)
    df_out.to_csv("electric_field_output.txt", index=False, sep='\t')


if __name__ == "__main__":
    main()
