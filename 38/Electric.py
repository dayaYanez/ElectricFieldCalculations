import numpy as np
import pandas as pd
from pathlib import Path
from MDAnalysis.analysis import distances

DISTANCE_CUTOFF = 37.794519772  # 20Å in Bohr
MO_CHARGE = -1.04
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
            for pos, charge in zip([mw, hw1, hw2], [MO_CHARGE, HW_CHARGE, HW_CHARGE]):
                r_vec = pos - ni_pos
                r = np.linalg.norm(r_vec)
                if r != 0:
                    r_hat = r_vec / r
                    ei_vector += (charge / r**2) * r_hat
    return ei_vector


def compute_bisector(n2, n3, n4):
    v1 = unit_vector(n3 - n2)
    v2 = unit_vector(n4 - n2)
    return unit_vector(v1 + v2)


def main():
    df = load_coordinates("bohr_coordinates.txt")

    ni_atoms = ["N4", "N3", "N2"]
    timesteps = df["timestep"].unique()
    results = []

    for ts in timesteps:
        frame = df[df["timestep"] == ts]
        ni_df = frame[frame["atomname"].isin(ni_atoms)]

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

        for _, ni_row in ni_df.iterrows():
            ni_pos = np.array([ni_row["x"], ni_row["y"], ni_row["z"]])
            ei = compute_electric_field(ni_pos, water_molecules)

            if ni_row["atomname"] == "N2":
                try:
                    n3 = frame[frame["atomname"] == "N3"][['x', 'y', 'z']].values[0]
                    n4 = frame[frame["atomname"] == "N4"][['x', 'y', 'z']].values[0]
                    rb = compute_bisector(ni_pos, n3, n4)
                    e_out = np.dot(rb, ei)
                except IndexError:
                    e_out = np.nan
            else:
                e_out = np.nan

            results.append({
                "timestep": ts,
                "atom": ni_row["atomname"],
                "Ei_x": ei[0],
                "Ei_y": ei[1],
                "Ei_z": ei[2],
                "E_out": e_out
            })

    df_out = pd.DataFrame(results)
    df_out.to_csv("electric_field_output.txt", index=False, sep='\t')


if __name__ == "__main__":
    main()
