import numpy as np
import pandas as pd
from pathlib import Path
from MDAnalysis.analysis import distances

DISTANCE_CUTOFF = 37.794519772  # 20A in Bohr
O_CHARGE = -1.04
HW_CHARGE = 0.52


def load_coordinates(file_path):
    cols = ["resid", "resname", "atomname", "atom", "timestep", "x", "y", "z"]
    df = pd.read_csv(file_path, delim_whitespace=True, names=cols, skiprows=1)
    return df

def unit_vector(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else np.zeros_like(v)


def compute_electric_field(ni_pos, water_coords):
    ei_vector = np.zeros(3)
    for positions, charge in water_coords:
        dists = distances.distance_array(ni_pos.reshape(1, 3), positions)[0] 
        for idx, dist in enumerate(dists):
            if dist < DISTANCE_CUTOFF:
                r_vec = positions[idx] - ni_pos
                r_hat = r_vec / dist
                ei_vector += (charge / dist**2) * r_hat
    return ei_vector

# Compute bisector N2 bisector with N3, N4
def compute_bisector(n2, n3, n4):
    v1 = unit_vector(n3 - n2)
    v2 = unit_vector(n4 - n2)
    bisector = unit_vector(v1 + v2)
    return bisector

def main():
    df = load_coordinates("bohr_coordinates.txt")

    ni_atoms = ["N4", "N3", "N2"]
    timesteps = df["timestep"].unique()

    results = []

    for ts in timesteps:
        frame = df[df["timestep"] == ts]
        ni_df = frame[frame["atomname"].isin(ni_atoms)]

        ow_positions = frame[frame["atomname"] == "OW"][["x", "y", "z"]].values
        hw1_positions = frame[frame["atomname"] == "HW1"][["x", "y", "z"]].values
        hw2_positions = frame[frame["atomname"] == "HW2"][["x", "y", "z"]].values

        water_coords = [
            (ow_positions, O_CHARGE),
            (hw1_positions, HW_CHARGE),
            (hw2_positions, HW_CHARGE),
        ]

        for _, ni_row in ni_df.iterrows():
            ni_pos = np.array([ni_row["x"], ni_row["y"], ni_row["z"]])

            ei = compute_electric_field(ni_pos, water_coords)

            if ni_row["atomname"] == "N2":  # only compute projection once
                try:
                    n2 = ni_pos
                    n3 = frame[frame["atomname"] == "N3"][['x', 'y', 'z']].values[0]
                    n4 = frame[frame["atomname"] == "N4"][['x', 'y', 'z']].values[0]
                    rb = compute_bisector(n2, n3, n4)
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

