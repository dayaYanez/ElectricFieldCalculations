import os
import subprocess

phy_values = [-0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0,
              0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01]

def read_coordinates(filename):
    atoms = []
    coords = []
    with open(filename, 'r') as f:
        next(f)
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                atoms.append(parts[0])
                coords.append([float(x) for x in parts[1:4]])
    return atoms, coords

opt_atoms, R_opt = read_coordinates('R_opt.txt')
freq_atoms, R_freq = read_coordinates('R_freq.txt')

if opt_atoms != freq_atoms:
    raise ValueError("Atom lists in R_opt.txt and R_freq.txt don't match!")
if len(R_opt) != len(R_freq):
    raise ValueError("Number of coordinates in R_opt.txt and R_freq.txt don't match!")

template = """# wB97xD/Def2TZVP SCF(Tight) Int(Grid=299974)			

 Displacement_{phy}			

-1 1			
"""

def format_label(phy: float) -> str:
    return f"{phy:.2f}".replace('-', 'm').replace('.', '_')

for phy in phy_values:
    displaced_coords = [
        [
            opt[0] + phy * freq[0],
            opt[1] + phy * freq[1],
            opt[2] + phy * freq[2]
        ]
        for opt, freq in zip(R_opt, R_freq)
    ]

    label = format_label(phy)
    com_filename = f"R_dis_{label}.com"
    out_filename = f"R_dis_{label}.out"

    output_content = template.format(phy=phy)
    for atom, coord in zip(opt_atoms, displaced_coords):
        output_content += f"{atom}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
    output_content += "\n"  # Ensure empty line at the end

    with open(com_filename, 'w') as f:
        f.write(output_content)

    print(f"Generated {com_filename}")

    try:
        cmd = f"g16 < {com_filename} > {out_filename}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Ran Gaussian16: {com_filename} -> {out_filename}")
    except subprocess.CalledProcessError as e:
        print(f"g16 execution failed for {com_filename}: {e}")

print("All files processed and g16 jobs submitted.")


