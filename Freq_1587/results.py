import os
import pandas as pd
import matplotlib.pyplot as plt

phy_values = [
    -0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0,
     0.1,  0.09,  0.08,  0.07,  0.06,  0.05,  0.04,  0.03,  0.02,  0.01
]

def format_label(phy: float) -> str:
    return f"{phy:.2f}".replace('-', 'm').replace('.', '_')

results = []

for phy in phy_values:
    filename = f"R_dis_{format_label(phy)}.out"
    if not os.path.exists(filename):
        continue

    with open(filename, 'r') as file:
        lines = file.readlines()

    E = None
    Dipole_x = Dipole_y = Dipole_z = Dipole_tot = None

    for i, line in enumerate(lines):
        if "SCF Done:" in line:
            parts = line.split()
            E = float(parts[4])
        elif "Dipole moment (field-independent basis, Debye):" in line:
            if i + 1 < len(lines):
                dipole_line = lines[i + 1].split()
                Dipole_x = float(dipole_line[1])
                Dipole_y = float(dipole_line[3])
                Dipole_z = float(dipole_line[5])
                Dipole_tot = float(dipole_line[7])

    results.append({
        "phy": phy,
        "E": E,
        "Dipole_x": Dipole_x,
        "Dipole_y": Dipole_y,
        "Dipole_z": Dipole_z,
        "Dipole_tot": Dipole_tot
    })

# Convert to DataFrame and save
output_df = pd.DataFrame(results)
output_df.to_csv("results.txt", sep='\t', index=False)
print(output_df)

# Plot phy vs E 
plt.figure()
plt.scatter(output_df['phy'], output_df['E'])
plt.xlabel('phy')
plt.ylabel('Energy (E)')
plt.title('phy vs E')
plt.grid(True)
plt.savefig("phy_vs_E.png")
plt.show()

# Plot phy vs Dipole_x
plt.figure()
plt.scatter(output_df['phy'], output_df['Dipole_x'])
plt.xlabel('phy')
plt.ylabel('Dipole_x')
plt.title('phy vs Dipole_x')
plt.grid(True)
plt.savefig("phy_vs_Dipole_x.png")
plt.show()

