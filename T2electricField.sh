input_file="bohr_coordinates.txt"
output_file="electric_field.txt"
target_atoms="N4 C8 N3 N2 C7"
cutoff=37.794519772
charge=-1.04  # TIP4P charge on MW

awk -v cutoff="$cutoff" -v q="$charge" -v atoms="$target_atoms" -v out="$output_file" '
BEGIN {
    OFS=" "
    split(atoms, t_atoms, " ")
    print "timestep atomname Ei_x Ei_y Ei_z" > out
}
{
    if (NR == 1) next  # Skip header

    resid=$1
    resname=$2
    atomname=$3
    atom=$4
    timestep=$5
    x=$6
    y=$7
    z=$8

    if (atomname == "MW") {
        mw_coords[timestep][mw_count[timestep]++] = x " " y " " z
    }

    for (i in t_atoms) {
        if (atomname == t_atoms[i]) {
            idx = target_count[timestep]++
            target_names[timestep, idx] = atomname
            target_coords[timestep, idx] = x " " y " " z
        }
    }
}
END {
    for (t in target_count) {
        for (i = 0; i < target_count[t]; i++) {
            atomname = target_names[t, i]
            split(target_coords[t, i], Ni, " ")

            sum_ex = 0
            sum_ey = 0
            sum_ez = 0

            for (j = 0; j < mw_count[t]; j++) {
                split(mw_coords[t][j], MWj, " ")
                dx = Ni[1] - MWj[1]
                dy = Ni[2] - MWj[2]
                dz = Ni[3] - MWj[3]
                d = sqrt(dx*dx + dy*dy + dz*dz)

                if (d < cutoff) {
                    factor = q / (d*d)
                    sum_ex += factor * dx
                    sum_ey += factor * dy
                    sum_ez += factor * dz
                }
            }

            print t, atomname, sum_ex, sum_ey, sum_ez >> out
        }
    }
}
' "$input_file"

