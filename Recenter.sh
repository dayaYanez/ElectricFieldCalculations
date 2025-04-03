#!/bin/bash

set -euo pipefail

START=0
END=100
GMX_CMD="gmx_mpi"
TPR_FILE="mdout.tpr"
TRJ_IN="mdout.xtc"
TRJ_OUT="traj_centered.xtc"

run_trjconv() {
    local folder="$1"
    local work_dir="${folder}"

    if [[ ! -d "$work_dir" ]]; then
        printf "Skipping folder %s: Directory not found.\n" "$work_dir" >&2
        return
    fi

    if [[ ! -f "$work_dir/$TPR_FILE" || ! -f "$work_dir/$TRJ_IN" ]]; then
        printf "Skipping folder %s: Required trajectory or TPR files missing.\n" "$work_dir" >&2
        return
    fi

    (
        cd "$work_dir" || exit 1
        echo -e "0\n2" | $GMX_CMD trjconv -s "$TPR_FILE" -f "$TRJ_IN" -o "$TRJ_OUT" -center -pbc mol
    )
}

main() {
    for ((i=START; i<=END; i++)); do
        run_trjconv "$i"
    done
}

main

