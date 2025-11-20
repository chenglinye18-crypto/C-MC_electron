#!/usr/bin/env bash
# run_postprocess.sh
# Helper to reproduce the usual plotting commands after a Monte Carlo run.
# Usage (from finfet/):
#   bash scripts/run_postprocess.sh 6999
# Optional environment variables:
#   MC_ENV   : conda environment name (default mc_viz)
#   PLANE    : 2D slice plane for potential/heat/real-space plots (default xy)
#   INDEX    : slice index along the orthogonal axis (default 15)
#   GRID     : mesh file (default lgrid.txt)
#   CURRENT_CONTACT : icont index for plot_current.py (default 0)
#   CURRENT_START   : first step to include (default 0)
#   CURRENT_YMIN/YMAX: optional y-axis bounds (blank = auto)

set -euo pipefail

STEP="${1:-6999}"
MC_ENV="${MC_ENV:-mc_viz}"
PLANE="${PLANE:-xy}"
INDEX="${INDEX:-15}"
GRID="${GRID:-lgrid.txt}"
CURRENT_CONTACT="${CURRENT_CONTACT:-0}"
CURRENT_START="${CURRENT_START:-0}"
CURRENT_FILE="${CURRENT_FILE:-data/current}"
CURRENT_YMIN="${CURRENT_YMIN:-}"
CURRENT_YMAX="${CURRENT_YMAX:-}"

cd "$(dirname "$0")/.."

run_py() {
  conda run -n "$MC_ENV" python "$@"
}

echo "Using step=$STEP, plane=$PLANE, index=$INDEX, grid=$GRID"

# Potential (2D & 3D)
run_py scripts/plot_2D_pot.py --pot_file "data/pot${STEP}" \
  --plane "$PLANE" --index "$INDEX" --grid "$GRID" \
  --output "pot_${PLANE}_${INDEX}.png"

run_py scripts/plot_3D_pot.py --pot_file "data/pot${STEP}" \
  --grid "$GRID" --output "pot${STEP}.vtr"

# Heat map (2D & 3D)
run_py scripts/plot_2D_heat.py --heat_file "data/heat${STEP}" \
  --plane "$PLANE" --index "$INDEX" --grid "$GRID" \
  --output "heat_${PLANE}_${INDEX}.png"

run_py scripts/plot_3D_heat.py --heat_file "data/heat${STEP}" \
  --grid "$GRID" --output "heat${STEP}.vtr"

# Real-space carrier density (electrons & holes)
run_py scripts/plot_2D_realspace.py --data_file "data/Electron${STEP}" \
  --plane "$PLANE" --index "$INDEX" --grid "$GRID" \
  --output "electron_${PLANE}_${INDEX}.png" --log

run_py scripts/plot_3D_realspace.py --data_file "data/Electron${STEP}" \
  --grid "$GRID" --output "Electron${STEP}_realspace.vtr"

run_py scripts/plot_2D_realspace.py --data_file "data/Hole${STEP}" \
  --plane "$PLANE" --index "$INDEX" --grid "$GRID" \
  --output "hole_${PLANE}_${INDEX}.png" --log

run_py scripts/plot_3D_realspace.py --data_file "data/Hole${STEP}" \
  --grid "$GRID" --output "Hole${STEP}_realspace.vtr"

# k-space scatter (electrons)
run_py scripts/plot_kspace.py "data/Electron${STEP}" \
  --plane xy --output "kspace_xy_${STEP}.png" \
  --vtk "Electron${STEP}_kspace.vtp"

# Drain/source current plots
current_args=(scripts/plot_current.py --current_file "$CURRENT_FILE"
  --contact "$CURRENT_CONTACT" --start_step "$CURRENT_START"
  --output "current_contact${CURRENT_CONTACT}.png")
[[ -n "$CURRENT_YMIN" ]] && current_args+=(--ymin="$CURRENT_YMIN")
[[ -n "$CURRENT_YMAX" ]] && current_args+=(--ymax="$CURRENT_YMAX")
run_py "${current_args[@]}"

echo "Post-processing completed."
