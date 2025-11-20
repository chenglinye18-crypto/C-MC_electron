#!/usr/bin/env python3
"""
plot_3D_heat.py
---------------
Convert heatXXXX files to VTK RectilinearGrid for visualising heat sources.
"""

import argparse
from pathlib import Path

import numpy as np

COLS = {
    "electron_per_charge": 3,
    "electron_density": 4,
    "hole_per_charge": 5,
    "hole_density": 6,
    "total_per_charge": None,
    "total_density": None,
}

UNITS = {
    "electron_per_charge": "W_per_electron",
    "electron_density": "W_per_cm3",
    "hole_per_charge": "W_per_hole",
    "hole_density": "W_per_cm3",
    "total_per_charge": "W_per_carrier",
    "total_density": "W_per_cm3",
}


def value_from_row(row, component):
    if component in ("electron_per_charge", "electron_density", "hole_per_charge", "hole_density"):
        return -row[COLS[component]]
    if component == "total_per_charge":
        return -(row[COLS["electron_per_charge"]] + row[COLS["hole_per_charge"]])
    if component == "total_density":
        return -(row[COLS["electron_density"]] + row[COLS["hole_density"]])
    raise ValueError(f"Unknown component {component}")


def read_grid(grid_path: Path):
    coords = []
    with grid_path.open("r", encoding="utf-8") as f:
        for _ in range(3):
            line = ""
            while line.strip() == "":
                line = f.readline()
                if not line:
                    raise ValueError("Unexpected end of lgrid.txt when reading counts")
            count = int(line.strip())
            values = []
            while len(values) < count:
                coord_line = f.readline()
                if not coord_line:
                    raise ValueError("Unexpected end of lgrid.txt when reading coordinates")
                stripped = coord_line.strip()
                if stripped == "":
                    continue
                values.append(float(stripped))
            coords.append(np.array(values))
    return coords


def accumulate_heat(files, dims, component, aggregate):
    nx, ny, nz = dims
    grid = np.zeros((nx, ny, nz))
    for path in files:
        raw = np.loadtxt(path)
        if raw.ndim == 1:
            raw = raw[None, :]
        for row in raw:
            i, j, k = map(int, row[:3])
            grid[i, j, k] += value_from_row(row, component)
    if aggregate == "mean" and len(files) > 0:
        grid /= len(files)
    return grid


def write_vtr(x, y, z, scalar, out_path: Path, name):
    nx, ny, nz = len(x), len(y), len(z)
    with out_path.open("w", encoding="utf-8") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write(f'  <RectilinearGrid WholeExtent="0 {nx-1} 0 {ny-1} 0 {nz-1}">\n')
        f.write(f'    <Piece Extent="0 {nx-1} 0 {ny-1} 0 {nz-1}">\n')
        f.write('      <PointData Scalars="{}">\n'.format(name))
        f.write('        <DataArray type="Float32" Name="{}" format="ascii">\n'.format(name))
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    f.write(f'          {scalar[i, j, k]:.8e}\n')
        f.write('        </DataArray>\n')
        f.write('      </PointData>\n')
        f.write('      <Coordinates>\n')
        for axis_name, coords in zip(
                ["X_COORDINATES", "Y_COORDINATES", "Z_COORDINATES"], [x, y, z]):
            f.write(f'        <DataArray type="Float32" Name="{axis_name}" NumberOfComponents="1" format="ascii">\n')
            for val in coords:
                f.write(f'          {val:.8e}\n')
            f.write('        </DataArray>\n')
        f.write('      </Coordinates>\n')
        f.write('    </Piece>\n')
        f.write('  </RectilinearGrid>\n')
        f.write('</VTKFile>\n')


def main():
    parser = argparse.ArgumentParser(description="Convert heatXXXX to VTK rectilinear grid.")
    parser.add_argument(
        "--heat_file",
        type=Path,
        nargs="+",
        dest="heat_files",
        required=True,
        help="One or more heat files (averaged if multiple).",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--grid", type=Path, default=Path("../lgrid.txt"))
    parser.add_argument("--component", choices=list(COLS.keys()), default="total_density")
    parser.add_argument("--aggregate", choices=["mean", "sum"], default="mean",
                        help="Combine multiple heat files via mean (default) or sum.")
    args = parser.parse_args()

    x, y, z = read_grid(args.grid)
    dims = (len(x), len(y), len(z))
    scalar = accumulate_heat(args.heat_files, dims, args.component, args.aggregate)
    name = f"{args.component}_{UNITS[args.component]}"
    write_vtr(x, y, z, scalar, args.output, name=name)
    print(f"Saved heat VTR to {args.output}")


if __name__ == "__main__":
    main()
