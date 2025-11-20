#!/usr/bin/env python3
"""
plot_3D_pot.py
---------------
Convert a data/potXXXX file into a VTK rectilinear grid (.vtr) for 3D visualisation.

Usage:
  python3 plot_3D_pot.py ../data/pot6999 pot6999.vtr
"""

import argparse
from pathlib import Path

import numpy as np


def read_grid(grid_path: Path):
    coords = []
    with grid_path.open("r", encoding="utf-8") as f:
        for _ in range(3):
            line = ""
            while line.strip() == "":
                line = f.readline()
                if not line:
                    raise ValueError("Unexpected end of lgrid.txt when reading counts")
            num = int(line.strip())
            vals = []
            while len(vals) < num:
                coord_line = f.readline()
                if not coord_line:
                    raise ValueError("Unexpected end of lgrid.txt when reading coordinates")
                stripped = coord_line.strip()
                if stripped == "":
                    continue
                vals.append(float(stripped))
            coords.append(np.array(vals))
    return coords  # x, y, z


def load_potential(pot_file: Path, dims):
    nx, ny, nz = dims
    arr = np.loadtxt(pot_file)
    if arr.ndim == 1:
        arr = arr[None, :]
    grid = np.zeros((nx, ny, nz))
    for row in arr:
        i, j, k, val = row
        grid[int(i), int(j), int(k)] = val
    return grid


def write_vtr(x, y, z, potentials, out_path: Path):
    nx, ny, nz = len(x), len(y), len(z)
    with out_path.open("w", encoding="utf-8") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write(f'  <RectilinearGrid WholeExtent="0 {nx-1} 0 {ny-1} 0 {nz-1}">\n')
        f.write(f'    <Piece Extent="0 {nx-1} 0 {ny-1} 0 {nz-1}">\n')
        f.write('      <PointData Scalars="potential">\n')
        f.write('        <DataArray type="Float32" Name="potential" format="ascii">\n')
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    f.write(f'          {potentials[i, j, k]:.8e}\n')
        f.write('        </DataArray>\n')
        f.write('      </PointData>\n')
        f.write('      <Coordinates>\n')
        f.write('        <DataArray type="Float32" Name="X_COORDINATES" NumberOfComponents="1" format="ascii">\n')
        for val in x:
            f.write(f'          {val:.8e}\n')
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="Float32" Name="Y_COORDINATES" NumberOfComponents="1" format="ascii">\n')
        for val in y:
            f.write(f'          {val:.8e}\n')
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="Float32" Name="Z_COORDINATES" NumberOfComponents="1" format="ascii">\n')
        for val in z:
            f.write(f'          {val:.8e}\n')
        f.write('        </DataArray>\n')
        f.write('      </Coordinates>\n')
        f.write('    </Piece>\n')
        f.write('  </RectilinearGrid>\n')
        f.write('</VTKFile>\n')


def main():
    parser = argparse.ArgumentParser(description="Convert pot data to VTK rectilinear grid.")
    parser.add_argument("--pot_file", type=Path, help="Path to data/potXXXX file")
    parser.add_argument("--output", type=Path, help="Output .vtr filename")
    parser.add_argument("--grid", type=Path, default=Path("../lgrid.txt"), help="Path to lgrid.txt")
    args = parser.parse_args()

    x, y, z = read_grid(args.grid)
    potentials = load_potential(args.pot_file, (len(x), len(y), len(z)))
    write_vtr(x, y, z, potentials, args.output)
    print(f"Saved VTK rectilinear grid to {args.output}")


if __name__ == "__main__":
    main()
