#!/usr/bin/env python3
"""
plot_3D_realspace.py
--------------------
Convert Electron/Hole data files to VTK RectilinearGrid storing charge
or charge density.
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


def load_volume_map(pvolume_file: Path, dims):
    nx, ny, nz = dims
    vol = np.zeros((nx, ny, nz))
    arr = np.loadtxt(pvolume_file)
    if arr.ndim == 1:
        arr = arr[None, :]
    for row in arr:
        i, j, k, value = row
        vol[int(i), int(j), int(k)] = value
    return vol


def accumulate_charge(data, dims, volume=None):
    nx, ny, nz = dims
    grid = np.zeros((nx, ny, nz))
    for row in data:
        i, j, k = map(int, row[:3])
        grid[i, j, k] += row[7]
    if volume is not None:
        mask = volume > 0
        grid[mask] = grid[mask] / volume[mask]
    return grid


def write_vtr(x, y, z, scalar, out_path: Path, name="charge"):
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
    parser = argparse.ArgumentParser(description="Convert Electron/Hole data to VTK rectilinear grid.")
    parser.add_argument("--data_file", type=Path, help="Path to data/ElectronXXXX or data/HoleXXXX")
    parser.add_argument("--output", type=Path, help="Output .vtr filename")
    parser.add_argument("--grid", type=Path, default=Path("./lgrid.txt"))
    parser.add_argument("--pvolume", type=Path, default=None, help="Optional pvolume file to convert to density")
    args = parser.parse_args()

    x, y, z = read_grid(args.grid)
    dims = (len(x), len(y), len(z))

    raw = np.loadtxt(args.data_file)
    if raw.ndim == 1:
        raw = raw[None, :]

    volume = load_volume_map(args.pvolume, dims) if args.pvolume else None
    scalar = accumulate_charge(raw, dims, volume=volume)
    name = "charge_density" if volume is not None else "charge"
    write_vtr(x, y, z, scalar, args.output, name=name)
    print(f"Saved {name} rectilinear grid to {args.output}")


if __name__ == "__main__":
    main()
