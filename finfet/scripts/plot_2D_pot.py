#!/usr/bin/env python3
"""
plot_2D_pot.py
---------------
Visualise a 2D slice of the potential stored in data/potXXXX files.

Usage example:
  python3 plot_2D_pot.py ../data/pot6999 --plane xy --index 10 --output pot_xy_k10.png
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_grid(grid_path: Path):
    """Return (x, y, z) coordinate arrays from lgrid.txt."""
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


def load_potential(pot_file: Path):
    data = np.loadtxt(pot_file)
    if data.ndim == 1:
        data = data[None, :]
    return data.astype(float)


def extract_plane(data, plane: str, index: int, dims):
    nx, ny, nz = dims
    grid = np.full(dims, np.nan)
    for row in data:
        i, j, k, val = row
        grid[int(i), int(j), int(k)] = val

    if plane == "xy":
        if not (0 <= index < nz):
            raise ValueError(f"k index {index} out of range (0..{nz-1})")
        plane_data = grid[:, :, index]
    elif plane == "xz":
        if not (0 <= index < ny):
            raise ValueError(f"j index {index} out of range (0..{ny-1})")
        plane_data = grid[:, index, :]
    elif plane == "yz":
        if not (0 <= index < nx):
            raise ValueError(f"i index {index} out of range (0..{nx-1})")
        plane_data = grid[index, :, :]
    else:
        raise ValueError("plane must be one of xy, xz, yz")
    return plane_data


def main():
    parser = argparse.ArgumentParser(description="Plot a 2D slice of potential data.")
    parser.add_argument("pot_file", type=Path, help="Path to data/potXXXX file")
    parser.add_argument("--grid", type=Path, default=Path("../lgrid.txt"), help="Path to lgrid.txt")
    parser.add_argument("--plane", choices=["xy", "xz", "yz"], required=True, help="Plane to visualise")
    parser.add_argument("--index", type=int, required=True, help="Index along the orthogonal axis (i/j/k)")
    parser.add_argument("--output", type=Path, default=Path("pot_slice.png"), help="Output image file")
    parser.add_argument("--title", type=str, default=None, help="Optional plot title")
    args = parser.parse_args()

    x, y, z = read_grid(args.grid)
    nx, ny, nz = len(x), len(y), len(z)
    pot_data = load_potential(args.pot_file)

    plane_array = extract_plane(pot_data, args.plane, args.index, (nx, ny, nz))

    if args.plane == "xy":
        X, Y = np.meshgrid(y, x)
        coord_label = ("Y (µm)", "X (µm)")
    elif args.plane == "xz":
        X, Y = np.meshgrid(z, x)
        coord_label = ("Z (µm)", "X (µm)")
    else:  # yz
        X, Y = np.meshgrid(z, y)
        coord_label = ("Z (µm)", "Y (µm)")

    plt.figure(figsize=(8, 5))
    cmap = plt.get_cmap("viridis")
    pcm = plt.pcolormesh(X, Y, plane_array if args.plane != "yz" else plane_array.T, shading="auto", cmap=cmap)
    plt.xlabel(coord_label[0])
    plt.ylabel(coord_label[1])
    title = args.title or f"{args.pot_file.name} - {args.plane}-plane idx={args.index}"
    plt.title(title)
    plt.colorbar(pcm, label="Potential (arb. units)")
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved plot to {args.output}")


if __name__ == "__main__":
    main()
