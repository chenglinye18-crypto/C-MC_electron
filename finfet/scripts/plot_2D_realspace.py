#!/usr/bin/env python3
"""
plot_2D_realspace.py
--------------------
Plot a 2D slice of the particle real-space charge (electron or hole)
based on data/ElectronXXXX or data/HoleXXXX files.

Example:
  python plot_2D_realspace.py ../data/Electron6999 --plane xy --index 10 \
      --grid ../lgrid.txt --pvolume ../data/pvolume --output electron_xy.png
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
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
    return coords  # x, y, z


def load_scalar_from_points(data, dims, column=7, volume=None):
    nx, ny, nz = dims
    grid = np.zeros((nx, ny, nz))
    for row in data:
        i, j, k = map(int, row[:3])
        grid[i, j, k] += row[column]
    if volume is not None:
        mask = volume > 0
        grid[mask] = grid[mask] / volume[mask]
    return grid


def load_volume_map(pvolume_file: Path, dims):
    nx, ny, nz = dims
    volume = np.zeros((nx, ny, nz))
    arr = np.loadtxt(pvolume_file)
    if arr.ndim == 1:
        arr = arr[None, :]
    for row in arr:
        i, j, k, val = row
        volume[int(i), int(j), int(k)] = val
    return volume


def extract_plane(grid, plane, index):
    if plane == "xy":
        return grid[:, :, index]
    if plane == "xz":
        return grid[:, index, :]
    if plane == "yz":
        return grid[index, :, :]
    raise ValueError("plane must be xy/xz/yz")


def main():
    parser = argparse.ArgumentParser(description="Plot 2D real-space slice from Electron/Hole data.")
    parser.add_argument("--data_file", type=Path, help="Path to data/ElectronXXXX or HoleXXXX")
    parser.add_argument("--plane", choices=["xy", "xz", "yz"], required=True)
    parser.add_argument("--index", type=int, required=True, help="Index along orthogonal axis")
    parser.add_argument("--grid", type=Path, default=Path("../lgrid.txt"))
    parser.add_argument("--pvolume", type=Path, default=None, help="Optional pvolume file to compute density")
    parser.add_argument("--output", type=Path, default=Path("realspace_slice.png"))
    parser.add_argument("--title", type=str, default=None)
    parser.add_argument(
        "--log",
        action="store_true",
        help="Use log10 color scale (values <= 0 clipped to eps).",
    )
    args = parser.parse_args()

    x, y, z = read_grid(args.grid)
    nx, ny, nz = len(x), len(y), len(z)

    raw = np.loadtxt(args.data_file)
    if raw.ndim == 1:
        raw = raw[None, :]

    volume = None
    if args.pvolume:
        volume = load_volume_map(args.pvolume, (nx, ny, nz))

    scalar_grid = load_scalar_from_points(raw, (nx, ny, nz), volume=volume)
    plane_array = extract_plane(scalar_grid, args.plane, args.index)

    if args.plane == "xy":
        X, Y = np.meshgrid(y, x)
        xlabel, ylabel = "Y (µm)", "X (µm)"
    elif args.plane == "xz":
        X, Y = np.meshgrid(z, x)
        xlabel, ylabel = "Z (µm)", "X (µm)"
    else:
        X, Y = np.meshgrid(z, y)
        xlabel, ylabel = "Z (µm)", "Y (µm)"

    plt.figure(figsize=(8, 5))
    plot_data = plane_array if args.plane != "yz" else plane_array.T
    if args.log:
        eps = np.max(plot_data[plot_data > 0]) * 1e-12 if np.any(plot_data > 0) else 1e-30
        plot_data = np.log10(np.maximum(plot_data, eps))
        cbar_label = "log10(Charge density)" if volume is not None else "log10(Charge)"
    else:
        cbar_label = "Charge density (arb. units)" if volume is not None else "Charge (arb. units)"

    im = plt.pcolormesh(X, Y, plot_data, shading="auto")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.colorbar(im, label=cbar_label)
    title = args.title or f"{args.data_file.name} {args.plane} idx={args.index}"
    plt.title(title)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved slice to {args.output}")


if __name__ == "__main__":
    main()
