#!/usr/bin/env python3
"""
plot_2D_heat.py
---------------
Visualise a 2D slice of heat data (electron or hole) stored in data/heatXXXX.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
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
    "electron_per_charge": "W / electron",
    "electron_density": "W / cm³",
    "hole_per_charge": "W / hole",
    "hole_density": "W / cm³",
    "total_per_charge": "W / carrier",
    "total_density": "W / cm³",
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
            num = int(line.strip())
            values = []
            while len(values) < num:
                coord_line = f.readline()
                if not coord_line:
                    raise ValueError("Unexpected end of lgrid.txt when reading coordinates")
                stripped = coord_line.strip()
                if stripped == "":
                    continue
                values.append(float(stripped))
            coords.append(np.array(values))
    return coords


def reshape_heat(files, dims, component, aggregate):
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


def extract_plane(grid, plane, index):
    if plane == "xy":
        return grid[:, :, index]
    if plane == "xz":
        return grid[:, index, :]
    if plane == "yz":
        return grid[index, :, :]
    raise ValueError("plane must be xy/xz/yz")


def main():
    parser = argparse.ArgumentParser(description="Plot 2D slice from heatXXXX files.")
    parser.add_argument(
        "--heat_file",
        type=Path,
        nargs="+",
        dest="heat_files",
        required=True,
        help="One or more heat files (averaged if multiple).",
    )
    parser.add_argument("--plane", choices=["xy", "xz", "yz"], required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--grid", type=Path, default=Path("../lgrid.txt"))
    parser.add_argument("--component", choices=list(COLS.keys()), default="total_density")
    parser.add_argument("--aggregate", choices=["mean", "sum"], default="mean",
                        help="Combine multiple heat files via mean (default) or sum.")
    parser.add_argument("--cmap", type=str, default="magma", help="Matplotlib colormap (default: magma)")
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument("--output", type=Path, default=Path("heat_slice.png"))
    parser.add_argument("--title", type=str, default=None)
    args = parser.parse_args()

    x, y, z = read_grid(args.grid)
    dims = (len(x), len(y), len(z))
    scalar_grid = reshape_heat(args.heat_files, dims, args.component, args.aggregate)
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
    vmin = args.vmin if args.vmin is not None else plot_data.min()
    vmax = args.vmax if args.vmax is not None else plot_data.max()
    norm = None
    if args.vmin is None and args.vmax is None and vmin < 0 < vmax:
        norm = TwoSlopeNorm(vcenter=0, vmin=vmin, vmax=vmax)

    im = plt.pcolormesh(
        X,
        Y,
        plot_data,
        shading="auto",
        cmap=args.cmap,
        vmin=None if norm else vmin,
        vmax=None if norm else vmax,
        norm=norm,
    )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    label = f"{args.component.replace('_', ' ').title()} ({UNITS[args.component]})"
    plt.colorbar(im, label=label)
    title_suffix = "mean" if args.aggregate == "mean" else "sum"
    plt.title(args.title or f"{args.component} ({title_suffix})")
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved heat slice to {args.output}")


if __name__ == "__main__":
    main()
