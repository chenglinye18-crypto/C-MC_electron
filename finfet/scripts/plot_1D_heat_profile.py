#!/usr/bin/env python3
"""
plot_1D_heat_profile.py
-----------------------
Integrate the heatXXXX data along X/Z to obtain a 1D power profile vs Y.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_lgrid(path: Path):
    coords = []
    with path.open("r", encoding="utf-8") as f:
        for _ in range(3):
            line = ""
            while line.strip() == "":
                line = f.readline()
                if not line:
                    raise ValueError("Unexpected end of lgrid.txt")
            count = int(float(line.strip()))
            values = []
            while len(values) < count:
                line = f.readline()
                if not line:
                    raise ValueError("Unexpected end of lgrid.txt (coords)")
                stripped = line.strip()
                if stripped == "":
                    continue
                values.extend(float(x) for x in stripped.split())
            coords.append(np.array(values) * 1e-6)  # µm -> m
    return coords  # x, y, z


def load_heat_nodes(path: Path, nx, ny, nz, include, dt, stat):
    raw = np.loadtxt(path)
    if raw.ndim == 1:
        raw = raw[None, :]
    node_field = np.zeros((nx, ny, nz), dtype=float)
    scale = 1.602176634e-19 * 1e6 / (dt * stat)

    for row in raw:
        i, j, k = row[:3].astype(int)
        if not (0 <= i < nx and 0 <= j < ny and 0 <= k < nz):
            raise ValueError(f"Index ({i},{j},{k}) out of range")
        value = 0.0
        if include in ("both", "electron"):
            value += row[4]
        if include in ("both", "hole"):
            value += row[6]
        node_field[i, j, k] += -value * scale  # flip sign so + => lattice gain
    return node_field


def nodes_to_cells(field):
    return (
        field[:-1, :-1, :-1]
        + field[1:, :-1, :-1]
        + field[:-1, 1:, :-1]
        + field[:-1, :-1, 1:]
        + field[1:, 1:, :-1]
        + field[1:, :-1, 1:]
        + field[:-1, 1:, 1:]
        + field[1:, 1:, 1:]
    ) / 8.0


def compute_profile(
    cell_field, x_nodes, y_nodes, z_nodes, xrange_um, zrange_um, avg=True
):
    x_centers = 0.5 * (x_nodes[:-1] + x_nodes[1:])
    y_centers = 0.5 * (y_nodes[:-1] + y_nodes[1:])
    z_centers = 0.5 * (z_nodes[:-1] + z_nodes[1:])
    Xc, Yc, Zc = np.meshgrid(x_centers, y_centers, z_centers, indexing="ij")

    xrange_m = np.array(xrange_um) * 1e-6
    zrange_m = np.array(zrange_um) * 1e-6
    mask = (
        (Xc >= min(xrange_m))
        & (Xc <= max(xrange_m))
        & (Zc >= min(zrange_m))
        & (Zc <= max(zrange_m))
    )
    if not mask.any():
        raise ValueError(
            f"Selected X/Z range ({xrange_um} µm, {zrange_um} µm) contains no cells. "
            f"Domain X: [{x_nodes[0]*1e6:.3f}, {x_nodes[-1]*1e6:.3f}] µm, "
            f"Z: [{z_nodes[0]*1e6:.3f}, {z_nodes[-1]*1e6:.3f}] µm."
        )

    dx = np.diff(x_nodes)
    dy = np.diff(y_nodes)
    dz = np.diff(z_nodes)
    DX, DY, DZ = np.meshgrid(dx, dy, dz, indexing="ij")
    volumes = DX * DY * DZ

    power_density = cell_field * mask
    slice_power = np.sum(power_density * volumes, axis=(0, 2))
    if avg:
        slice_volume = np.sum(mask * volumes, axis=(0, 2))
        slice_density = np.divide(
            slice_power,
            slice_volume,
            out=np.zeros_like(slice_power),
            where=slice_volume > 0,
        )
        if not np.any(slice_volume):
            raise ValueError(
                "Requested X/Z range resulted in zero volume slices; "
                "check that the bounds overlap with the mesh."
            )
        return y_centers * 1e6, slice_density, slice_volume
    total_vol = np.sum(mask * volumes, axis=(0, 2))
    if not np.any(total_vol):
        raise ValueError(
            "Requested X/Z range resulted in zero volume; "
            "check that the bounds overlap with the mesh."
        )
    return y_centers * 1e6, slice_power, total_vol


def main():
    parser = argparse.ArgumentParser(
        description="Integrate heatXXXX along X/Z to get a 1D profile vs Y."
    )
    parser.add_argument(
        "--heat_file",
        type=Path,
        nargs="+",
        required=True,
        help="One or more heat files (averaged if multiple).",
    )
    parser.add_argument("--grid", type=Path, default=Path("lgrid.txt"))
    parser.add_argument("--dt", type=float, default=1e-16)
    parser.add_argument("--stat", type=float, default=5000)
    parser.add_argument(
        "--include", choices=["both", "electron", "hole"], default="both"
    )
    parser.add_argument(
        "--xrange", type=float, nargs=2, help="X range in µm (default: full domain)"
    )
    parser.add_argument(
        "--zrange", type=float, nargs=2, help="Z range in µm (default: full domain)"
    )
    parser.add_argument("--csv", type=Path, default=Path("heat_profile_y.csv"))
    parser.add_argument("--plot", type=Path, default=Path("heat_profile_y.png"))
    parser.add_argument(
        "--raw",
        action="store_true",
        help="Output total power (W) instead of average volumetric density.",
    )
    parser.add_argument("--title", type=str, default=None)
    args = parser.parse_args()

    x_nodes, y_nodes, z_nodes = read_lgrid(args.grid)
    nx, ny, nz = len(x_nodes), len(y_nodes), len(z_nodes)
    accum = np.zeros((nx, ny, nz), dtype=float)
    for hf in args.heat_file:
        accum += load_heat_nodes(hf, nx, ny, nz, args.include, args.dt, args.stat)
    heat_nodes = accum / len(args.heat_file)
    cell_field = nodes_to_cells(heat_nodes)

    xrange_um = (
        args.xrange
        if args.xrange
        else (x_nodes[0] * 1e6, x_nodes[-1] * 1e6)
    )
    zrange_um = (
        args.zrange
        if args.zrange
        else (z_nodes[0] * 1e6, z_nodes[-1] * 1e6)
    )
    y_um, profile, volumes = compute_profile(
        cell_field,
        x_nodes,
        y_nodes,
        z_nodes,
        xrange_um,
        zrange_um,
        avg=not args.raw,
    )

    header_unit = "W_per_m3" if not args.raw else "power_W"
    data = np.column_stack((y_um, profile))
    np.savetxt(
        args.csv,
        data,
        delimiter=",",
        header=f"y_um,{header_unit}",
        comments="",
    )

    plt.figure(figsize=(8, 4))
    plt.plot(y_um, profile, marker="o")
    plt.xlabel("Y (µm)")
    plt.title(
        args.title
        or f"Heat profile along Y (x={xrange_um}, z={zrange_um})"
    )
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.plot, dpi=300)
    print(f"Used {len(args.heat_file)} heat file(s)")
    print(f"Saved CSV to {args.csv}")
    print(f"Saved plot to {args.plot}")


if __name__ == "__main__":
    main()
