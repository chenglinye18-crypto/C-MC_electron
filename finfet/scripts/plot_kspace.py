#!/usr/bin/env python3
"""
plot_kspace.py
---------------
Visualise particle distribution in k-space (using velocity components
as proxies for crystal momentum).

Example:
  python plot_kspace.py ../data/Electron6999 --output electron_kspace.png
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def label_valley(vx, vy, vz):
    vec = np.array([vx, vy, vz])
    axis = np.argmax(np.abs(vec))
    sign = 1 if vec[axis] >= 0 else 0
    return axis * 2 + sign


def main():
    parser = argparse.ArgumentParser(description="Plot k-space distribution from Electron/Hole files.")
    parser.add_argument("data_file", type=Path, help="Path to data/ElectronXXXX or HoleXXXX")
    parser.add_argument("--max_points", type=int, default=100000, help="Subsample for plotting speed")
    parser.add_argument("--output", type=Path, default=Path("kspace.png"))
    parser.add_argument("--vtk", type=Path, default=None, help="Optional PolyData (.vtp) output")
    parser.add_argument("--title", type=str, default=None)
    args = parser.parse_args()

    raw = np.loadtxt(args.data_file)
    if raw.ndim == 1:
        raw = raw[None, :]

    velocities = raw[:, 3:6]
    if args.max_points and velocities.shape[0] > args.max_points:
        idx = np.random.choice(velocities.shape[0], args.max_points, replace=False)
        velocities = velocities[idx]

    labels = np.array([label_valley(*v) for v in velocities])
    cmap = plt.get_cmap("tab10")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    planes = [("vx", "vy"), ("vx", "vz"), ("vy", "vz")]
    for ax, (a, b) in zip(axes, planes):
        idx_map = {"vx": 0, "vy": 1, "vz": 2}
        ax.scatter(velocities[:, idx_map[a]], velocities[:, idx_map[b]],
                   c=labels, cmap=cmap, s=5, alpha=0.6)
        ax.set_xlabel(f"{a}")
        ax.set_ylabel(f"{b}")
        ax.set_title(f"{a}-{b} projection")
        ax.grid(True, alpha=0.3)

    title = args.title or f"k-space distribution: {args.data_file.name}"
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(args.output, dpi=300)
    print(f"Saved k-space plot to {args.output}")

    if args.vtk:
        write_vtp(velocities, labels, args.vtk)
        print(f"Saved VTK PolyData to {args.vtk}")


def write_vtp(points, labels, path: Path):
    with path.open("w", encoding="utf-8") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <PolyData>\n')
        f.write(f'    <Piece NumberOfPoints="{points.shape[0]}" NumberOfVerts="{points.shape[0]}">\n')
        f.write('      <PointData Scalars="valley">\n')
        f.write('        <DataArray type="Int32" Name="valley" format="ascii">\n')
        for label in labels:
            f.write(f'          {label}\n')
        f.write('        </DataArray>\n')
        f.write('      </PointData>\n')
        f.write('      <Points>\n')
        f.write('        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        for vx, vy, vz in points:
            f.write(f'          {vx:.8e} {vy:.8e} {vz:.8e}\n')
        f.write('        </DataArray>\n')
        f.write('      </Points>\n')
        f.write('      <Verts>\n')
        f.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for idx in range(points.shape[0]):
            f.write(f'          {idx}\n')
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
        for idx in range(1, points.shape[0] + 1):
            f.write(f'          {idx}\n')
        f.write('        </DataArray>\n')
        f.write('      </Verts>\n')
        f.write('    </Piece>\n')
        f.write('  </PolyData>\n')
        f.write('</VTKFile>\n')


if __name__ == "__main__":
    main()
