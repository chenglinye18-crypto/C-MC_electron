#!/usr/bin/env python3
"""
plot_current.py
---------------
Parse data/current and plot Current vs step for a chosen contact.

Usage:
  python plot_current.py ../data/current --contact 1 --output Ids.png
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt


def parse_current_file(path: Path):
    steps = []
    currents = {}
    current_step = None
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("step"):
                parts = line.split("=")
                current_step = int(parts[1])
                steps.append(current_step)
            elif line.startswith("icont"):
                icont = int(line.split("=")[1])
                currents.setdefault(icont, [])
            elif line.startswith("Current"):
                value = float(line.split("=")[1])
                currents.setdefault(icont, [])
                currents[icont].append((current_step, value))
    return currents


def main():
    parser = argparse.ArgumentParser(description="Plot contact current vs step from data/current.")
    parser.add_argument("--current_file", type=Path)
    parser.add_argument("--contact", type=int, default=1, help="Contact index (icont)")
    parser.add_argument(
        "--start_step",
        type=int,
        default=0,
        help="Ignore records with step < start_step (skip transient phase).",
    )
    parser.add_argument("--output", type=Path, default=Path("current.png"))
    parser.add_argument("--title", type=str, default=None)
    parser.add_argument("--ymin", type=float, default=None, help="Lower bound for the Y axis")
    parser.add_argument("--ymax", type=float, default=None, help="Upper bound for the Y axis")
    args = parser.parse_args()

    data = parse_current_file(args.current_file)
    if args.contact not in data or len(data[args.contact]) == 0:
        raise ValueError(f"No entries for contact {args.contact} in {args.current_file}")

    filtered = [(s, v) for (s, v) in data[args.contact] if s >= args.start_step]
    if not filtered:
        raise ValueError(f"No entries for contact {args.contact} after start_step={args.start_step}")

    steps, values = zip(*filtered)

    plt.figure(figsize=(7, 4))
    plt.plot(steps, values, marker="o", lw=1)
    plt.xlabel("Step")
    plt.ylabel("Current (A)")
    if args.ymin is not None or args.ymax is not None:
        plt.ylim(args.ymin, args.ymax)
    title = args.title or f"Contact {args.contact} current"
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved current plot to {args.output}")


if __name__ == "__main__":
    main()
