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
        "--contacts",
        type=int,
        nargs="+",
        help="Optional list of icont indices to plot together (overrides --contact).",
    )
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
    contacts = args.contacts if args.contacts else [args.contact]

    plt.figure(figsize=(7, 4))
    for ic in contacts:
        if ic not in data or len(data[ic]) == 0:
            print(f"[WARN] No entries for contact {ic} in {args.current_file}")
            continue
        filtered = [(s, v) for (s, v) in data[ic] if s >= args.start_step]
        if not filtered:
            print(f"[WARN] No entries for contact {ic} after start_step={args.start_step}")
            continue
        steps, values = zip(*filtered)
        plt.plot(steps, values, marker="o", lw=1, label=f"Contact {ic}")

    if not plt.gca().has_data():
        raise ValueError("No valid contacts to plot.")

    plt.xlabel("Step")
    plt.ylabel("Current (A)")
    if args.ymin is not None or args.ymax is not None:
        plt.ylim(args.ymin, args.ymax)
    title = args.title or "Contact current"
    plt.title(title)
    if len(contacts) > 1:
        plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved current plot to {args.output}")


if __name__ == "__main__":
    main()
