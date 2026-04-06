#!/usr/bin/env python3
"""
Render a Newick tree as a PNG image using BioPython + matplotlib.

Usage:
    python render_tree_png.py <input.nwk> <output.png> [--width 12] [--height 8] [--dpi 300]
"""

import sys
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo


def main():
    parser = argparse.ArgumentParser(
        description="Render a Newick tree as PNG."
    )
    parser.add_argument("tree", help="Input Newick tree (.nwk)")
    parser.add_argument("output", help="Output PNG file")
    parser.add_argument("--width", type=float, default=12, help="Figure width in inches (default: 12)")
    parser.add_argument("--height", type=float, default=8, help="Figure height in inches (default: 8)")
    parser.add_argument("--dpi", type=int, default=300, help="Resolution (default: 300)")
    args = parser.parse_args()

    tree = Phylo.read(args.tree, "newick")

    n_tips = len(tree.get_terminals())
    height = max(args.height, n_tips * 0.3)

    fig, ax = plt.subplots(figsize=(args.width, height))
    Phylo.draw(tree, axes=ax, do_show=False)
    ax.set_title(f"Neighbor-Joining Tree ({n_tips} sequences)", fontsize=14)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"[INFO] Tree PNG written to {args.output} ({args.dpi} dpi)", file=sys.stderr)


if __name__ == "__main__":
    main()
