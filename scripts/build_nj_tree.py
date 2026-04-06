#!/usr/bin/env python3
"""
Build a Neighbor-Joining tree from a MUSCLE alignment.
Output: Newick (.nwk) file compatible with Jalview.

Usage:
    python build_nj_tree.py <alignment.afa> <output.nwk> [--model identity]
"""

import sys
import argparse
from pathlib import Path

from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import write as write_tree


def main():
    parser = argparse.ArgumentParser(
        description="Build a Neighbor-Joining tree from an aligned FASTA."
    )
    parser.add_argument("alignment", help="Input aligned FASTA (.afa)")
    parser.add_argument("output", help="Output Newick tree (.nwk)")
    parser.add_argument(
        "--model",
        default="blosum62",
        help="Distance model for DistanceCalculator (default: blosum62)",
    )
    args = parser.parse_args()

    # Read alignment
    aln = AlignIO.read(args.alignment, "fasta")
    print(f"[INFO] Loaded alignment: {len(aln)} sequences, {aln.get_alignment_length()} columns", file=sys.stderr)

    if len(aln) < 3:
        print("[WARN] Fewer than 3 sequences — NJ needs at least 3. Writing star tree.", file=sys.stderr)

    # Distance matrix
    calculator = DistanceCalculator(args.model)
    dm = calculator.get_distance(aln)

    # Neighbor-Joining
    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(dm)

    # Root at midpoint for a cleaner Jalview display
    nj_tree.root_at_midpoint()

    # Ladderise so the tree looks tidy in Jalview
    nj_tree.ladderize()

    # Write Newick
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    write_tree(nj_tree, args.output, "newick")

    print(f"[INFO] NJ tree written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
