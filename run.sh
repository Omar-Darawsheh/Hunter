#!/bin/bash
set -euo pipefail

if [ -z "${1:-}" ]; then
    echo "Usage: ./run.sh /path/to/pdbs"
    exit 1
fi

snakemake --config pdb_dir="$1" --use-conda -j "$(nproc)"
