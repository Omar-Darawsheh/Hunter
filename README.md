# Hunter — Enzyme Property Prediction Pipeline

A Snakemake pipeline that takes PDB structure files and predicts key biochemical properties of enzymes:

- **Enzyme activity** — filters inputs using [EnzyMM](https://github.com/your-link) to identify likely enzymes
- **Solubility & usability** — predicted by [NetSolP-1.0](https://services.healthtech.dtu.dk/services/NetSolP-1.0/)
- **Optimal pH** — predicted by [pHoptNN](https://github.com/kuenzelab/pHoptNN)
- **Optimal temperature (Topt) & melting temperature (Tm)** — predicted by [Seq2Topt](https://github.com/SizheQiu/Seq2Topt)

The pipeline also generates a multiple sequence alignment (MUSCLE), a neighbor-joining phylogenetic tree, and a merged summary table of all predictions.

## Requirements

- Linux (tested on CentOS/RHEL)
- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/)
- [Snakemake](https://snakemake.readthedocs.io/) ≥ 7.0
- Internet access (for first run — downloads external tools and model weights automatically)

## Quick Start

```bash
git clone https://github.com/Omar-Darawsheh/Hunter.git
cd Hunter
chmod +x run.sh
./run.sh /path/to/your/pdb/files/
```

Or run Snakemake directly with more control:

```bash
snakemake --config pdb_dir=/path/to/your/pdb/files/ --use-conda -j 4
```

## Input

A directory containing one or more `.pdb` files. Each file should contain a protein structure.

## Output

All results are written to `predictions_output/`:

- **`all_predictions.tsv`** — one row per PDB with all predicted properties (if all inputs are single-chain)

If any input has multiple unique chains, two files are produced instead:

- **`all_predictions_structure.tsv`** — one row per PDB (pH at structure level; other properties per-chain marked as NA for multi-chain entries)
- **`all_predictions_chains.tsv`** — one row per unique chain of multi-chain PDBs

### Predicted Properties

| Column | Description | Tool |
|---|---|---|
| `predicted_solubility` | Predicted protein solubility | NetSolP |
| `predicted_usability` | Predicted protein usability | NetSolP |
| `predicted_ph_opt` | Predicted optimal pH | pHoptNN |
| `predicted_topt_C` | Predicted optimal temperature (°C) | Seq2Topt |
| `predicted_tm_C` | Predicted melting temperature (°C) | Seq2Topt |

### Additional Outputs

- `data/alignments/muscle.afa` — MUSCLE multiple sequence alignment
- `data/trees/nj_tree.nwk` — Neighbor-joining tree (Newick format, compatible with Jalview)
- `data/trees/nj_tree.png` — Tree rendered as PNG

## External Tools

On the first run, the pipeline automatically downloads and sets up:

- **NetSolP-1.0** (~5.6 GB) — DTU Health Tech
- **pHoptNN** — cloned from GitHub
- **Seq2Topt** — cloned from GitHub with model weights from GitHub Releases

These are installed into the `external/` directory and managed by Snakemake rules.

## Configuration

The pipeline is configured via the command line. Optional settings can be added to `config/config.yaml`:

```yaml
netsolp_model: "ESM1b"   # Options: ESM12, ESM1b, Distilled, Both
```

## License

MIT
