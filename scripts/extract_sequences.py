#!/usr/bin/env python3

import sys
from pathlib import Path
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def main():
    pdb_path = sys.argv[1]
    out_fasta = sys.argv[2]

    pdb_id = Path(pdb_path).stem

    # Group chains by sequence
    seq_to_chains = OrderedDict()
    for rec in SeqIO.parse(pdb_path, "pdb-seqres"):
        chain = str(rec.annotations["chain"]).strip()
        s = str(rec.seq)
        if s not in seq_to_chains:
            seq_to_chains[s] = []
        seq_to_chains[s].append(chain)

    if not seq_to_chains:
        print("[ERROR] No SEQRES sequences found in the PDB.", file=sys.stderr)
        sys.exit(1)

    records = []
    for i, (seq, chains) in enumerate(seq_to_chains.items(), 1):
        chain_str = ", ".join(chains)
        header = f"{pdb_id}_{i}|Chains {chain_str}|length={len(seq)}"
        rec = SeqRecord(Seq(seq), id=header, description="")
        records.append(rec)

    Path(Path(out_fasta).parent).mkdir(parents=True, exist_ok=True)

    with open(out_fasta, "w") as fh:
        SeqIO.write(records, fh, "fasta")

    print(f"[INFO] Extracted {len(records)} unique sequence(s) from {sum(len(c) for c in seq_to_chains.values())} chains.", file=sys.stderr)

if __name__ == "__main__":
    main()
