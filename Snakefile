import os
import glob
from pathlib import Path

configfile: "config/config.yaml"

# User must provide the directory containing PDB files
PDB_DIR = config.get("pdb_dir", None)
if PDB_DIR is None:
    raise ValueError(
        "ERROR: You must provide a PDB directory.\n"
        "Usage example:\n"
        "  snakemake --config pdb_dir=/path/to/pdbs --use-conda -j 4"
    )

PDB_DIR = os.path.abspath(PDB_DIR)

if not os.path.isdir(PDB_DIR):
    raise ValueError(f"ERROR: PDB directory does not exist:\n{PDB_DIR}")

PDBS = [
    os.path.splitext(os.path.basename(p))[0]
    for p in glob.glob(os.path.join(PDB_DIR, "*.pdb"))
]

if len(PDBS) == 0:
    raise ValueError(f"ERROR: No .pdb files found in directory: {PDB_DIR}")

for d in [
    "data/enzymm",
    "data/sequences",
    "logs/enzymm",
    "logs/extract_sequences",
    "data/alignments",
    "logs/muscle",
    "data/predictions",
    "logs/netsolp",
    "logs/phoptnn",
    "logs/seq2topt",
]:
    Path(d).mkdir(parents=True, exist_ok=True)

NETSOLP_DIR = "external/NetSolP-1.0"
NETSOLP_PS  = os.path.join(NETSOLP_DIR, "PredictionServer")
PHOPTNN_DIR = "external/pHoptNN"
SEQ2TOPT_DIR = "external/Seq2Topt"

# Read the list of PDBs that passed EnzyMM (have ≥1 hit)
def get_hit_pdbs(wildcards=None):
    """
    Called after the checkpoint resolves.  Returns list of PDB basenames
    whose EnzyMM TSV contained at least one data row.
    """
    ck = checkpoints.collect_enzymm_hits.get()
    with open(ck.output.hitlist) as fh:
        return [line.strip() for line in fh if line.strip()]


rule all:
    input:
        "data/alignments/muscle.afa",
        "data/trees/nj_tree.nwk",
        "data/trees/nj_tree.png",
        "data/predictions/netsolp.tsv",
        "data/predictions/phoptnn.tsv",
        "data/predictions/seq2topt_topt.tsv",
        "data/predictions/seq2topt_tm.tsv",
        "predictions_output/.done",
# Run EnzyMM on every PDB
rule run_enzymm:
    input:
        pdb=lambda wc: os.path.join(PDB_DIR, f"{wc.pdb}.pdb"),
    output:
        tsv="data/enzymm/{pdb}.tsv",
    log:
        "logs/enzymm/{pdb}.log",
    conda:
        "envs/enzymm.yaml"
    threads: 4
    shell:
        r"""
        enzymm \
            --input {input.pdb} \
            --output {output.tsv} \
            -j {threads} \
            --verbose > {log} 2>&1
        """


# Checkpoint: collect which PDBs had ≥1 EnzyMM hit
checkpoint collect_enzymm_hits:
    input:
        expand("data/enzymm/{pdb}.tsv", pdb=PDBS),
    output:
        hitlist="data/enzymm/hit_pdbs.txt",
    run:
        hits = []
        for tsv_path in input:
            pdb_name = os.path.splitext(os.path.basename(tsv_path))[0]
            has_hit = False
            with open(tsv_path) as fh:
                for line in fh:
                    # Skip comment lines and the header
                    if line.startswith("#"):
                        continue
                    if line.startswith("query_id"):
                        continue
                    # Any remaining non-empty line is a data row
                    if line.strip():
                        has_hit = True
                        break
            if has_hit:
                hits.append(pdb_name)
        hits.sort()
        with open(output.hitlist, "w") as out:
            for h in hits:
                out.write(h + "\n")


# Extract sequences — only for PDBs that passed EnzyMM
rule extract_sequences:
    input:
        pdb=lambda wc: os.path.join(PDB_DIR, f"{wc.pdb}.pdb"),
    output:
        fasta="data/sequences/{pdb}.fasta",
    log:
        "logs/extract_sequences/{pdb}.log",
    conda:
        "envs/extract_sequences.yaml"
    shell:
        r"""
        python scripts/extract_sequences.py {input.pdb} {output.fasta} 2> {log}
        """

rule merge_sequences:
    input:
        fastas=lambda wc: expand(
            "data/sequences/{pdb}.fasta",
            pdb=get_hit_pdbs(wc),
        ),
    output:
        merged="data/sequences/all_hits.fasta",
    log:
        "logs/merge_sequences.log",
    run:
        from Bio import SeqIO
        records = []
        for fasta in input.fastas:
            for rec in SeqIO.parse(fasta, "fasta"):
                records.append(rec)
        with open(output.merged, "w") as fh:
            SeqIO.write(records, fh, "fasta")
        with open(log[0], "w") as lf:
            lf.write(f"Merged {len(input.fastas)} files -> {len(records)} sequences\n")

rule setup_netsolp:
    output:
        sentinel=os.path.join(NETSOLP_PS, ".setup_done"),
    log:
        "logs/netsolp/setup_netsolp.log",
    shell:
        r"""
        bash scripts/setup_netsolp.sh {NETSOLP_DIR} > {log} 2>&1
        """
rule setup_phoptnn:
    output:
        sentinel=os.path.join(PHOPTNN_DIR, ".setup_done"),
    log:
        "logs/phoptnn/setup_phoptnn.log",
    shell:
        r"""
        bash scripts/setup_phoptnn.sh {PHOPTNN_DIR} > {log} 2>&1
        """

rule phoptnn:
    input:
        hitlist="data/enzymm/hit_pdbs.txt",
        setup=os.path.join(PHOPTNN_DIR, ".setup_done"),
    output:
        tsv="data/predictions/phoptnn.tsv",
    log:
        "logs/phoptnn/phoptnn.log",
    params:
        phoptnn_dir=PHOPTNN_DIR,
        pdb_dir=PDB_DIR,
    conda:
        "envs/phoptnn.yaml"
    shell:
        r"""
        python scripts/run_phoptnn.py \
            --pdb_dir {params.pdb_dir} \
            --hit_list {input.hitlist} \
            --phoptnn_dir {params.phoptnn_dir} \
            --output {output.tsv} \
            2> {log}
        """

rule setup_seq2topt:
    output:
        sentinel=os.path.join(SEQ2TOPT_DIR, ".setup_done"),
    log:
        "logs/seq2topt/setup_seq2topt.log",
    shell:
        r"""
        bash scripts/setup_seq2topt.sh {SEQ2TOPT_DIR} > {log} 2>&1
        """

rule seq2topt:
    input:
        fasta="data/sequences/all_hits.fasta",
        setup=os.path.join(SEQ2TOPT_DIR, ".setup_done"),
    output:
        topt="data/predictions/seq2topt_topt.tsv",
        tm="data/predictions/seq2topt_tm.tsv",
    log:
        "logs/seq2topt/seq2topt.log",
    params:
        seq2topt_dir=SEQ2TOPT_DIR,
    conda:
        "envs/seq2topt.yaml"
    shell:
        r"""
        python scripts/run_seq2topt.py \
            --fasta {input.fasta} \
            --seq2topt_dir {params.seq2topt_dir} \
            --output_topt {output.topt} \
            --output_tm {output.tm} \
            2> {log}
        """

rule netsolp:
    input:
        fasta="data/sequences/all_hits.fasta",
        setup=os.path.join(NETSOLP_PS, ".setup_done"),
    output:
        tsv="data/predictions/netsolp.tsv",
    log:
        "logs/netsolp/netsolp.log",
    params:
        netsolp_dir=NETSOLP_PS,
        model_type=config.get("netsolp_model", "ESM1b"),
    conda:
        "envs/netsolp.yaml"
    shell:
        r"""
        python scripts/run_netsolp.py \
            --fasta {input.fasta} \
            --output {output.tsv} \
            --netsolp_dir {params.netsolp_dir} \
            --model_type {params.model_type} \
            2> {log}
        """

rule muscle_alignment:
    input:
        "data/sequences/all_hits.fasta",
    output:
        "data/alignments/muscle.afa",
    log:
        "logs/muscle/muscle_alignment.log",
    threads: 4
    conda:
        "envs/muscle.yaml"
    shell:
        r"""
        muscle -align {input} -output {output} -threads {threads} > {log} 2>&1
        """

rule build_nj_tree:
    input:
        "data/alignments/muscle.afa",
    output:
        "data/trees/nj_tree.nwk",
    log:
        "logs/build_nj_tree.log",
    conda:
        "envs/extract_sequences.yaml"   # already has biopython
    shell:
        r"""
        python scripts/build_nj_tree.py {input} {output} 2> {log}
        """

rule tree_png:
    input:
        "data/trees/nj_tree.nwk",
    output:
        "data/trees/nj_tree.png",
    log:
        "logs/tree_png.log",
    conda:
        "envs/extract_sequences.yaml"
    shell:
        r"""
        python scripts/render_tree_png.py {input} {output} --dpi 300 2> {log}
        """

rule merge_predictions:
     input:
         hitlist="data/enzymm/hit_pdbs.txt",
         netsolp="data/predictions/netsolp.tsv",
         phoptnn="data/predictions/phoptnn.tsv",
         seq2topt="data/predictions/seq2topt_topt.tsv",
         seq2tm="data/predictions/seq2topt_tm.tsv",
     output:
         sentinel="predictions_output/.done",
     log:
         "logs/merge_predictions.log",
     conda:
         "envs/extract_sequences.yaml"
     shell:
         r"""
         python scripts/merge_predictions.py \
             --hit_list {input.hitlist} \
             --netsolp  {input.netsolp} \
             --phoptnn  {input.phoptnn} \
             --seq2topt {input.seq2topt} \
             --seq2tm   {input.seq2tm} \
             --output_dir predictions_output \
             2> {log} && touch {output.sentinel}
         """
