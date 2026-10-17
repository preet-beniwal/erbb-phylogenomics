#!/usr/bin/env python3
"""
Parse HyPhy FEL JSON. Handles nested wrappers, column/row orientation, 
and alpha=0 edge cases.
"""
import json
import pandas as pd
from Bio import AlignIO
import sys
import re

FEL_JSON = "results/selection/fel_results.json"
ALIGNMENT = "results/asr/alignment.fasta"

print("Loading FEL results...")
with open(FEL_JSON) as f:
    data = json.load(f)

content = data["MLE"]["content"]

# 1. Unwrap the "0" container if it exists
if isinstance(content, dict) and "0" in content and len(content) == 1:
    content = content["0"]

# 2. Determine if it's column-oriented or row-oriented
rows = []
if isinstance(content, dict):
    # Dictionary of rows
    rows = [content[k] for k in sorted(content.keys(), key=lambda x: int(x) if x.isdigit() else 0)]
elif isinstance(content, list):
    if len(content) == 6:  # There are 6 headers, so this must be a list of columns
        # Zip columns into rows: [(col1[0], col2[0], ...), (col1[1], col2[1], ...)]
        rows = list(zip(*content))
    else:
        # Already a list of rows
        rows = content

# 3. Helper to safely extract floats from nested lists or strings
def get_float(val):
    while isinstance(val, (list, tuple)):
        if len(val) == 0: return 0.0
        val = val[0]
    if isinstance(val, (int, float)): return float(val)
    # Extract numbers from strings like "Neg. p = 0.0038"
    match = re.search(r"[-+]?\d*\.\d+|\d+", str(val))
    if match:
        return float(match.group())
    return 0.0

# 4. Map alignment codons to human EGFR residues
aln = AlignIO.read(ALIGNMENT, "fasta")
human_seq = next(rec.seq for rec in aln if rec.id == "EGFR_human")

codon_to_res = {}
human_res_id = 1
for i in range(0, len(human_seq), 3):
    codon = human_seq[i:i+3]
    if "-" not in codon and "N" not in codon:
        codon_to_res[i // 3] = human_res_id
        human_res_id += 1

positive_sites = []

# 5. Iterate through rows
for i, row in enumerate(rows):
    # Based on the headers, indices are:
    # 0: alpha, 1: beta, 2: alpha=beta, 3: LRT, 4: p-value, 5: branch length
    alpha = get_float(row[0])
    beta  = get_float(row[1])
    pval  = get_float(row[4])
    
    # FEL positive selection: beta > alpha and p-value <= 0.1
    if beta > alpha and pval <= 0.1:
        human_res = codon_to_res.get(i, "unknown")
        positive_sites.append({
            "alignment_codon": i + 1,
            "human_residue": human_res,
            "alpha": alpha,
            "beta": beta,
            "p_value": pval
        })

# Report
print(f"\nFound {len(positive_sites)} sites under positive selection:\n")
for s in positive_sites:
    print(f"  Codon {s['alignment_codon']:>4} -> Human EGFR Residue {s['human_residue']:<4} (alpha={s['alpha']:.3f}, beta={s['beta']:.3f}, p={s['p_value']})")

# Save
df = pd.DataFrame(positive_sites)
out_path = "results/selection/positive_selection_sites.tsv"
df.to_csv(out_path, sep="\t", index=False)
print(f"\nWrote {out_path}")
