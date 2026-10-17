#!/usr/bin/env python3
"""
Parse HyPhy MEME JSON output and map episodic sites to human ErbB genes.
"""
import json
import pandas as pd
from Bio import AlignIO
import sys
import re

MEME_JSON = "results/selection/meme_results.json"
ALIGNMENT = "results/asr/alignment.fasta"

print("Loading MEME results...")
with open(MEME_JSON) as f:
    data = json.load(f)

# Unwrap the data
content = data.get("MLE", {}).get("content", {})
if isinstance(content, dict) and "0" in content and len(content) == 1:
    content = content["0"]

# Determine structure (Row vs Column)
rows = []
if isinstance(content, dict):
    keys = list(content.keys())
    if keys[0].isdigit():
        rows = [content[k] for k in sorted(keys, key=int)]
    else:
        rows = list(zip(*[content[k] for k in keys]))
elif isinstance(content, list):
    if len(content) > 0 and isinstance(content[0], list):
        if len(content) == len(data["MLE"]["headers"]):
            rows = list(zip(*content))
        else:
            rows = content

# Headers
headers = [h[0].lower() if isinstance(h, list) else str(h).lower() for h in data["MLE"]["headers"]]

# Find columns
idx_episodic = next((i for i, h in enumerate(headers) if "episodic" in h or "p+" in h), None)
idx_pval = next((i for i, h in enumerate(headers) if "p-value" in h), None)

if idx_episodic is None:
    idx_episodic = 6 # Default based on MEME table layout
if idx_pval is None:
    idx_pval = 4 # Default based on MEME table layout

def get_float(val):
    while isinstance(val, (list, tuple)):
        if not val: return 0.0
        val = val[0]
    if isinstance(val, (int, float)): return float(val)
    match = re.search(r"[-+]?\d*\.\d+|\d+", str(val))
    if match: return float(match.group())
    return 0.0

# Map to human genes
aln = AlignIO.read(ALIGNMENT, "fasta")
human_genes = {
    "EGFR_human": "EGFR", "ERBB2_human": "ERBB2", "ERBB3_human": "ERBB3", "ERBB4_human": "ERBB4"
}
gene_maps = {}
for seq_id, name in human_genes.items():
    seq = next(rec.seq for rec in aln if rec.id == seq_id)
    mapping = {}
    res = 1
    for i in range(0, len(seq), 3):
        if "-" not in seq[i:i+3] and "N" not in seq[i:i+3]:
            mapping[i//3] = res
            res += 1
    gene_maps[name] = mapping

sites = []
for i, row in enumerate(rows):
    pval = get_float(row[idx_pval])
    episodic_val = str(row[idx_episodic]).lower()
    
    # Check if it's a 'Yes' or a numeric p-value <= 0.1
    is_episodic = "yes" in episodic_val or (pval <= 0.1 and pval > 0)
    
    if is_episodic:
        matches = [f"{g} Res {m[i]}" for g, m in gene_maps.items() if i in m]
        if matches:
            sites.append({
                "codon": i + 1,
                "p_value": pval,
                "human_residues": " | ".join(matches)
            })

print(f"\nFound {len(sites)} sites under episodic selection:\n")
for s in sites:
    print(f"  Codon {s['codon']:>4} (p={s['p_value']:.4f}) -> {s['human_residues']}")

pd.DataFrame(sites).to_csv("results/selection/episodic_selection_sites.tsv", sep="\t", index=False)
print("\nWrote results/selection/episodic_selection_sites.tsv")

