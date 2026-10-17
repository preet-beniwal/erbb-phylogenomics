#!/usr/bin/env python3
"""
Map the 7 positively selected alignment codons to ALL human ErbB genes.
"""
from Bio import AlignIO
import pandas as pd

ALIGNMENT = "results/asr/alignment.fasta"

# The 7 positively selected alignment codons found by FEL
# (1-based indexing from the output)
TARGET_CODONS = [83, 103, 1401, 1435, 1506, 1553, 1823]

# Human genes in the alignment
HUMAN_GENES = {
    "EGFR_human": "EGFR",
    "ERBB2_human": "ERBB2 (HER2)",
    "ERBB3_human": "ERBB3 (HER3)",
    "ERBB4_human": "ERBB4 (HER4)"
}

aln = AlignIO.read(ALIGNMENT, "fasta")

# For each human gene, map alignment codon -> residue number
gene_maps = {}
for seq_id, gene_name in HUMAN_GENES.items():
    seq = next(rec.seq for rec in aln if rec.id == seq_id)
    mapping = {}
    res_id = 1
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if "-" not in codon and "N" not in codon:
            mapping[i // 3] = res_id
            res_id += 1
    gene_maps[gene_name] = mapping

print("Mapping the 7 positively selected sites to human ErbB genes:\n")

for target_1based in TARGET_CODONS:
    align_idx = target_1based - 1
    matches = []
    
    for gene_name, mapping in gene_maps.items():
        if align_idx in mapping:
            res = mapping[align_idx]
            matches.append(f"{gene_name} Residue {res}")
            
    if matches:
        print(f"Alignment Codon {target_1based:>4} -> " + " | ".join(matches))
    else:
        print(f"Alignment Codon {target_1based:>4} -> (Gap in all human sequences)")
