#!/usr/bin/env python3
"""
Extract ancestral states for the EGFR allosteric pocket (K745, F856, L858).
"""

from Bio import AlignIO
import pandas as pd

ALIGNMENT = "results/asr/alignment.fasta"
STATE_FILE = "results/asr/erbb_asr_nt.state"

# The residues we care about (human EGFR canonical numbering)
TARGET_RESIDUES = {
    745: "K (Lysine - catalytic)",
    856: "F (Phenylalanine - DFG motif)",
    858: "L (Leucine - allosteric pocket)"
}

# 1. Map human residue numbers to alignment columns
aln = AlignIO.read(ALIGNMENT, "fasta")
human_seq = None
for rec in aln:
    if rec.id == "EGFR_human":
        human_seq = rec.seq
        break

if human_seq is None:
    raise ValueError("EGFR_human not found in alignment")

res_to_cols = {}
res_idx = 1
nucl_idx = 0
for col_idx, base in enumerate(human_seq):
    if base != "-":
        nucl_idx += 1
        if (nucl_idx - 1) % 3 == 0:  # first position of a codon
            # Wait, if there are gaps, a codon might be split. 
            # Better logic: count non-gap nucleotides, residue = (nucl_idx // 3) + 1
            pass
        
        # Correct logic for codon alignment:
        # every 3 non-gap nucleotides = 1 residue
        current_res = (nucl_idx - 1) // 3 + 1
        if current_res in TARGET_RESIDUES:
            if current_res not in res_to_cols:
                res_to_cols[current_res] = []
            res_to_cols[current_res].append(col_idx)

print("Mapped human residues to alignment columns:")
for res, cols in res_to_cols.items():
    print(f"  Residue {res} ({TARGET_RESIDUES[res]}) -> columns {cols} (0-indexed)")

# 2. Load the IQ-TREE state file
print("\nLoading ASR states...")
df = pd.read_csv(STATE_FILE, sep="\t", comment="#")

# Find the root node. IQ-TREE usually names it 'Node1' or 'Root'
nodes = df['Node'].unique()
root_node = "Node1" if "Node1" in nodes else nodes[0]
print(f"Using root node: {root_node}")

# 3. Extract and translate
# Standard genetic code
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# Filter for the root node
root_df = df[df["Node"] == root_node].set_index("Site")

print("\nAncestral states at the vertebrate root:")
for res, cols in res_to_cols.items():
    # Columns are 0-indexed in python, but IQ-TREE sites are 1-indexed
    iqtree_sites = [c + 1 for c in cols]
    
    # Get the most likely state for each nucleotide in the codon
    codon = ""
    for site in iqtree_sites:
        if site in root_df.index:
            state = root_df.loc[site, "State"]
            codon += state
        else:
            codon += "?"
            
    aa = CODON_TABLE.get(codon.upper(), "?")
    print(f"  Residue {res} ({TARGET_RESIDUES[res]}):")
    print(f"    Ancestral Codon: {codon} -> Amino Acid: {aa}")
