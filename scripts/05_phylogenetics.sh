#!/usr/bin/env bash
# Phase 4 - phylogenetic reconstruction
#
# Two independent methods:
#   A. PHYLIP Neighbor-Joining with 1000 bootstrap replicates
#   B. IQ-TREE Maximum Likelihood (GTR+F+I+G4) with UFBoot + SH-aLRT
#
# PHYLIP is a 1980s program that reads from a file literally named
# "infile" and writes "outfile"/"outtree". Consense wants "intree" -
# yes, different from every other PHYLIP program. Just do the mv dance
# and move on.
#
# PAUP* was originally planned but is licensed-only, not available via
# conda. IQ-TREE is what modern papers actually use now anyway and
# runs in ~18 seconds for a 17-taxon dataset.

set -euo pipefail

PROJ="$HOME/erbb_project"

# ============================================================
# A. PHYLIP Neighbor-Joining
# ============================================================
mkdir -p "$PROJ/results/trees"
cd "$PROJ/results/trees"

# Classic PHYLIP truncates names to 10 chars and then chokes because
# "EGFR_chicken" -> "EGFR_chick" + leftover "en" reads as a bad base.
# Solution: write the alignment with short names (GgEGFR, DmEGFR, etc.)
echo "[$(date)] Writing short-name PHYLIP alignment..."
python3 << 'PYEOF'
from Bio import AlignIO
name_map = {
    "EGFR_human":"HsEGFR", "EGFR_mouse":"MmEGFR", "EGFR_chicken":"GgEGFR",
    "EGFR_zebrafish":"DrEGFR", "EGFR_drosophila":"DmEGFR",
    "ERBB2_human":"HsERBB2", "ERBB2_mouse":"MmERBB2", "ERBB2_chicken":"GgERBB2", "ERBB2_zebrafish":"DrERBB2",
    "ERBB3_human":"HsERBB3", "ERBB3_mouse":"MmERBB3", "ERBB3_chicken":"GgERBB3", "ERBB3_zebrafish":"DrERBB3",
    "ERBB4_human":"HsERBB4", "ERBB4_mouse":"MmERBB4", "ERBB4_chicken":"GgERBB4", "ERBB4_zebrafish":"DrERBB4",
}
aln = AlignIO.read("/home/preet/erbb_project/results/alignments/erbb_codon_aligned.fasta","fasta")
n_tax, n_char = len(aln), aln.get_alignment_length()
with open("infile.phy","w") as f:
    f.write(f"  {n_tax}  {n_char}\n")
    for rec in aln:
        f.write(f"{name_map.get(rec.id, rec.id[:10]):<10} {str(rec.seq).upper()}\n")
print(f"Wrote infile.phy: {n_tax} taxa x {n_char} sites")
PYEOF

# --- seqboot: 1000 bootstrap replicates ---
# menu order: R (change #replicates) -> 1000 -> Y (accept) -> 12345 (seed)
# seed must be ODD. Not kidding.
echo "[$(date)] seqboot..."
cp infile.phy infile
rm -f outfile
printf "R\n1000\nY\n12345\n" | seqboot > seqboot.log 2>&1
mv outfile infile.boot

# --- dnadist: pairwise distances (Kimura 2-parameter) ---
# menu order: D (change model) -> K (Kimura) -> M -> D (datasets) -> 1000 -> Y
echo "[$(date)] dnadist..."
cp infile.boot infile
rm -f outfile
printf "D\nK\nM\nD\n1000\nY\n" | dnadist > dnadist.log 2>&1
mv outfile infile.dist

# --- neighbor: NJ trees for each replicate ---
# outgroup is taxon 2 = DmEGFR (alphabetical order after chimp dropped)
# menu order: 2 (outgroup) -> M -> D (datasets) -> 1000 -> Y
echo "[$(date)] neighbor..."
cp infile.dist infile
rm -f outfile outtree
printf "2\nM\nD\n1000\nY\n" | neighbor > neighbor.log 2>&1
mv outfile infile.trees
mv outtree infile.treefile

# --- consense: majority-rule consensus ---
# reads from "intree" not "infile". yes this is different.
echo "[$(date)] consense..."
cp infile.treefile intree
rm -f outfile outtree
printf "Y\n" | consense > consense.log 2>&1
mv outtree erbb_consensus.treefile

echo "PHYLIP NJ tree: $(wc -c < erbb_consensus.treefile) bytes"

# ============================================================
# B. IQ-TREE Maximum Likelihood
# ============================================================
mkdir -p "$PROJ/results/trees_iqtree"
cd "$PROJ/results/trees_iqtree"

cp "$PROJ/results/alignments/erbb_codon_aligned.fasta" erbb_aln.fasta

# ModelFinder auto-selects the best-fit substitution model.
# -B 1000 = ultrafast bootstrap
# -alrt 1000 = SH-aLRT support (second metric)
# -o = outgroup (long name here, this is the fasta alignment)
echo "[$(date)] IQ-TREE..."
iqtree -s erbb_aln.fasta \
       -m MFP \
       -B 1000 \
       -alrt 1000 \
       -o EGFR_drosophila \
       -T 4 \
       --prefix erbb_ml

echo "[$(date)] IQ-TREE done."
grep "Best-fit model" erbb_ml.iqtree
