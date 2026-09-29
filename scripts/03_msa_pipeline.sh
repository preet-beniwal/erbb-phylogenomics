#!/usr/bin/env bash
# Phase 3 - cross-species MSA of the ErbB family
#
# took me a while to land on this approach:
# - full-length nucleotide alignment crashes (MUSCLE OOM'd at 24GB,
#   clustalo got killed on an 8GB laptop)
# - so we translate to protein, align proteins, then back-translate
#   with pal2nal. This is standard practice in phylogenetics anyway.
#
# IMPORTANT: the binary is `pal2nal.pl` not `pal2nal`, even though
# the conda package is called pal2nal. spent 20 minutes on that one.

set -euo pipefail

PROJ="$HOME/erbb_project"
cd "$PROJ/reference/orthologs"

# ---- translate each CDS to protein ----
# -f 1 = standard genetic code (not the vertebrate mitochondrial one)
mkdir -p proteins
for f in cds/*_cds.fasta; do
  base=$(basename "$f" _cds.fasta)
  seqkit translate -f 1 "$f" \
    | sed "s|^>.*|>${base}|" \
    > "proteins/${base}_protein.fasta"
done

# ---- rebuild CDS files with matching headers ----
# pal2nal matches sequences by name, so protein and CDS headers must
# line up exactly. easier to rewrite the CDS headers than to fight
# the multi-line FASTA headers NCBI gives us.
mkdir -p cds_norm
for f in cds/*_cds.fasta; do
  base=$(basename "$f" _cds.fasta)
  echo ">${base}" > "cds_norm/${base}.fasta"
  grep -v "^>" "$f" >> "cds_norm/${base}.fasta"
done

# ---- concatenate in consistent order ----
# alphabetical glob order is fine as long as BOTH the protein list and
# the CDS list use the same glob. They do.
cat proteins/*_protein.fasta > erbb_proteins_all.fasta
cat cds_norm/*.fasta        > erbb_cds_all.fasta

# ---- MAFFT ----
# L-INS-i is the slow/accurate mode. For 17 sequences it's a minute
# or so on 4 threads. Could use --auto if in a hurry.
mafft --localpair --maxiterate 1000 --thread 4 \
      erbb_proteins_all.fasta > erbb_proteins_aligned.fasta

# ---- back-translate to codons ----
pal2nal.pl erbb_proteins_aligned.fasta erbb_cds_all.fasta \
  -output fasta > erbb_codon_aligned.fasta

# verify - alignment should be a single number, divisible by 3
echo "alignment width check:"
awk '/^>/{if(seq!=""){print length(seq); seq=""}; next}{seq=seq$0} END{print length(seq)}' \
  erbb_codon_aligned.fasta | sort -u

# ---- format conversion for PAUP* and PHYLIP ----
# biopython's NEXUS writer wants an explicit molecule_type annotation
# that we don't have, so writing both formats by hand. it's just text.
python3 - <<'PYEOF'
from Bio import AlignIO

aln = AlignIO.read("erbb_codon_aligned.fasta", "fasta")
n_tax = len(aln)
n_char = aln.get_alignment_length()

with open("erbb_codon_aligned.nex", "w") as f:
    f.write("#NEXUS\nBEGIN DATA;\n")
    f.write(f"  DIMENSIONS NTAX={n_tax} NCHAR={n_char};\n")
    f.write("  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;\n")
    f.write("  MATRIX\n")
    for rec in aln:
        f.write(f"  {rec.id[:30]:<30} {str(rec.seq).upper()}\n")
    f.write("  ;\nEND;\n")

with open("erbb_codon_aligned.phy", "w") as f:
    f.write(f"  {n_tax}  {n_char}\n")
    for rec in aln:
        f.write(f"{rec.id[:30]:<30} {str(rec.seq).upper()}\n")
PYEOF

# ---- save outputs ----
mkdir -p "$PROJ/results/alignments"
cp erbb_codon_aligned.fasta erbb_codon_aligned.nex erbb_codon_aligned.phy \
   erbb_proteins_aligned.fasta erbb_proteins_all.fasta erbb_cds_all.fasta \
   "$PROJ/results/alignments/"
cp cds/*_cds.fasta "$PROJ/results/alignments/cds/"

echo "done. check $PROJ/results/alignments/"
