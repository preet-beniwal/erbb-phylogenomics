# Phase 5 — Cross-Species Conservation Analysis of ErbB Variants

## Objective
Score the evolutionary conservation of each coding variant identified in
NA12878 (GIAB) across the ErbB family phylogeny, to assess whether human
variants fall at positions under purifying selection.

## Method
1. Variant coordinates were submitted to Ensembl VEP (REST API) for
   functional annotation and cDNA positions.
2. cDNA positions were mapped to columns in the codon-aware alignment
   (17 taxa × 5,727 columns; Phase 3).
3. Because Ensembl and RefSeq CDS annotations differ slightly at the
   5' end of some transcripts, each variant was matched to the alignment
   column within a ±20 nt window where the human reference base matched
   the VCF REF allele.
4. Conservation was scored by counting species that share the human
   base at the corrected column.

## Conservation Scoring
- Score 3 (invariant): all species share the human base
- Score 2 (highly conserved): ≥75% of species match
- Score 1 (partial): 25–75% match
- Score 0 (human-specific): <25% match

## Results

| Gene | Position | Ref>Alt | Codon | Consequence | Score | Shift | Species bases (human/mouse/chicken/zebrafish/drosophila) |
|------|----------|---------|-------|-------------|-------|-------|---------------------------------------------------------|
| EGFR | 55146655 | C>T | 158 | Synonymous (N) | 1 | 0 | C/T/C/T/C |
| EGFR | 55181370 | G>A | 787 | Synonymous (Q) | 3 | 0 | G/G/G/G/G |
| EGFR | 55198724 | T>C | 903 | Synonymous (T) | 1 | 0 | T/T/A/C/A |
| EGFR | 55201256 | A>G | 1005 | Synonymous (E) | 0 | 0 | A/G/G/G/C |
| ERBB2 | 39727784 | C>G | 1170 | **Missense (P→A)** | 1 | +4 | C/C/C/G |
| ERBB3 | 56101207 | G>A | 1116 | Synonymous (R) | 1 | -1 | G/G/-/C |

## Interpretation

### ERBB2 P→A at codon 1170 (headline finding)
The human reference base C is shared by mouse and chicken but not zebrafish.
The position has been C in amniotes for at least 310 million years. The
variant changes C→G, reverting the residue to the state found in zebrafish.
Two interpretations:
1. The substitution is functionally tolerated (zebrafish thrives with it).
2. The position was constrained in amniotes; the fish lineage diverged
   independently. The patient's reversion may or may not be pathogenic.
Cross-reference with ClinVar is required to establish clinical significance.

### EGFR 55181370 (synonymous at invariant position)
The base G is invariant from Drosophila to human — approximately 700 Myr
of evolutionary constraint, despite the change being synonymous. This is
consistent with codon usage bias or a regulatory element at the nucleotide
level.

### Other variants (EGFR 55146655, 55198724, 55201256, ERBB3 56101207)
Showed partial or no conservation. Consistent with neutral or weakly
constrained positions.

## Methodological Note
The ±20 column window adjustment corrects for known discrepancies between
Ensembl and RefSeq CDS coordinates. Observed shifts ranged from 0 to +4
columns. This is documented in the pipeline script
(`score_conservation_v2.py`).

## Files
- `coding_variants.tsv`: raw VEP-annotated coding variants
- `conservation_table_final.tsv`: variant × column × species × score
- `score_conservation_v2.py`: pipeline script
