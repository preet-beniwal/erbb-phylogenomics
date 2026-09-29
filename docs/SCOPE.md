# Project Scope: ErbB/RTK Comparative Phylogenomics

## Target Gene Family
- EGFR (ErbB1) - NM_005228
- ERBB2 (HER2/Neu) - NM_004448
- ERBB3 (HER3) - NM_001982
- ERBB4 (HER4) - NM_005235

## Species Panel
- Homo sapiens (human, 0 MYA)
- Pan troglodytes (chimp, ~6 MYA)
- Mus musculus (mouse, ~90 MYA)
- Gallus gallus (chicken, ~310 MYA)
- Danio rerio (zebrafish, ~450 MYA)
- Drosophila melanogaster (outgroup, ~700 MYA)

## Human Data
- Exome: GIAB NA12878 (Garvan HiSeq 2x100 bp)
- Reference: GRCh38 (no-alt analysis set)
- Truth set: NIST v3.3 high-confidence VCF + BED

## Software
- DNASTAR Lasergene (SeqMan NGen, MegAlign Pro)
- PAUP*
- PHYLIP
- FastQC, Trimmomatic, bcftools
- FigTree, iTOL

## Analytical Question
Do human ErbB variants fall at evolutionarily conserved positions
across vertebrates, and does conservation correlate with ClinVar
pathogenicity classifications?

## Deliverables
1. Filtered annotated VCF for ErbB genes
2. NEXUS and PHYLIP alignments per ErbB paralog
3. PAUP* and PHYLIP consensus trees with bootstrap
4. Conservation score table per human SNP
5. Annotated tree figure
6. GitHub repository with README and pipeline diagram
