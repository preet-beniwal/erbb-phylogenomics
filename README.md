# Comparative Phylogenomics of the ErbB/RTK Gene Family

**A variant interpretation pipeline combining exome sequencing, cross-species multiple sequence alignment, and Bayesian/parsimony phylogenetic reconstruction to assess evolutionary conservation at human mutation sites.**

---

## Project Overview

This project investigates whether human single nucleotide variants (SNVs) in the ErbB receptor tyrosine kinase family (EGFR, ERBB2, ERBB3, ERBB4) occur at evolutionarily conserved positions across vertebrates, and whether conservation status correlates with known clinical pathogenicity annotations.

## Target Gene Family

| Gene | Protein | Clinical Relevance |
|------|---------|-------------------|
| EGFR | ErbB1 | Lung adenocarcinoma, glioblastoma |
| ERBB2 | HER2/Neu | Breast cancer, gastric cancer |
| ERBB3 | HER3 | Colorectal cancer, therapeutic resistance |
| ERBB4 | HER4 | Neural development, tumor suppressor |

## Species Panel

*Homo sapiens*, *Pan troglodytes*, *Mus musculus*, *Gallus gallus*, *Danio rerio*, with *Drosophila melanogaster* as outgroup.

## Analytical Pipeline

1. **Variant calling** — GIAB NA12878 exome aligned to GRCh38 in DNASTAR SeqMan NGen
2. **Ortholog retrieval** — NCBI RefSeq CDS for each ErbB paralog across six species
3. **Multiple sequence alignment** — MUSCLE/Clustal Omega via DNASTAR MegAlign Pro
4. **Phylogenetic reconstruction** — Parallel analysis in PAUP* (Maximum Parsimony + 1000 bootstrap) and PHYLIP (Neighbor-Joining with distance matrix bootstrapping)
5. **Conservation analysis** — Human SNPs mapped onto consensus tree; conservation scored 0–3; correlation with ClinVar/gnoMAD assessed

## Status

- [x] Phase 0: Environment setup
- [ ] Phase 1: Data acquisition and quality control
- [ ] Phase 2: Variant calling and filtering
- [ ] Phase 3: Multiple sequence alignment
- [ ] Phase 4: Phylogenetic reconstruction
- [ ] Phase 5: Conservation analysis
- [ ] Phase 6: Manuscript-style writeup and figures

## Author

**Preet Beniwal** — Project portfolio for PhD applications in Bioinformatics / Computational Biology

## License

See `LICENSE`.
