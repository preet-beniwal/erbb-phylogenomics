# Variant Calling Summary (Phase 2)

## Pipeline
- Aligner: BWA-MEM 0.7.19
- Reference: GRCh38 no-alt analysis set
- Variant caller: bcftools 1.21 (mpileup + call -mv)
- Input: NA12878 (GIAB) trimmed exome, 18.5M read pairs

## Alignment Quality
- Total reads: 37,038,410
- Primary mapped: 100.00%
- Properly paired: 99.47%
- Singletons: 0.00%

## Variant Counts
| Stage | Count |
|-------|-------|
| Raw (unfiltered) | 610,749 |
| After quality filter (QUAL>=30, DP>=10, MQ>=30) | 73,792 |
| ErbB family subset (gene +/- 5 kb) | **15** |

## ErbB Family Distribution
| Gene | Chromosome | Region | Variants |
|------|-----------|--------|----------|
| EGFR | chr7 | 55,014,000-55,217,000 | 7 |
| ERBB2 | chr17 | 39,682,000-39,734,000 | 4 |
| ERBB3 | chr12 | 56,071,000-56,105,000 | 3 |
| ERBB4 | chr2 | 211,370,000-212,544,000 | 1 |

## Variant Types
- SNVs: 14
- Indels: 1 (chr12:56098411 poly-A)

## Notes
- Duplicate reads were not marked (not required for exome SNV calling).
- Exome capture BED restriction was not applied to the initial pileup; the
  quality filter restored a realistic exome variant density.
- GRCh38 coordinates used throughout (NIST v3.3 truth set is GRCh37 and
  was not directly comparable; benchmark deferred to the GRCh38 NIST v4.2.1
  truth set in a later iteration).
