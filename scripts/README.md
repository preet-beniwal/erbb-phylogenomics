cat > ~/erbb_project/scripts/README.md << 'EOF'
# Pipeline Scripts

Executable scripts for the full ErbB phylogenomics pipeline.

## Order of execution

| # | Script | Phase | Runtime |
|---|--------|-------|---------|
| 00 | (setup) | Environment | 15 min |
| 01 | `01_download_data.sh` | Data acquisition | 1-4 hrs (network-bound) |
| 02 | `02_qc_and_trim.sh` | QC + trimming | 15 min |
| 03 | `03_msa_pipeline.sh` | MSA + format conversion | 2 min |
| 04 | `04_variant_calling.sh` | BWA + bcftools | 2.5 hrs |
| 05 | `05_phylogenetics.sh` | PHYLIP + IQ-TREE | 30 min |
| 06 | `06_conservation_analysis.sh` | VEP + conservation scoring | 1 min |

Helper files:
- `_score_conservation.py` - coordinate correction and scoring logic used by script 06

## Prerequisites

Environment set up as described in `docs/SCOPE.md`. Briefly:

```bash
conda create -n erbb python=3.11 -y
conda activate erbb
conda install -c bioconda -c conda-forge \
  fastqc trimmomatic bwa samtools bcftools tabix \
  muscle clustalo mafft pal2nal seqkit emboss phylip iqtree \
  -y

