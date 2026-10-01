#!/usr/bin/env bash
# Phase 2 - variant calling pipeline
#
# BWA alignment + bcftools variant calling. About 2.5 hours total
# on the 16GB laptop, mostly unattended.
#
# Prerequisites:
#   - conda env "erbb" active with bwa, samtools, bcftools
#   - Trimmed FASTQ in data/
#   - GRCh38 reference in reference/
#
# Notes from the actual run:
#   - bwa index on GRCh38 took 88 min (BWT construction is the long part)
#   - alignment of 18.5M read pairs: 40 min at 4 threads
#   - bcftools mpileup: 21 min, produced 4.6 GB intermediate BCF
#   - final ErbB variant count: 15
#
# Learned the hard way: do NOT chain "bcftools mpileup | bcftools call"
# under nohup. The pipe gets detached and call tries to open "-" as a
# file. Split into two steps with an intermediate BCF - costs 4.6GB
# disk but works every time.

set -euo pipefail

PROJ="$HOME/erbb_project"
cd "$PROJ"
mkdir -p results/bwa results/vcf

# --- BWA index (skip if already built) ---
if [ ! -f reference/GRCh38_no_alt.fna.bwt ]; then
  echo "[$(date)] Building BWA index (this takes ~90 min)..."
  bwa index reference/GRCh38_no_alt.fna
fi

# --- Alignment + sort in one pass ---
# streaming into samtools sort avoids a 40GB intermediate SAM.
# 1GB memory buffer is plenty for a 2.6GB sorted BAM.
echo "[$(date)] Running bwa mem + samtools sort..."
bwa mem -t 4 \
  reference/GRCh38_no_alt.fna \
  data/R1_trimmed.fastq.gz \
  data/R2_trimmed.fastq.gz \
  2> results/bwa/bwa_mem.log \
| samtools sort -@ 2 -m 1G \
    -o results/bwa/NA12878_erbb.sorted.bam \
    -T results/bwa/sort_tmp -

samtools index results/bwa/NA12878_erbb.sorted.bam
samtools flagstat results/bwa/NA12878_erbb.sorted.bam > results/bwa/flagstat.txt

# --- Pileup (separate step - see header note) ---
echo "[$(date)] Running bcftools mpileup..."
bcftools mpileup \
  -f reference/GRCh38_no_alt.fna \
  -q 20 -Q 20 -d 250 \
  -a FORMAT/AD,FORMAT/DP,INFO/AD \
  -Ob -o results/vcf/NA12878_erbb.pileup.bcf \
  results/bwa/NA12878_erbb.sorted.bam

# --- Variant calling ---
# -mv = multiallelic caller, variants only
# the "assuming diploid" note is expected - human autosomes are diploid
echo "[$(date)] Running bcftools call..."
bcftools call -mv -Oz \
  -o results/vcf/NA12878_erbb.vcf.gz \
  results/vcf/NA12878_erbb.pileup.bcf

bcftools index -t results/vcf/NA12878_erbb.vcf.gz

# Free the 4.6GB intermediate - regenerable from the BAM
rm -f results/vcf/NA12878_erbb.pileup.bcf

# --- Filter ---
# unfiltered call has 610K variants because we never restricted to the
# exome capture BED. Quality filter drops to ~74K, then ErbB subsetting
# gives the biologically meaningful 15.
echo "[$(date)] Filtering..."
bcftools filter \
  -e 'QUAL<30 || FORMAT/DP<10 || MQ<30' \
  -Oz -o results/vcf/NA12878_filtered.vcf.gz \
  results/vcf/NA12878_erbb.vcf.gz
bcftools index -t results/vcf/NA12878_filtered.vcf.gz

# --- ErbB subset ---
# gene body +/- 5kb, GRCh38 coordinates from Ensembl
echo "[$(date)] Subsetting to ErbB genes..."
bcftools view \
  -r chr7:55014000-55217000,chr17:39682000-39734000,chr12:56071000-56105000,chr2:211370000-212544000 \
  -Oz -o results/vcf/NA12878_erbb_family.vcf.gz \
  results/vcf/NA12878_filtered.vcf.gz
bcftools index -t results/vcf/NA12878_erbb_family.vcf.gz

# Save the variant table for downstream analysis
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\n' \
  results/vcf/NA12878_erbb_family.vcf.gz > results/vcf/erbb_variants.tsv

echo "[$(date)] Done. ErbB variants: $(bcftools view -H results/vcf/NA12878_erbb_family.vcf.gz | wc -l)"
