#!/usr/bin/env bash
# Phase 1 - grab everything we need
# Run from anywhere but keep ~/erbb_project as the project root
# NOTE: switched everything from ftp to https because the ftp-trace PASV
#       responses were getting blocked on my home network. Lost an
#       evening to that. If you are on a corporate network you will hit
#       the same thing - stick with https URLs.

set -euo pipefail

PROJ="$HOME/erbb_project"
mkdir -p "$PROJ/data" "$PROJ/reference/orthologs/cds"

# ---- 1. NA12878 exome FASTQ (Garvan, 2x100bp) ----
cd "$PROJ/data"
for R in R1 R2; do
  wget -c --tries=10 --timeout=60 --waitretry=10 \
    "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_${R}_001.fastq.gz" \
    -O "NIST7035_${R}.fastq.gz"
done

# ---- 2. GRCh38 no-alt analysis set ----
# using the UCSC-style "no_alt" version, not the full primary assembly
# (the alt contigs add junk hits for a project this small)
cd "$PROJ/reference"
wget -c --tries=10 --timeout=60 \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" \
  -O GRCh38_no_alt.fna.gz

# ---- 3. GIAB truth set for NA12878 (NIST v3.3, GRCh37 coords) ----
# we later realised this is GRCh37 and got a GRCh38 version separately.
# keeping both in the repo just in case.
wget -c \
  "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.gz" \
  -O GIAB_truth_v3.3.vcf.gz
wget -c \
  "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.bed" \
  -O GIAB_truth_regions.bed

# ---- 4. ErbB family CDS orthologs ----
# Originally I hardcoded accessions but half of them were wrong genes
# (got a zinc-finger protein for Drosophila Egfr, a claudin for
# zebrafish Erbb2, etc). Now querying NCBI by gene symbol + taxid
# and letting the database pick the best RefSeq mRNA.
#
# Also: chimp EGFR is not curated at full length on NCBI - every
# transcript has N gaps. Dropped chimp from the panel. Not worth
# fighting; 6 Mya divergence adds nothing to a conservation analysis
# that already spans 450 Myr.

cd "$PROJ/reference/orthologs/cds"

declare -A TAX=(
  [human]=9606
  [mouse]=10090
  [chicken]=9031
  [zebrafish]=7955
)

for gene in EGFR ERBB2 ERBB3 ERBB4; do
  for sp in human mouse chicken zebrafish; do
    tx=${TAX[$sp]}
    q="${gene}%5BGene%5D+AND+txid${tx}%5BOrganism%5D+AND+refseq%5BFilter%5D+AND+biomol_mrna%5BPROP%5D"
    id=$(wget -qO- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=${q}&retmax=1&sort=relevance" \
         | grep -oP '(?<=<Id>)[0-9]+' | head -1)
    wget -qO "${gene}_${sp}_cds.fasta" \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${id}&rettype=fasta_cds_na&retmode=text"
    sleep 0.4
  done
done

# Drosophila Egfr as outgroup. txid 7227.
q="Egfr%5BGene%5D+AND+txid7227%5BOrganism%5D+AND+refseq%5BFilter%5D+AND+biomol_mrna%5BPROP%5D"
id=$(wget -qO- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=${q}&retmax=1&sort=relevance" \
     | grep -oP '(?<=<Id>)[0-9]+' | head -1)
wget -qO "EGFR_drosophila_cds.fasta" \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${id}&rettype=fasta_cds_na&retmode=text"

# sanity check - should be 17
echo "CDS files fetched: $(ls *_cds.fasta | wc -l)"
