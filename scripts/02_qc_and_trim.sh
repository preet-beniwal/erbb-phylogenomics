#!/usr/bin/env bash
# Phase 1 - FastQC + Trimmomatic
# needs the erbb conda env active
# 18.5M read pairs survived trimming, ~91.7% retention which is fine

set -euo pipefail

PROJ="$HOME/erbb_project"
cd "$PROJ/data"

# ---- FastQC ----
# reports saved as html + zip, open the html in a browser to look
mkdir -p "$PROJ/results/qc"
fastqc NIST7035_R1.fastq.gz NIST7035_R2.fastq.gz \
  -o "$PROJ/results/qc/" -t 4

# ---- Trimmomatic ----
# this is the standard TruSeq3-PE adapter file that ships with trimmomatic
# parameters:
#   LEADING:3 / TRAILING:3  - chop bases below Q3 at the ends
#   SLIDINGWINDOW:4:20      - scan 4bp windows, cut when avg Q<20
#   MINLEN:36               - drop reads shorter than 36bp
# adapter file was at $CONDA_PREFIX/share/trimmomatic-0.41-0/adapters/
trimmomatic PE -threads 4 \
  NIST7035_R1.fastq.gz NIST7035_R2.fastq.gz \
  R1_trimmed.fastq.gz R1_unpaired.fastq.gz \
  R2_trimmed.fastq.gz R2_unpaired.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# should see "Both Surviving: ~91%" in the output
