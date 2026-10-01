#!/usr/bin/env bash
# Phase 5 - conservation analysis of ErbB coding variants
#
# Queries Ensembl VEP for functional annotations, maps the coding
# variants to alignment columns, scores conservation across the
# 17-taxon alignment.
#
# Gotcha: Ensembl and RefSeq assign different CDS start positions for
# some transcripts, so VEP's cds_start does not always map 1:1 to
# columns in the RefSeq-based alignment. Solved by searching +/-20
# columns around the expected position for a match to the VCF REF
# allele. See _score_conservation.py for the implementation.

set -euo pipefail

PROJ="$HOME/erbb_project"
SCRIPTS="$PROJ/scripts"
mkdir -p "$PROJ/results/conservation"
cd "$PROJ/results/conservation"

# --- Extract SNVs (skip poly-A indel, VEP REST doesn't handle it well) ---
echo "[$(date)] Extracting SNVs..."
bcftools view -v snps "$PROJ/results/vcf/NA12878_erbb_family.vcf.gz" | \
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > snvs.tsv
echo "SNVs to annotate: $(wc -l < snvs.tsv)"

# --- Query Ensembl VEP REST API ---
# free, no install, no key. Just POST JSON and get JSON back.
echo "[$(date)] Querying Ensembl VEP..."
python3 - << 'PYEOF'
import json
variants=[]
with open('snvs.tsv') as f:
    for line in f:
        c,p,r,a=line.strip().split('\t')
        variants.append(f'{c} {p} . {r} {a} . . .')
json.dump({'variants':variants}, open('vep_payload.json','w'))
print(f"Prepared {len(variants)} variants")
PYEOF

curl -s 'https://rest.ensembl.org/vep/human/region' \
  -H 'Content-Type: application/json' \
  -H 'Accept: application/json' \
  -X POST \
  -d @vep_payload.json > vep_response.json

echo "VEP response: $(wc -c < vep_response.json) bytes"

# --- Score conservation ---
echo "[$(date)] Scoring conservation..."
python3 "$SCRIPTS/_score_conservation.py"

echo ""
echo "Conservation table:"
head -15 conservation_table_final.tsv

echo ""
echo "[$(date)] Phase 5 complete."
