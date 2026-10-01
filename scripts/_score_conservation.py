#!/usr/bin/env python3
"""
Phase 5 helper - conservation scoring.

Maps variant cDNA positions from VEP to alignment columns and scores
each variant's conservation across the ErbB family tree.

The coordinate correction:
    Ensembl and RefSeq annotate CDS starts differently for some
    transcripts. This means VEP's cds_start does not always land on
    the correct alignment column. We search +/-20 columns around the
    expected position for a match to the VCF REF allele and pick the
    closest one.
"""

from Bio import AlignIO
import os

PROJ = os.path.expanduser("~/erbb_project")
ALIGNMENT = os.path.join(PROJ, "results/alignments/erbb_codon_aligned.fasta")
OUT_TABLE = "conservation_table_final.tsv"

# ----- load alignment -----
aln = AlignIO.read(ALIGNMENT, "fasta")
seqs = {rec.id: str(rec.seq).upper() for rec in aln}

def build_col_map(seq):
    """CDS position (1-based, ignoring gaps) -> alignment column index (0-based)."""
    m = {}
    cds_pos = 0
    for i, b in enumerate(seq):
        if b != '-':
            cds_pos += 1
            m[cds_pos] = i
    return m

col_maps = {sid: build_col_map(s) for sid, s in seqs.items()}

def find_column(human_seq, expected_col, ref_base, window=20):
    """Search +/-window for the column where human matches ref_base."""
    best = None
    for delta in range(0, window + 1):
        for col in (expected_col - delta, expected_col + delta):
            if 0 <= col < len(human_seq) and human_seq[col] == ref_base:
                if best is None or abs(col - expected_col) < abs(best - expected_col):
                    best = col
        if best is not None:
            return best
    return None

# ----- variant list (from VEP output) -----
# gene, chrom, pos, ref, alt, cds_start, codon, aa_change
VARIANTS = [
    ('EGFR',  'chr7',  55146655, 'C', 'T', 474,  158,  'N (syn)'),
    ('EGFR',  'chr7',  55181370, 'G', 'A', 2361, 787,  'Q (syn)'),
    ('EGFR',  'chr7',  55198724, 'T', 'C', 2709, 903,  'T (syn)'),
    ('EGFR',  'chr7',  55201256, 'A', 'G', 3015, 1005, 'E (syn)'),
    ('ERBB2', 'chr17', 39727784, 'C', 'G', 3508, 1170, 'P->A (missense)'),
    ('ERBB3', 'chr12', 56101207, 'G', 'A', 3348, 1116, 'R (syn)'),
]

GENE_SPECIES = {
    'EGFR':  ['EGFR_human', 'EGFR_mouse', 'EGFR_chicken', 'EGFR_zebrafish', 'EGFR_drosophila'],
    'ERBB2': ['ERBB2_human', 'ERBB2_mouse', 'ERBB2_chicken', 'ERBB2_zebrafish'],
    'ERBB3': ['ERBB3_human', 'ERBB3_mouse', 'ERBB3_chicken', 'ERBB3_zebrafish'],
    'ERBB4': ['ERBB4_human', 'ERBB4_mouse', 'ERBB4_chicken', 'ERBB4_zebrafish'],
}
GENE_HUMAN = {'EGFR':'EGFR_human','ERBB2':'ERBB2_human','ERBB3':'ERBB3_human','ERBB4':'ERBB4_human'}

def score_cons(human_base, bases):
    n_match = sum(1 for b in bases.values() if b == human_base) - 1
    n_total = len(bases) - 1
    if n_match == n_total: return 3
    if n_match >= 3:       return 2
    if n_match >= 1:       return 1
    return 0

# ----- process -----
rows = []
for gene, chrom, pos, ref, alt, cds_start, codon, aachange in VARIANTS:
    human_id = GENE_HUMAN[gene]
    human_seq = seqs[human_id]
    expected_col = col_maps[human_id].get(cds_start)
    if expected_col is None:
        print(f"WARNING: cds_start {cds_start} out of range for {human_id}")
        continue

    col = find_column(human_seq, expected_col, ref, window=20)
    if col is None:
        print(f"WARNING: no column matching {ref} near {cds_start} for {gene}")
        continue

    shift = col - expected_col
    bases = {sid: seqs[sid][col] for sid in GENE_SPECIES[gene] if sid in seqs}
    human_base = bases[human_id]
    score = score_cons(human_base, bases)

    rows.append({
        'gene': gene, 'chrom': chrom, 'pos': pos,
        'ref': ref, 'alt': alt, 'codon': codon,
        'aa_change': aachange, 'human_base': human_base,
        'score': score, 'shift': shift, 'species_bases': bases,
    })

# ----- write output -----
with open(OUT_TABLE, 'w') as f:
    f.write("gene\tchrom\tpos\tref\talt\tcodon\taa_change\thuman_base\tscore\tshift\tspecies_bases\n")
    for r in rows:
        sb = '|'.join(f"{sid}={b}" for sid, b in r['species_bases'].items())
        f.write(f"{r['gene']}\t{r['chrom']}\t{r['pos']}\t{r['ref']}\t{r['alt']}\t"
                f"{r['codon']}\t{r['aa_change']}\t{r['human_base']}\t{r['score']}\t"
                f"{r['shift']}\t{sb}\n")

# ----- summary -----
print(f"Scored {len(rows)} variants. Written to {OUT_TABLE}")
by_score = {}
for r in rows:
    by_score[r['score']] = by_score.get(r['score'], 0) + 1
print("Score distribution:")
for s in (3, 2, 1, 0):
    print(f"  Score {s}: {by_score.get(s, 0)}")
