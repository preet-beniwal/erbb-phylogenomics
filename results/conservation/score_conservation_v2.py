from Bio import AlignIO

aln = AlignIO.read('/home/preet/erbb_project/results/alignments/erbb_codon_aligned.fasta', 'fasta')
seqs = {rec.id: str(rec.seq).upper() for rec in aln}

def find_column(human_seq, expected_col, ref_base, window=20):
    """Search ±window around expected_col for a column where human matches ref_base.
       Prefer the closest match; return None if no match found."""
    best = None
    for delta in range(0, window+1):
        for col in (expected_col - delta, expected_col + delta):
            if 0 <= col < len(human_seq) and human_seq[col] == ref_base:
                if best is None or abs(col - expected_col) < abs(best - expected_col):
                    best = col
        if best is not None:
            return best
    return None

# Build CDS position -> column maps per species
def build_col_map(seq):
    m = {}
    cds_pos = 0
    for i, b in enumerate(seq):
        if b != '-':
            cds_pos += 1
            m[cds_pos] = i
    return m

col_maps = {sid: build_col_map(s) for sid, s in seqs.items()}

variants = [
    ('EGFR',  'chr7',  55146655, 'C', 'T', 474,  158,  'N (syn)'),
    ('EGFR',  'chr7',  55181370, 'G', 'A', 2361, 787,  'Q (syn)'),
    ('EGFR',  'chr7',  55198724, 'T', 'C', 2709, 903,  'T (syn)'),
    ('EGFR',  'chr7',  55201256, 'A', 'G', 3015, 1005, 'E (syn)'),
    ('ERBB2', 'chr17', 39727784, 'C', 'G', 3508, 1170, 'P->A (missense)'),
    ('ERBB3', 'chr12', 56101207, 'G', 'A', 3348, 1116, 'R (syn)'),
]

gene_species = {
    'EGFR':  ['EGFR_human', 'EGFR_mouse', 'EGFR_chicken', 'EGFR_zebrafish', 'EGFR_drosophila'],
    'ERBB2': ['ERBB2_human', 'ERBB2_mouse', 'ERBB2_chicken', 'ERBB2_zebrafish'],
    'ERBB3': ['ERBB3_human', 'ERBB3_mouse', 'ERBB3_chicken', 'ERBB3_zebrafish'],
}
gene_to_human = {'EGFR':'EGFR_human','ERBB2':'ERBB2_human','ERBB3':'ERBB3_human','ERBB4':'ERBB4_human'}

def score_cons(human_base, bases):
    n_match = sum(1 for b in bases.values() if b == human_base) - 1
    n_total = len(bases) - 1
    if n_match == n_total: return 3
    if n_match >= 3:       return 2
    if n_match >= 1:       return 1
    return 0

print("=" * 105)
print(f"{'Gene':<6} {'Pos':<11} {'Ref>Alt':<8} {'Codon':<6} {'AaChange':<17} {'Score':<6} {'Shift':<7} Species bases")
print("=" * 105)

results = []
for gene, chrom, pos, ref, alt, cds_start, codon, aachange in variants:
    human_id = gene_to_human[gene]
    human_seq = seqs[human_id]

    # Try coordinate first
    expected_col = col_maps[human_id].get(cds_start)
    if expected_col is None:
        expected_col = len(human_seq) // 2  # fallback

    # Search for ref base within ±20 columns
    col = find_column(human_seq, expected_col, ref, window=20)

    if col is None:
        print(f"{gene:<6} {pos:<11} {ref+'>'+alt:<8} {codon:<6} {aachange:<17} UNKNOWN (no column matching {ref} near {cds_start})")
        continue

    shift = col - expected_col

    bases = {}
    for sid in gene_species[gene]:
        if sid in seqs:
            bases[sid] = seqs[sid][col]

    human_base = bases[human_id]
    score = score_cons(human_base, bases)

    base_str = ', '.join(f"{sid.split('_')[1][:3]}={b}" for sid, b in bases.items())
    shift_str = f"+{shift}" if shift > 0 else str(shift)
    print(f"{gene:<6} {pos:<11} {ref+'>'+alt:<8} {codon:<6} {aachange:<17} {score:<6} {shift_str:<7} {base_str}")

    results.append((gene, chrom, pos, ref, alt, codon, aachange, human_base, score, bases, shift))

with open('conservation_table_final.tsv', 'w') as f:
    f.write("gene\tchrom\tpos\tref\talt\tcodon\taa_change\thuman_base\tscore\tshift\tspecies_bases\n")
    for r in results:
        gene, chrom, pos, ref, alt, codon, aachange, hb, score, bases, shift = r
        sb = '|'.join(f"{sid}={b}" for sid, b in bases.items())
        f.write(f"{gene}\t{chrom}\t{pos}\t{ref}\t{alt}\t{codon}\t{aachange}\t{hb}\t{score}\t{shift}\t{sb}\n")

print()
print("=" * 105)
print(f"{len(results)} variants scored with coordinate correction. Saved to conservation_table_final.tsv")
print("=" * 105)

by_score = {}
for r in results:
    by_score[r[8]] = by_score.get(r[8], 0) + 1
print()
print("Score distribution:")
for s in (3, 2, 1, 0):
    print(f"  Score {s}: {by_score.get(s, 0)}")
