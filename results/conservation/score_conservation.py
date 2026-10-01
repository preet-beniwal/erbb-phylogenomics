from Bio import AlignIO

aln = AlignIO.read('/home/preet/erbb_project/results/alignments/erbb_codon_aligned.fasta', 'fasta')

# Map: sequence ID -> {cds_position -> alignment_column}
# CDS position = number of non-gap bases seen up to that alignment column
seq_to_cols = {}
for rec in aln:
    seq = str(rec.seq).upper()
    mapping = {}
    cds_pos = 0
    for col_idx, base in enumerate(seq):
        if base != '-':
            cds_pos += 1
            mapping[cds_pos] = col_idx
    seq_to_cols[rec.id] = mapping

# Gene -> human alignment ID
gene_to_human = {
    'EGFR': 'EGFR_human',
    'ERBB2': 'ERBB2_human',
    'ERBB3': 'ERBB3_human',
    'ERBB4': 'ERBB4_human',
}

# Variants from VEP (gene, chrom, pos, ref, alt, cds_start, codon, aa_change)
variants = [
    ('EGFR',  'chr7',  55146655, 'C', 'T', 474,  158,  'N (syn)'),
    ('EGFR',  'chr7',  55181370, 'G', 'A', 2361, 787,  'Q (syn)'),
    ('EGFR',  'chr7',  55198724, 'T', 'C', 2709, 903,  'T (syn)'),
    ('EGFR',  'chr7',  55201256, 'A', 'G', 3015, 1005, 'E (syn)'),
    ('ERBB2', 'chr17', 39727784, 'C', 'G', 3508, 1170, 'P->A (missense)'),
    ('ERBB3', 'chr12', 56101207, 'G', 'A', 3348, 1116, 'R (syn)'),
]

# Which species belong to which gene
gene_species = {
    'EGFR':  ['EGFR_human', 'EGFR_mouse', 'EGFR_chicken', 'EGFR_zebrafish', 'EGFR_drosophila'],
    'ERBB2': ['ERBB2_human', 'ERBB2_mouse', 'ERBB2_chicken', 'ERBB2_zebrafish'],
    'ERBB3': ['ERBB3_human', 'ERBB3_mouse', 'ERBB3_chicken', 'ERBB3_zebrafish'],
    'ERBB4': ['ERBB4_human', 'ERBB4_mouse', 'ERBB4_chicken', 'ERBB4_zebrafish'],
}

# Score function
def conservation_score(human_base, species_bases):
    # species_bases is dict: species_id -> base
    # Verify human's base matches
    matches = {sid: b for sid, b in species_bases.items() if b == human_base}
    n_match = len(matches) - 1  # exclude human itself
    n_total = len(species_bases) - 1  # exclude human
    if n_match == n_total:
        return 3   # invariant
    if n_match >= 3:
        return 2   # highly conserved
    if n_match >= 1:
        return 1   # partially conserved
    return 0       # human-specific

# Run
print("=" * 90)
print(f"{'Gene':<6} {'Pos':<11} {'Ref>Alt':<8} {'Codon':<6} {'AaChange':<17} {'ConsScore':<10} {'Species bases'}")
print("=" * 90)

results = []
for gene, chrom, pos, ref, alt, cds_start, codon, aachange in variants:
    human_id = gene_to_human[gene]
    if human_id not in seq_to_cols:
        print(f"  ERROR: {human_id} not found in alignment")
        continue

    col_idx = seq_to_cols[human_id].get(cds_start)
    if col_idx is None:
        print(f"  ERROR: cds_start {cds_start} out of range for {human_id}")
        continue

    # Extract each species' base at this column
    species_bases = {}
    for sid in gene_species[gene]:
        if sid in seq_to_cols:
            # find the sequence
            for rec in aln:
                if rec.id == sid:
                    species_bases[sid] = str(rec.seq[col_idx]).upper()
                    break

    human_base = species_bases.get(human_id, '?')
    score = conservation_score(human_base, species_bases)

    # Format species bases as compact string
    base_str = ', '.join(f"{sid.split('_')[1][:3]}={b}" for sid, b in species_bases.items())
    print(f"{gene:<6} {pos:<11} {ref+'>'+alt:<8} {codon:<6} {aachange:<17} {score:<10} {base_str}")

    results.append((gene, chrom, pos, ref, alt, codon, aachange, human_base, species_bases, score))

# Write results
with open('conservation_table.tsv', 'w') as f:
    f.write("gene\tchrom\tpos\tref\talt\tcodon\taa_change\thuman_base\tscore\tspecies_bases\n")
    for row in results:
        gene, chrom, pos, ref, alt, codon, aachange, human_base, sb, score = row
        sb_str = '|'.join(f"{sid}={b}" for sid, b in sb.items())
        f.write(f"{gene}\t{chrom}\t{pos}\t{ref}\t{alt}\t{codon}\t{aachange}\t{human_base}\t{score}\t{sb_str}\n")

print()
print("=" * 90)
print(f"{len(results)} variants scored. Results saved to conservation_table.tsv")
print("=" * 90)

# Summary
by_score = {}
for row in results:
    s = row[9]
    by_score[s] = by_score.get(s, 0) + 1
print()
print("Score distribution:")
print(f"  Score 3 (invariant):   {by_score.get(3, 0)}")
print(f"  Score 2 (conserved):   {by_score.get(2, 0)}")
print(f"  Score 1 (partial):     {by_score.get(1, 0)}")
print(f"  Score 0 (human-specific): {by_score.get(0, 0)}")
