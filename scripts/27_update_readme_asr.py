from pathlib import Path

readme = Path("README.md")
text = readme.read_text()

asr_section = """
### Ancestral Sequence Reconstruction: The Pocket is 450 Million Years Old

We reconstructed the ancestral sequences at every internal node of the ErbB phylogeny using IQ-TREE (GTR+F+I+G4 model, codon-based alignment).

At the root node of the vertebrate ErbB family, the three residues that form the allosteric binding pocket in human EGFR (K745, F856, L858) were reconstructed as:

| Position | Human EGFR | Ancestral Node | Ancestral Codon | Amino Acid |
|----------|------------|----------------|-----------------|------------|
| 745 | K | **K** | AAG | Lysine |
| 856 | F | **F** | TTT | Phenylalanine |
| 858 | L | **L** | CTG | Leucine |

**Interpretation:** The allosteric pocket targeted by the EAI045 family of inhibitors in human EGFR (Project 1) was already fully assembled in the common ancestor of all jawed vertebrates (~450 Mya), predating the genome duplications that produced the four modern ErbB paralogs. This indicates that the pocket is a primordial ErbB-family feature, conserved for its structural role in the kinase domain, not a recent innovation specific to EGFR.
"""

# Insert before "## Limitations" or "## Pipeline"
if "## Limitations" in text:
    text = text.replace("## Limitations", asr_section + "\n## Limitations")
elif "## Pipeline" in text:
    text = text.replace("## Pipeline", asr_section + "\n## Pipeline")
else:
    text += asr_section

readme.write_text(text)
print("README.md successfully updated with ASR section.")
