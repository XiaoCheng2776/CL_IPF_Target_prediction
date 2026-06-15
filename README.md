# CL_IPF_Target_prediction

scRNA-seq target-discovery pipeline for **idiopathic pulmonary fibrosis (IPF)**.  
Goal: rank cell-surface receptors on AT2 alveolar epithelium whose **agonism** promotes successful AT2 → Krt8⁺ADI → AT1 regeneration and resolves fibrosis.

> Full results and interpretation: **[RESULTS_SUMMARY.md](RESULTS_SUMMARY.md)**

---

## Background

In IPF, alveolar type-2 (AT2) progenitor cells fail to complete regeneration and instead arrest in a Krt8⁺ basaloid/ADI state that drives progressive fibrosis. The therapeutic goal is to identify receptors on AT2 cells that, when stimulated, allow cells to exit the arrested state and reach the AT1 destination.

**Positive control:** FGF10 → FGFR2b axis. Fgf10-deficient mice develop fibrosis more readily; FGF10 agonism rescues AT2 regeneration in vivo.

**Biology is inverted vs a proliferation/EMT screen:** EMT and Krt8⁺ expansion are disease features here, not therapeutic goals. Target phenotype = trajectory completion (AT2 → AT1), not activation magnitude.

---

## Pipeline

```
01_download.sh          — idempotent data download (mouse + human datasets)
02_qc.R                 — per-dataset QC, doublet removal, ambient RNA
03_phenotype.R          — epithelial subset, Slingshot pseudotime, fate labels
04_nichenet.R           — NicheNet v2 ligand-activity scoring (mouse)
05_cellchat.R           — CellChat v2 receptor communication probabilities
06_geneformer.py        — Geneformer V1-10M fine-tune + in silico perturbation
07_crossspecies.R       — composite scoring, cross-species conservation
08_human_nichenet.R     — NicheNet v2 human (Habermann GSE135893)
```

Run scripts in order; inspect outputs before proceeding. See [CLAUDE.md](CLAUDE.md) for per-step instructions and guardrails.

---

## Data

| Dataset | Species | Accession | Role |
|---|---|---|---|
| Strunz bleomycin time course | Mouse | GSE141259 | Primary — sender characterization + trajectory |
| Kobayashi PATS | Mouse | GSE141634 | AT2 receiver (sorted, epithelial-enriched) |
| Choi DATPs | Mouse | GSE145031 | AT2 receiver (sorted) |
| Habermann IPF | Human | GSE135893 | Cross-species validation |
| Adams IPF/Control | Human | GSE136831 | Optional validation (large; pre-filter required) |

Raw data → `data/` (git-ignored). Results → `results/`.

---

## Key Results

### Fate labels (03_phenotype.R)
- **completing** (n = 2,685): 36% AT1 cells + AT2/AT2-activated with high Ager/Pdpn/Hopx. RegenSuccess score = +0.976.
- **arrested** (n = 2,332): 85% Krt8⁺ ADI. Transitional score = +0.652.
- geneset_oi = completing vs arrested (bleomycin cells only; d14_PBS excluded).

### Top agonist candidates

| Receptor | Composite | Conservation | Key ligands |
|---|---|---|---|
| Sdc1 | 0.850 | mouse + human | ANGPTL4, LPL (FGFR2b co-receptor) |
| Erbb2 | 0.821 | mouse + human | AREG, HBEGF |
| Egfr | 0.796 | mouse + human | AREG, HBEGF, CDCP1 |
| Fzd2 | 0.722 | mouse-only | WNT (top Geneformer OE) |
| Sdc4 | 0.710 | mouse + human | ANGPTL4 (FGFR2b co-receptor) |
| **Fgfr2** | 0.603 | mouse + human | FGF7, FGF10 — **positive control** |
| Cd9 | 0.551 | mouse + human | HBEGF |
| Nrp1 | 0.437 | mouse-only | FGF7, VEGFA (FGFR2b co-receptor) |

### Block candidates (TGFβ-driven; antagonize to release arrest)

| Receptor | Rationale |
|---|---|
| App | Top ligand = TGFB1 (both species); TGFβ→APP promotes AT2 senescence |
| Tgfbr2 | Direct TGFβ receptor → SMAD2/3 → Krt8/pEMT/fibrosis |
| Tgfbr1 | Direct TGFβ type-I receptor (ALK5) |

> Note: TGFB1 ranks #1 in NicheNet because it is the most active sender signal in fibrosis — but it drives the **arrested** state. TGFβ receptors are block candidates, not agonist targets.

See [RESULTS_SUMMARY.md](RESULTS_SUMMARY.md) for full ranked tables, scoring details, and caveats.

---

## Environment

- **R 4.3+**: Seurat v5, nichenetr, CellChat, slingshot, tidyverse
- **Python 3.10+** (conda env `gf`): scanpy, geneformer, transformers, torch
- Geneformer fine-tuning requires GPU ≥ 16 GB VRAM

---

## Caveats

- **FGFR2 isoform:** Standard Fgfr2 quantification combines IIIb (epithelial, pro-regenerative) and IIIc (mesenchymal). Therapeutic strategies must use FGFR2b-selective agents.
- **NicheNet v2 Fgf10 prior:** Fgf10's top predicted targets in the mouse lt_matrix are prolactin-family mammary genes (not AT2 biology). Fgf10 ranks 27; Fgf7 (same receptor) ranks 3. Gate adjusted accordingly.
- **Habermann is IPF-only:** No healthy AT2 reference; FGF10 is nearly absent from established-IPF fibroblasts (disease biology, not artifact).
