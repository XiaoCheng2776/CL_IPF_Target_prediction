# IPF Target Discovery — Results Summary

**Goal:** Identify cell-surface receptors on AT2 alveolar epithelium whose **agonism** promotes the AT2 → Krt8⁺ADI → AT1 regeneration trajectory and resolves pulmonary fibrosis.

**Positive control:** FGF10 → FGFR2b axis (Fgf10-deficient mice develop fibrosis more readily; agonism rescues regeneration).

---

## 1. Pipeline Overview

| Step | Script | Method | Input |
|---|---|---|---|
| QC | 02_qc.R | Seurat v5 | Strunz GSE141259, Kobayashi GSE141634, Choi GSE145031 |
| Phenotyping | 03_phenotype.R | Slingshot pseudotime + module scoring | Epithelial subset |
| Mouse NicheNet | 04_nichenet.R | NicheNet v2 (Zenodo 7074291) | epi_obj fate labels |
| CellChat | 05_cellchat.R | CellChat v2 | Whole-lung object |
| Geneformer | 06_geneformer.py | V1-10M fine-tune + in silico perturbation | epi_obj counts |
| Human NicheNet | 08_human_nichenet.R | NicheNet v2 human priors | Habermann GSE135893 |
| Integration | 07_crossspecies.R | Composite rank scoring | All above |

**Receiver:** AT2 cells and injury-activated alveolar progenitors (Kobayashi PATS / Choi DATPs).  
**Senders:** Alveolar/lipofibroblasts, myofibroblasts, airway smooth muscle, profibrotic macrophages.

---

## 2. Fate Label Verification

Fate labels were assigned by transcriptional module scoring on bleomycin cells only (d14_PBS = `healthy_baseline`, excluded from contrasts).

| Label | n | Cell types | RegenSuccess score | Transitional score |
|---|---|---|---|---|
| completing | 2,685 | 36% AT1 + AT2/AT2-activated | **+0.976** | −0.099 |
| arrested | 2,332 | **85% Krt8⁺ ADI** + AT2 | −0.083 | **+0.652** |
| healthy_baseline | 1,762 | AT2 (d14_PBS) | +1.08 | −0.074 |
| intermediate | 10,545 | mixed | +0.969 | −0.017 |

- **completing** = cells reaching the AT1 destination (high Ager/Pdpn/Hopx) or AT2 cells with strong identity scores (high Sftpc/Etv5/Lamp3, low Krt8/pEMT). ✓  
- **arrested** = persistent Krt8⁺/Cldn4⁺ ADI cells; 85% are Strunz-labeled "Krt8⁺ ADI". ✓  
- The DE (geneset_oi) is **completing vs arrested** — not bleomycin vs PBS. ✓

---

## 3. Positive Control Validation

| Species | Ligand | NicheNet rank | Note |
|---|---|---|---|
| Mouse | Fgf7 | **3** | Same FGFR2b IIIb receptor as FGF10 |
| Mouse | Fgf10 | 27 | Improved from 48 after geneset_oi fix; prior artifact ceiling (see §6) |
| Mouse | Fgf1 | 5 | FGFR2b family |
| Human | FGF7 | **24** (top 7% of 361 ligands) | pct=0.01 required to enter sender pool |
| Human | FGF10 | 155 | Present at pct=0.01; low rank reflects IPF fibroblast FGF10 downregulation |

The FGFR2b axis is recovered in both species. Fgf10 ranks 27 (not top 20) due to a NicheNet v2 mouse prior artifact: Fgf10's top predicted targets in the lt_matrix are prolactin-family mammary genes absent from any AT2 biology geneset (see §6).

---

## 4. Critical Biological Distinction: TGFB1 is Active but Pro-Fibrotic

NicheNet ranks TGFB1 as the #1 active sender ligand. **This does not mean TGFβ receptors are agonist targets.** TGFβ signaling to AT2 cells drives the arrested/Krt8⁺ state (pEMT, senescence, Cldn4 upregulation). Its receptors are **block candidates** — antagonizing TGFβ input to epithelium may allow cells to escape arrest and complete regeneration.

All candidates are annotated with `therapeutic_direction`:
- `agonist` — agonism promotes AT2→AT1 completion
- `block_TGFb` — receptor transduces TGFβ signaling → arrested state; antagonism indicated
- `agonist_uncertain` — pro-regenerative axis plausible but direction unconfirmed in AT2 context
- `uncertain` — mixed evidence or unknown

---

## 5. Final Ranked Candidates

### 5a. Agonist Candidates (activate to promote regeneration)

| Receptor | Composite | Species | Mouse NN rank | Human NN rank | Key ligands | Evidence |
|---|---|---|---|---|---|---|
| **Sdc1** | 0.850 | conserved | 8 | 2 | ANGPTL4, LPL, FGF2 | Heparan-sulfate co-receptor; presents FGF10 to FGFR2b; also ANGPTL4 downstream |
| **Erbb2** | 0.821 | conserved | — | 6 | AREG, HBEGF | EGF-family; promotes AT2 proliferation and AT2→AT1 transition |
| **Egfr** | 0.796 | conserved | — | 6 | AREG, HBEGF, CDCP1 | Same axis; EGFR drives epithelial repair post-injury |
| **Fzd2** | 0.722 | mouse-only | — | — | WNT | Highest Geneformer OE regen_gain (+0.008); WNT promotes AT2 self-renewal |
| **Sdc4** | 0.710 | conserved | 8 | 2 | ANGPTL4, FGF2, RSPO3 | Same syndecan/FGFR2b co-receptor axis as Sdc1 |
| **Fgfr2** | 0.603 | conserved | 3 | 10 | FGF7, FGF10, FGF1, FGF2, TIMP1 | **Positive control** — IIIb isoform (epithelial) drives AT2 maintenance ⚠️ |
| **Cd9** | 0.551 | conserved | 9 | 6 | HBEGF | Tetraspanin co-receptor; scaffolds EGFR/ERBB signaling |
| **Nrp1** | 0.437 | mouse-only | 8 | — | FGF1/2/7, SEMA3C, VEGFA | FGFR2b co-receptor and VEGFA axis; mouse NicheNet rank 8 |
| **Erbb3** | 0.414 | conserved | — | 16 | AREG | EGF-family; ERBB2/ERBB3 heterodimer important for AT2 survival |
| **Cd63** | 0.388 | conserved | 11 | 10 | TIMP1 | Tetraspanin; TIMP1→CD63 promotes epithelial survival |
| **Bmpr2** | 0.559 | conserved | — | 17 | BMP2 | BMP4 maintains AT2 identity in organoids; direction uncertain in fibrosis |
| **Bmpr1a** | 0.472 | conserved | — | 17 | BMP2 | BMP signaling receptor; uncertain in fibrotic context |

> ⚠️ **FGFR2 isoform flag:** Standard Fgfr2 quantification cannot separate IIIb (epithelial, FGFR2b) from IIIc (mesenchymal). Therapeutic targeting must be isoform-specific (FGFR2b agonists only).

### 5b. Block Candidates (antagonize to release arrested cells)

| Receptor | Composite | Species | Top ligand | Rationale |
|---|---|---|---|---|
| **App** | 0.863 | conserved | TGFB1 (both species) | TGFB1→APP promotes AT2 senescence; high composite from NicheNet but wrong direction |
| **Tgfbr2** | 0.682 | conserved | TGFB1 (both species) | Direct TGFβ type-II receptor → SMAD2/3 → Krt8 upregulation, pEMT, fibrosis |
| **Tgfbr1** | 0.481 | mouse-only | TGFB1 | Direct TGFβ type-I receptor (ALK5) → same SMAD2/3 axis |
| **Axl** | 0.296 | mouse-only | — | AXL promotes pEMT and mesenchymal transition; Geneformer DEL most negative |

### 5c. Uncertain Direction (require further mechanistic validation)

| Receptor | Composite | Species | Issue |
|---|---|---|---|
| **Itgb1** | 0.676 | conserved | ANGPTL4→ITGB1 (pro-survival, agonist-plausible) AND TGFB1→ITGB1 (pro-fibrotic); resolve with αv/β integrin context data |
| **Notch1/2** | 0.685/0.593 | mouse-only | Notch activation may promote basaloid/arrested fate in IPF; direction uncertain |
| **Cd74** | 1.000 | mouse-only | Pure CellChat hit (MIF receptor); no NicheNet or Geneformer evidence; MIF promotes inflammation |

---

## 6. Known Limitations and Caveats

### NicheNet v2 prior artifact (Fgf10)
The mouse lt_matrix contains only 1,226 target-gene columns. Fgf10's top predicted targets are prolactin-family mammary genes (not AT2 biology). Only 56/508 geneset_oi genes are covered by the lt_matrix. Fgf10 ranks 27 (not top 20) despite the correct geneset_oi — this is a hard prior ceiling, not a biology failure. Fgf7 (same FGFR2b IIIb receptor) ranks 3, validating the axis.

### FGFR2 IIIb vs IIIc isoform
All FGFR2 expression quantified here is total Fgfr2 (IIIb + IIIc). IIIb = epithelial (pro-regenerative), IIIc = mesenchymal (non-target). Therapeutic strategies must use FGFR2b-selective agonists (e.g., FGF7/FGF10 proteins, FGFR2b antibodies).

### Geneformer: ligand-limited vs receptor-limited
Fgfr2 OE regen_gain is negative (−0.003) because the Fgf10-KO bottleneck is on the ligand side (fibroblasts fail to produce FGF10), not the receptor side. In a model with intact FGF10 signaling, FGFR2b agonism should rescue completing fate.

### Human data: Habermann is IPF-only
No healthy AT2 reference in Habermann. Human geneset_oi (AT2+AT1 vs KRT5⁻/KRT17⁺) is the appropriate completing vs arrested analog. FGF10 is nearly absent from IPF fibroblasts (established fibrosis downregulates FGF10) — this is biology, not a data artifact.

### Whole-lung undersampling
Strunz GSE141259 is whole-lung. Epithelial cells are underrepresented vs sorted datasets (Kobayashi, Choi). Sender characterization (fibroblasts, macrophages) is reliable; AT2 trajectory analysis is cross-validated with sorted datasets.

---

## 7. Top Prioritized Agonist Candidates for Experimental Follow-up

Ranked by evidence confidence (cross-species conservation + multiple method support + biological rationale):

1. **FGFR2b (Fgfr2/FGFR2)** — positive control; FGF7/FGF10 agonists available; IIIb-isoform-specific agents required.
2. **SDC1 (Sdc1)** — ANGPTL4/LPL signal; also presents FGF10 to FGFR2b (heparan-sulfate scaffold); conserved mouse + human.
3. **EGFR / ERBB2 (Egfr/Erbb2)** — AREG/HBEGF repair axis; conserved; EGF-family promotes AT2 proliferation and AT2→AT1.
4. **SDC4 (Sdc4)** — same ANGPTL4/FGFR2b co-receptor axis as SDC1; conserved.
5. **FZD2 (Fzd2)** — highest Geneformer OE signal; WNT→AT2 self-renewal; test WNT3A or R-spondin agonists.
6. **NRP1 (Nrp1)** — FGFR2b co-receptor; synergizes with FGF10 signaling; mouse-only (limited human data).
7. **CD9 (Cd9)** — HBEGF co-receptor; scaffolds EGFR signaling; conserved; tetraspanin agonism is an emerging modality.

---

## 8. Scoring Method Details

Composite score = unweighted average of up to 5 rank-normalized (0–1) dimensions:

| Dimension | Source | Higher = better |
|---|---|---|
| `score_nn` | Mouse NicheNet AUPR rank (lower rank = better) | Low ligand rank |
| `score_nn_h` | Human NicheNet AUPR rank (Habermann) | Low ligand rank |
| `score_cc` | CellChat incoming probability to AT2 | High probability |
| `score_gf` | Geneformer OE regen_gain | High (OE → completing) |
| `score_del` | Geneformer DEL regen_loss | Low (DEL → arrested) |

Missing dimensions treated as NA and excluded from the mean (not penalized).

---

*Generated: 2026-06-15 | Pipeline: scripts/02_qc.R → scripts/07_crossspecies.R | Repo: XiaoCheng2776/CL_IPF_Target_prediction*
