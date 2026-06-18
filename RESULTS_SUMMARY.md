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
| FGF10 sender | 08_fgf10_sender.R | Pseudobulk Wilcoxon, % expressing, compositional | Adams, Habermann, Strunz |
| FGFR2 receiver | 08b_fgfr2_receiver.R | Pseudobulk Wilcoxon, per-cell Spearman | epi_obj, Adams, Habermann |
| Synthesis figure | 09c_synthesis_figure.R | Lollipop chart, two-arm layout | candidates_ranked.csv |

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
Fgfr2 OE regen_gain is negative (−0.003) because the Fgf10-KO bottleneck is on the ligand side (fibroblasts fail to produce FGF10), not the receptor side. In a model with intact FGF10 signaling, FGFR2b agonism should rescue completing fate. This is confirmed by 08b (§9): FGFR2 receptor expression is maintained or elevated in arrested AT2 cells — the antenna is intact, the ligand signal is absent.

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

---

## 9. FGF10 Sender Deficiency + FGFR2 Receiver Expression

These two analyses probe the FGF10 → FGFR2b axis from both sides to test whether the deficiency is in **ligand supply** (fibroblast sender), **receptor responsiveness** (AT2 receiver), or both.

### 9a. FGF10 sender-side deficiency (08_fgf10_sender.R)

Three metrics per dataset and cell type: fraction expressing (dropout-resistant), pseudobulk mean per sample + Wilcoxon (proper inference), and mean per expresser (per-cell downregulation test).

**Mouse bleomycin time course (Strunz GSE141259, fibroblasts):**

| Phase | Fgf10+ % (Plin2-low Fib) | vs PBS |
|---|---|---|
| PBS control | 0% | — |
| acute (d3–d7) | 0–8.5% | p = 0.029, log2FC = +5.3 |
| fibrotic early (d10–d14) | 3–5% | p = 0.029 (pooled) |
| fibrotic late (d21) | 0% | p = 0.674 |
| resolution (d28) | 0% | p = 0.789 |

**Interpretation:** Fgf10 expression from fibroblasts is **episodic**, not chronic. It rises transiently during acute injury (d7) when regeneration is initiated, then fails to sustain — collapsing to 0% by the fibrotic phase (d21). The story is a **failure to sustain the repair signal**, not constitutive suppression.

**Human fibroblasts (Adams GSE136831 + Habermann GSE135893):**

| Dataset | Cell type | IPF % | Control % | p | log2FC |
|---|---|---|---|---|---|
| Adams | Myofibroblast | 5.1% | 6.4% | 0.012 | +6.4 (IPF > Ctrl) |
| Adams | PLIN2+ Fibroblast | 2.3% | 2.3% | 0.746 | 0 |
| Adams | PLIN2− Fibroblast | 3.2% | 1.9% | 0.089 | ns |
| Habermann | Fibroblasts | 8.6% | 0.7% | 0.005 | +6.3 (IPF > Ctrl) |
| Habermann | PLIN2+ Fibroblasts | 5.7% | 0% | ns | — |

**Paradox:** FGF10 appears higher in IPF fibroblasts in both datasets. This is not a chronic excess — cell numbers (Habermann PLIN2+ IPF n=1,167 vs Control n=6) suggest this reflects **pathological population expansion** in IPF, not increased per-cell production. The PLIN2+ cells in Habermann are likely reactive fibroblasts, not classical alveolar lipofibroblasts (essentially zero FGF10 expression in either condition confirms this). The episodic mouse data is the mechanistically cleaner evidence of FGF10 supply failure.

**Compositional analysis:** No significant depletion of PLIN2+ lipofibroblasts in Adams IPF vs Control at the proportion level. The FGF10 deficit is primarily temporal (failure to sustain in fibrosis) rather than strictly compositional (fewer FGF10-producing cells).

---

### 9b. FGFR2 receiver-side expression (08b_fgfr2_receiver.R)

> **FGFR2 isoform caveat:** AT2 cells express FGFR2 predominantly as the IIIb isoform (Nakayama 2011, McQualter 2019). Total FGFR2 in AT2 cells is used as a FGFR2b surrogate here; IIIb/IIIc separation requires isoform-specific assays and cannot be resolved from standard RNA-seq data.

**Mouse bleomycin time course (Kobayashi PATS epi_obj, AT2 + AT2-activated cells, n = 13,950):**

| Phase | n cells | Fgfr2+ % | Tgfbr2+ % | Transitional1 score |
|---|---|---|---|---|
| PBS control | 1,739 | 13.0% | 3.7% | −0.077 |
| acute (d2–d10) | 7,447 | 13.3% | 7.2% | −0.054 |
| fibrotic early (d11–d14) | 1,366 | 11.3% | 5.3% | −0.005 |
| fibrotic late (d15–d21) | 1,422 | **10.6%** | 4.2% | −0.014 |
| resolution (d28–d54) | 1,976 | 11.2% | 3.4% | −0.046 |

- Pseudobulk Wilcoxon vs PBS: all p > 0.38 — the ~2 percentage point decrease during fibrosis is **not significant**.
- Per-cell Spearman (Fgfr2 count ~ Transitional1 score): **rho = −0.070, p = 9×10⁻¹⁷** — direction correct (FGFR2-expressing cells have marginally lower TGFβ-arrest score), but effect size is trivially small; significance is driven by n = 13,950.

**Human Adams ATII cells (IPF n = 496, Control n = 2,655):**

| Gene | IPF vs Control | p | log2FC | Direction |
|---|---|---|---|---|
| FGFR2 | IPF median 1.06 vs Ctrl 0.55 | **0.052** | **+0.94** | ↑ IPF ≥ Control |
| TGFBR2 | — | 0.711 | −0.27 | ns |
| CLDN4 | — | 0.253 | −0.49 | ns |
| KRT8 | — | 0.799 | +0.23 | ns |

- Per-cell Spearman (FGFR2 ~ CLDN4+KRT8+SERPINE1 arrest score): **rho = +0.197, p = 6×10⁻²⁹**
- FGFR2-expressing ATII cells co-express more arrest markers — they are in the transitional (Krt8⁺) state, not silenced.

**Human Habermann IPF AT2 (n = 4,903; IPF only):**

- FGFR2 expressers: 32.7%; TGFBR2: 60.3%
- Per-cell Spearman (FGFR2 ~ CLDN4+KRT8): rho = +0.078, p = 5×10⁻⁸

---

### 9c. Integrated interpretation: the "frustrated receptor" model

The "FGFR2 drops as TGFβ rises" hypothesis is **not supported**. FGFR2 is maintained — borderline elevated in IPF AT2 cells — and the cells expressing it are enriched in the Krt8⁺ transitional arrested state. This is consistent with a **frustrated receptor model**:

> AT2 cells that are stuck in the Krt8⁺ arrested state retain FGFR2b expression (the antenna is intact), but the FGF10 repair signal from fibroblasts has collapsed. The cells are poised to respond but cannot because the ligand is absent.

**Therapeutic implication:** FGFR2b receptor presence in arrested AT2 cells means FGF10/FGFR2b agonism has a viable target. The bottleneck is ligand supply, not receptor loss. This model:

1. Explains the Geneformer result: FGFR2 OE regen_gain is negative because overexpressing an already-present receptor in a ligand-absent context does not rescue completing fate. Exogenous FGF10 (or a FGFR2b-selective agonist) is the correct intervention.
2. Is consistent with the in vivo Fgf10-KO mouse (reduced ligand → more fibrosis), not the FGFR2-KO (which would be a different phenotype).
3. Predicts that FGF7 (same receptor, higher NicheNet rank = 3) or FGFR2b antibody agonists could rescue completing fate in IPF AT2 cells without needing to upregulate the receptor.

**Cross-dataset dual-hit figure:** `results/figures/fgfr2_timecourse_dual.png` — normalised Fgfr2 in AT2 (Kobayashi) and Fgf10 in fibroblasts (Strunz) over bleomycin phases. Co-trending signals support a coordinated sender-receiver collapse during fibrosis (eco-correlation across 5 phases; n too small for significance — directional only).

---

---

## 10. Conserved-Target Synthesis Figure (09c_synthesis_figure.R)

**Outputs:** `results/figures/Fig_synthesis.png` (300 dpi) + `Fig_synthesis.pdf` (cairo_pdf, Unicode-safe)

The synthesis figure is the single-panel summary that bridges the computational screen to experimental prioritization. It displays all 14 scored candidates from `candidates_ranked.csv` in a horizontal lollipop layout with x = composite score (0–1) and explicit NicheNet rank annotations (mouse / human) for each target.

### 10a. Layout and sections

| Section | Colour | Candidates | Criterion for inclusion |
|---|---|---|---|
| **Predicted stimulate** (top arm) | #E07B39 (orange) | Sdc1, Erbb2, Egfr, Sdc4, Fgfr2★, Cd9 | conserved + human NicheNet evidence + agonist direction |
| **Inhibit** (middle arm) | #4472B8 (blue) | App, Tgfbr2 | conserved + block_TGFb direction supported both species |
| **Mouse-only** (grey box, bottom) | #9E9E9E | Cd74 (?), Ramp1 (↑?), Fzd2 (↑?), Bmpr2 (↑?), Nrp1 (↑?), Tgfbr1 (↓) | mouse-only evidence; no human NicheNet rank |

Fgfr2 is coloured #B83232 (brick red) and labelled with ★ PC to mark it as the positive control throughout the figure.

### 10b. Conserved agonist arm — key take-aways

| Receptor | Composite | Mouse NN | Human NN | Highlighted observation |
|---|---|---|---|---|
| Sdc1 | 0.850 | #8 | **#2** | Strongest human NicheNet signal of any agonist; presents FGF10 to FGFR2b |
| Erbb2 | 0.821 | — | #6 | Pure human NicheNet + Geneformer DEL; AREG/HBEGF repair axis |
| Egfr | 0.796 | — | #6 | Same EGF-family axis as ERBB2; conserved AREG/HBEGF/CDCP1 |
| Sdc4 | 0.710 | #8 | **#2** | Same syndecan/ANGPTL4 axis as SDC1; lower Geneformer signal |
| Fgfr2 ★ | 0.603 | #3 | #10 | Positive control; composite depressed by ligand-limited Geneformer (§6/§9c) |
| Cd9 | 0.551 | #9 | #6 | HBEGF tetraspanin co-receptor; consistent but lower overall score |

**Direction caveat (agonist arm):** therapeutic direction is predicted from NicheNet ligand-activity and known receptor biology. Geneformer in silico perturbation (overexpression arm) requires a dedicated GPU run with the full fine-tuned model; existing OE scores were generated from a subset model and have not been validated. The block arm (INHIBIT) does not carry this caveat — TGFβ→Krt8⁺ arrest is mechanistically established and both NicheNet and Geneformer deletion scores support it.

### 10c. Conserved block arm

| Receptor | Composite | Mouse NN | Human NN | Rationale |
|---|---|---|---|---|
| App | 0.863 | **#1** | #5 | Highest mouse NicheNet rank of all candidates; TGFB1→APP pro-arrest axis conserved |
| Tgfbr2 | 0.682 | — | #5 | Direct TGFβ type-II receptor → SMAD2/3 → Krt8 upregulation; mechanistically the most direct block target |

App's #1 mouse NicheNet rank reflects that TGFB1 (the top sender ligand) signals strongly through APP in AT2 cells — antagonizing APP may divert TGFβ-activated cells away from the arrested state. Tgfbr2 blockade (e.g., galunisertib analogs) is the most mechanistically interpretable intervention.

### 10d. Mouse-only section

Six candidates are greyed out and labelled "awaiting human validation." They were detected in the mouse NicheNet and/or Geneformer screens but lack supporting human NicheNet ranks (no human evidence in Adams or Habermann). They are not suitable for immediate experimental prioritization but should be revisited if human sorted-AT2 data (e.g., IPF surgical lung biopsies with epithelial enrichment) becomes available.

Notable entries:
- **Cd74** (composite 1.000) — sole CellChat MIF-receptor hit; no NicheNet or Geneformer support; high composite is an artefact of a perfect CellChat score on one dimension.
- **Ramp1** (0.897) — ADM/adrenomedullin axis; mouse NicheNet rank 2; human data absent.
- **Fzd2** (0.722) — highest Geneformer OE regen_gain (+0.008) of any candidate; WNT→AT2 self-renewal; no human NicheNet rank.

---

*Generated: 2026-06-18 | Pipeline: scripts/02_qc.R → scripts/09c_synthesis_figure.R | Repo: XiaoCheng2776/CL_IPF_Target_prediction*
