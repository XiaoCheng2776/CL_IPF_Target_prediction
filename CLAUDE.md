# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project
scRNA-seq target-discovery pipeline for pulmonary fibrosis. Goal: rank cell-surface
receptors on alveolar epithelium (AT2 cells / injury-activated alveolar progenitors)
whose **agonism** promotes successful AT2 → AT1 regeneration and resolves fibrosis.

**FGF10 → FGFR2b is the POSITIVE CONTROL** (validated in vivo: Fgf10-deficient mice
develop fibrosis more readily). Full method is in
`docs/pulmonary_fibrosis_target_discovery_pipeline.txt` — read it first.

## CRITICAL: the biology is INVERTED vs a cardiac/epicardial pipeline
- EMT and proliferation are **disease/confounder** features here, **not** therapeutic goals.
- Define the target phenotype by **trajectory completion** (AT2 → Krt8⁺ transitional → AT1),
  NEVER by activation/EMT/proliferation magnitude.
- "Completing regeneration" = GOOD (reaches AT1). "Arrested in Krt8⁺/basaloid state" = BAD.
- In silico perturbation sign is flipped: OVEREXPRESS to find pro-regenerative receptors;
  DELETE to confirm direction. FGFR2 deletion pushing cells toward disease IS the
  Fgf10-KO mouse, in silico.

## Environment
- R 4.3+ : Seurat v5, harmony, nichenetr, CellChat, slingshot/monocle3, tidyverse
- Python 3.10+ : scanpy, geneformer, transformers, datasets, torch (GPU build)
- Geneformer fine-tuning requires a GPU with ≥16 GB VRAM — run that step on the
  GPU box/cluster; write & debug the script anywhere, then submit the job.

## Running the pipeline

```bash
# Step 1 — download raw data (idempotent; skips already-downloaded files)
bash scripts/01_download.sh

# Step 2-7 — one at a time; inspect output before proceeding
Rscript scripts/02_qc.R
Rscript scripts/03_phenotype.R
Rscript scripts/04_nichenet.R
Rscript scripts/05_cellchat.R
python scripts/06_geneformer.py          # GPU required for fine-tuning
Rscript scripts/07_crossspecies.R
```

> **Status:** Only `scripts/01_download.sh` exists. Scripts 02–07 are stubs to be
> scaffolded. Reference the pipeline doc's per-section code blocks when writing them.

## Data (accessions)
- Human: GSE136831 (Adams), GSE135893 (Habermann); healthy ref = HLCA via CELLxGENE
- Mouse: GSE141259 (Strunz, bleomycin time course d3–d28 — **primary**), GSE141634
  (Kobayashi PATS), GSE145031 (Choi DATPs)
- Raw data under `data/`, all outputs under `results/`. Both are git-ignored (`data/`,
  `results/`, `*.rds`, `.Rhistory`, `__pycache__/`).

## Output directories
`results/` already contains: `qc/`, `markers/`, `nichenet/`, `cellchat/`,
`geneformer/`, `figures/`. Geneformer sub-paths: `geneformer/oe/` (overexpression) and
`geneformer/del/` (deletion).

## Sender / Receiver framework (NicheNet & CellChat)
- **Senders** (niche signal source): alveolar (lipo)fibroblasts, airway smooth muscle
  cells, myofibroblasts, profibrotic macrophages
- **Receivers** (target population): AT2 cells and injury-activated alveolar progenitors

## Key marker gene sets

```r
# Mouse (lowercase); human = uppercase equivalents
at2_markers        <- c("Sftpc","Sftpb","Sftpa1","Lamp3","Etv5","Napsa","Slc34a2")
at1_markers        <- c("Ager","Pdpn","Hopx","Aqp5","Col4a3","Akap5")
transitional_mouse <- c("Krt8","Cldn4","Krt18","Sfn","Cdkn1a","Lgals3","Tpm1")   # Krt8+ ADI
basaloid_human     <- c("KRT17","TP63","LAMB3","LAMC2","ITGB6","MMP7","CDKN2A","SOX9","VIM","FN1")
                       # NOTE: KRT5 / KRT15 NEGATIVE in aberrant basaloid

# Signatures for module scoring
regen_success <- c(at2_markers, at1_markers)           # "completing" phenotype — GOOD
arrested      <- transitional_mouse                    # + senescence/pEMT confounders
```

## Conventions
- Mouse genes lowercase (`Fgf10`, `Fgfr2`); human uppercase (`FGF10`, `FGFR2`).
- Epithelium is fragile — keep QC permissive: nFeature_RNA 200–7000, percent.mt < 20%.
- Mito pattern: mouse `^mt-`, human `^MT-`.
- Figures: PNG, 300 dpi, into `results/figures/`.
- Scripts are numbered `01_…`→`07_…` under `scripts/`; keep them idempotent.

## Guardrails
- **Positive-control gate:** if FGF10/Fgf10 does NOT rank near the top in NicheNet,
  STOP and flag it — debug thresholds / geneset_oi, do not proceed silently.
- Never auto-delete anything under `data/`. Ask before any compute step >30 min.
- Pin package versions. Show me diffs before accepting edits.
- Flag the FGFR2 IIIb (epithelial) vs IIIc (mesenchymal) isoform issue wherever
  FGFR2 expression is quantified.
- Whole-lung dissociation undersamples epithelium; prefer sorted/epithelial-enriched
  data (Kobayashi/Choi) for the receiver trajectory; Strunz whole-lung is best for
  sender characterization.
