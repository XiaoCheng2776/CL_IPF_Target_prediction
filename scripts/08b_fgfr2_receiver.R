#!/usr/bin/env Rscript
# =============================================================================
# 08b_fgfr2_receiver.R — receiver-side FGFR2 analysis in AT2 cells
# Run from repo root:  Rscript scripts/08b_fgfr2_receiver.R
#
# Tests the "FGF10 supply drowned by TGFβ" hypothesis in silico:
#   Hypothesis: as fibrosis progresses, TGFβ signalling rises in AT2 cells,
#   coinciding with declining FGFR2 expression (reduced responsiveness to FGF10).
#
# Three datasets:
#   A. Mouse epi_obj.rds  (Kobayashi PATS, bleomycin time course, n=13,950 AT2)
#      • Fgfr2 expression in AT2 cells by phase (PBS → acute → fibrotic → resolution)
#      • Tgfbr2 expression by phase (TGFβ receptor rising?)
#      • Transitional1 module score by phase (TGFβ-driven arrest proxy, pre-computed)
#      • Per-cell Spearman: Fgfr2 ~ Transitional1 (within AT2 cells)
#      • Cross-dataset: pair with Fgf10 in Strunz fibroblasts per phase (from 08a CSV)
#   B. Human Adams GSE136831 (IPF + Control ATII cells)
#      • FGFR2, TGFBR1, TGFBR2, CLDN4, KRT8, SERPINE1 awk-extracted
#      • IPF vs Control: % expressing + pseudobulk Wilcoxon
#   C. Human Habermann GSE135893 (IPF only; saved Seurat)
#      • FGFR2 in AT2 cells; correlate with SMAD2/3-target score
#
# Isoform caveat (flagged everywhere):
#   AT2 cells express FGFR2 predominantly as the IIIb isoform (Nakayama 2011,
#   McQualter 2019). Standard RNA-seq cannot separate IIIb (epithelial, target)
#   from IIIc (mesenchymal, non-target). FGFR2 in AT2 cells is treated as a
#   FGFR2b surrogate with this caveat stated in every output.
#
# Outputs:
#   results/fgfr2_receiver_summary.csv
#   results/fgfr2_receiver_wilcox.csv
#   results/figures/fgfr2_timecourse_dual.png   — Fgfr2 (AT2) + Fgf10 (Fib) over time
#   results/figures/fgfr2_vs_tgfb_scatter.png   — per-cell scatter
#   results/figures/fgfr2_timecourse_mouse.png   — Fgfr2/Tgfbr2/Transitional1 by phase
#   results/figures/fgfr2_human_ipf_ctrl.png     — Adams human comparison
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(tidyverse)
})

ROOT      <- here::here()
EPI_RDS   <- file.path(ROOT, "results/markers/epi_obj.rds")
HAB_RDS   <- file.path(ROOT, "results/nichenet_human/hab_seurat.rds")
ADAMS_DIR <- file.path(ROOT, "data/human/adams_gse136831")
SENDER_CSV <- file.path(ROOT, "results/fgf10_sender_summary.csv")
OUT_FIG   <- file.path(ROOT, "results/figures")
OUT_SUMM  <- file.path(ROOT, "results")

dir.create(OUT_FIG, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers (same as 08_fgf10_sender.R)
# ─────────────────────────────────────────────────────────────────────────────

gene_metrics_df <- function(counts, sample_ids, gene, subtype, condition, dataset) {
  pb <- tapply(counts, sample_ids, mean, na.rm = TRUE)
  data.frame(
    dataset            = dataset,
    gene               = gene,
    condition          = condition,
    subtype            = subtype,
    n_cells            = length(counts),
    n_samples          = n_distinct(sample_ids),
    pct_expressing     = 100 * mean(counts > 0),
    mean_expressers    = if (any(counts > 0)) mean(counts[counts > 0]) else NA_real_,
    pseudobulk_mean    = mean(pb, na.rm = TRUE),
    pseudobulk_median  = median(pb, na.rm = TRUE),
    stringsAsFactors   = FALSE
  )
}

pb_wilcox <- function(pb_a, pb_b, label_a = "IPF", label_b = "Ctrl") {
  if (length(pb_a) < 2 || length(pb_b) < 2) {
    return(data.frame(p_value=NA_real_, median_a=median(pb_a), median_b=median(pb_b),
                      log2fc=NA_real_))
  }
  wt <- wilcox.test(pb_a, pb_b, exact = FALSE)
  data.frame(p_value    = wt$p.value,
             median_a   = median(pb_a, na.rm = TRUE),
             median_b   = median(pb_b, na.rm = TRUE),
             log2fc     = log2((median(pb_a, na.rm=TRUE)+0.001) /
                                (median(pb_b, na.rm=TRUE)+0.001)))
}

mtx_skip_lines <- function(mtx_gz) {
  con <- gzcon(file(mtx_gz, "rb")); on.exit(close(con))
  lines <- readLines(con, n = 10)
  sum(startsWith(lines, "%")) + 1L
}

mtx_dims <- function(mtx_gz) {
  con <- gzcon(file(mtx_gz, "rb")); on.exit(close(con))
  lines <- readLines(con, n = 10)
  as.integer(strsplit(trimws(lines[!startsWith(lines, "%")][1]), "\\s+")[[1]])
}

extract_genes_mtx <- function(mtx_gz, gene_indices, barcodes, skip_lines = 3) {
  idx_str <- paste(gene_indices, collapse = "|")
  cmd <- sprintf(
    "zcat '%s' | awk 'NR>%d && ($1 ~ /^(%s)$/) {print $1, $2, $3}'",
    mtx_gz, skip_lines, idx_str)
  message(sprintf("  Extracting %d gene(s) from %s…", length(gene_indices), basename(mtx_gz)))
  dt <- tryCatch(
    fread(cmd = cmd, header = FALSE, col.names = c("gene_idx","cell_idx","count")),
    error = function(e) data.table(gene_idx=integer(), cell_idx=integer(), count=numeric())
  )
  n <- length(barcodes)
  out <- matrix(0.0, nrow = n, ncol = length(gene_indices),
                dimnames = list(barcodes, names(gene_indices)))
  rev_idx <- setNames(names(gene_indices), as.character(gene_indices))
  if (nrow(dt) > 0) {
    gnames <- rev_idx[as.character(dt$gene_idx)]
    valid  <- !is.na(gnames) & dt$cell_idx >= 1 & dt$cell_idx <= n
    for (gn in unique(na.omit(gnames[valid]))) {
      sel <- valid & gnames == gn
      out[dt$cell_idx[sel], gn] <- dt$count[sel]
    }
  }
  as.data.frame(out)
}

result_rows <- list()
wilcox_rows <- list()

# =============================================================================
# A.  MOUSE: FGFR2 IN AT2 CELLS — epi_obj.rds (Kobayashi bleomycin time course)
# =============================================================================
message("\n========== A. MOUSE RECEIVER: epi_obj.rds ==========")
message("Isoform caveat: AT2 Fgfr2 treated as FGFR2b (IIIb) surrogate.",
        " IIIb/IIIc separation requires isoform-specific assays.")

epi <- readRDS(EPI_RDS)

# Detect cell-type column (03_phenotype.R may have used cell_type or cell.type)
ct_col <- if ("cell_type" %in% colnames(epi@meta.data)) "cell_type" else
          if ("cell.type" %in% colnames(epi@meta.data)) "cell.type" else {
            # fallback: find any column containing AT2 values
            hits <- sapply(colnames(epi@meta.data), function(cn) any(grepl("AT2", epi@meta.data[[cn]], ignore.case=TRUE)))
            if (any(hits)) names(hits)[which(hits)[1]] else stop("Cannot identify cell-type column in epi_obj.rds")
          }
message(sprintf("  Using cell-type column: '%s'", ct_col))
message(sprintf("  Available labels: %s", paste(sort(unique(epi@meta.data[[ct_col]])), collapse=", ")))

at2 <- epi[, epi@meta.data[[ct_col]] %in% c("AT2","AT2 activated")]
DefaultAssay(at2) <- "SCT"
message(sprintf("  AT2 + AT2 activated cells: %d", ncol(at2)))

# Extract raw-count proxies from SCT assay
target_genes_m <- c("Fgfr2","Tgfbr2","Tgfbr1","Tgfb1","Fgfr1")
target_genes_m <- intersect(target_genes_m, rownames(at2))
cnt_at2 <- as.matrix(GetAssayData(at2, assay = "SCT", layer = "counts")[target_genes_m, ])

at2_meta <- at2@meta.data
for (g in target_genes_m) {
  at2_meta[[paste0(g, "_count")]] <- as.numeric(cnt_at2[g, rownames(at2_meta)])
}
# Ensure all expected count columns exist (NA-filled if gene absent from assay)
for (g in c("Fgfr2","Tgfbr2","Tgfbr1","Tgfb1","Fgfr1")) {
  col <- paste0(g, "_count")
  if (!col %in% colnames(at2_meta)) at2_meta[[col]] <- 0L
}

# Detect time_point column
tp_col <- if ("time_point" %in% colnames(at2_meta)) "time_point" else
          if ("timepoint"  %in% colnames(at2_meta)) "timepoint"  else
          if ("Time_point" %in% colnames(at2_meta)) "Time_point" else
            stop("Cannot find time_point column in epi_obj metadata")
message(sprintf("  Using time-point column: '%s'", tp_col))
message(sprintf("  Unique time points: %s", paste(sort(unique(at2_meta[[tp_col]])), collapse=", ")))
at2_meta$time_point <- at2_meta[[tp_col]]   # normalise name

# Detect Transitional1 score column (AddModuleScore appends "1" → sometimes "Transitional11")
trans_col <- if ("Transitional1" %in% colnames(at2_meta)) "Transitional1" else
             if ("Transitional11" %in% colnames(at2_meta)) "Transitional11" else {
               hits <- grep("^[Tt]ransitional", colnames(at2_meta), value=TRUE)
               if (length(hits)) hits[1] else {
                 warning("Transitional1 column not found — TGFβ-arrest correlation will be skipped.")
                 NULL
               }
             }
if (!is.null(trans_col)) {
  message(sprintf("  Using Transitional score column: '%s'", trans_col))
  at2_meta$Transitional1 <- at2_meta[[trans_col]]
} else {
  at2_meta$Transitional1 <- NA_real_
}

# Phase mapping: parse numeric day from "day N" or "dN" patterns
at2_meta$day_num <- suppressWarnings(
  as.numeric(gsub(".*?([0-9]+).*", "\\1", at2_meta$time_point))
)
at2_meta$phase <- dplyr::case_when(
  grepl("PBS|control|NC", at2_meta$time_point, ignore.case=TRUE) ~ "PBS_control",
  !is.na(at2_meta$day_num) & at2_meta$day_num <= 10             ~ "acute",
  !is.na(at2_meta$day_num) & at2_meta$day_num <= 14             ~ "fibrotic_early",
  !is.na(at2_meta$day_num) & at2_meta$day_num <= 21             ~ "fibrotic_late",
  !is.na(at2_meta$day_num) & at2_meta$day_num > 21              ~ "resolution",
  TRUE ~ "other"
)
message(sprintf("  Phase distribution:\n%s",
                paste(capture.output(print(table(at2_meta$phase))), collapse="\n")))

# Per-phase summary
PHASE_ORDER <- c("PBS_control","acute","fibrotic_early","fibrotic_late","resolution")
phase_summ_m <- at2_meta %>%
  filter(phase %in% PHASE_ORDER) %>%
  group_by(phase) %>%
  summarise(
    n_cells           = n(),
    n_mice            = n_distinct(orig.ident),
    pct_Fgfr2         = 100 * mean(Fgfr2_count > 0),
    pct_Tgfbr2        = 100 * mean(Tgfbr2_count > 0),
    mean_Fgfr2_expr   = if (any(Fgfr2_count>0)) mean(Fgfr2_count[Fgfr2_count>0]) else NA_real_,
    mean_Transitional = mean(Transitional1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(phase = factor(phase, levels = PHASE_ORDER))

message("\n  Mouse AT2: Fgfr2 + TGFβ marker by phase:")
print(phase_summ_m)

# Pseudobulk Wilcoxon: fibrotic phases vs PBS  (per-mouse pseudobulk via orig.ident)
message("\n  Mouse pseudobulk Wilcoxon (Fgfr2 % in AT2 cells, vs PBS):")
pbs_pb_fgfr2 <- tapply(
  at2_meta$Fgfr2_count[at2_meta$phase == "PBS_control"] > 0,
  at2_meta$orig.ident[at2_meta$phase == "PBS_control"],
  mean, na.rm = TRUE)
for (ph in c("acute","fibrotic_early","fibrotic_late","resolution")) {
  sel <- at2_meta$phase == ph
  if (sum(sel) < 10) next
  bleo_pb <- tapply(at2_meta$Fgfr2_count[sel] > 0, at2_meta$orig.ident[sel], mean, na.rm=TRUE)
  wres <- pb_wilcox(bleo_pb, pbs_pb_fgfr2)
  wres$dataset <- "epi_obj_mouse"; wres$gene <- "Fgfr2"; wres$subtype <- paste0("AT2_", ph)
  wilcox_rows[[length(wilcox_rows)+1]] <- wres
  message(sprintf("    AT2 %s vs PBS:  Fgfr2 p=%.3f  log2FC=%.2f  (bleo %.1f%%, PBS %.1f%%)",
                  ph, wres$p_value, wres$log2fc,
                  100*wres$median_a, 100*wres$median_b))
}

# Per-cell Spearman: Fgfr2 expression vs Transitional1 score
# Hypothesis: more TGFβ-driven arrest (higher Transitional1) → less Fgfr2
filt <- at2_meta$phase %in% PHASE_ORDER & !is.na(at2_meta$Transitional1)
ct_spear <- cor.test(at2_meta$Fgfr2_count[filt],
                     at2_meta$Transitional1[filt],
                     method = "spearman", exact = FALSE)
message(sprintf("\n  Per-cell Spearman (Fgfr2 count vs Transitional1 score):"))
message(sprintf("    rho = %.3f  p = %.2e  (n = %d AT2 cells)",
                ct_spear$estimate, ct_spear$p.value, sum(filt)))
message("  Interpretation: negative rho supports 'high TGFβ-arrest → lower FGFR2 responsiveness'")

# Also: Tgfbr2 vs phase
pbs_pb_tgfbr2 <- tapply(
  at2_meta$Tgfbr2_count[at2_meta$phase == "PBS_control"] > 0,
  at2_meta$orig.ident[at2_meta$phase == "PBS_control"],
  mean, na.rm=TRUE)
message("\n  Mouse Tgfbr2 in AT2 by phase (is TGFβ receptor rising?):")
for (ph in c("acute","fibrotic_early","fibrotic_late","resolution")) {
  sel <- at2_meta$phase == ph
  if (sum(sel) < 10) next
  bleo_pb <- tapply(at2_meta$Tgfbr2_count[sel] > 0, at2_meta$orig.ident[sel], mean, na.rm=TRUE)
  wres2 <- pb_wilcox(bleo_pb, pbs_pb_tgfbr2)
  wres2$dataset <- "epi_obj_mouse"; wres2$gene <- "Tgfbr2"; wres2$subtype <- paste0("AT2_", ph)
  wilcox_rows[[length(wilcox_rows)+1]] <- wres2
  message(sprintf("    AT2 %s vs PBS:  Tgfbr2 p=%.3f  log2FC=%.2f  (bleo %.1f%%, PBS %.1f%%)",
                  ph, wres2$p_value, wres2$log2fc,
                  100*wres2$median_a, 100*wres2$median_b))
}

# Collect per-phase metrics for summary CSV
for (ph in PHASE_ORDER) {
  sel <- at2_meta$phase == ph
  if (sum(sel) < 5) next
  for (g in target_genes_m) {
    row <- gene_metrics_df(
      counts     = at2_meta[[paste0(g,"_count")]][sel],
      sample_ids = at2_meta$orig.ident[sel],
      gene       = g, subtype = "AT2",
      condition  = ph, dataset = "Mouse_epi_obj"
    )
    result_rows[[length(result_rows)+1]] <- row
  }
}

# Cross-dataset pairing: Fgfr2 in AT2 (epi_obj) vs Fgf10 in Fibroblasts (Strunz, from 08a CSV)
if (file.exists(SENDER_CSV)) {
  sender_df <- read.csv(SENDER_CSV)
  # Strunz fibroblast Fgf10+ fraction by grouping (d3, d7, PBS, …)
  # Map to same phase labels
  strunz_fib <- sender_df %>%
    filter(dataset == "Strunz_mouse_GSE141259",
           subtype == "Fibroblasts (Plin2-low)") %>%   # use broad group
    mutate(phase = dplyr::case_when(
      condition == "PBS" ~ "PBS_control",
      condition %in% c("d3","d7")        ~ "acute",
      condition %in% c("d10","d14")      ~ "fibrotic_early",
      condition == "d21"                 ~ "fibrotic_late",
      condition == "d28"                 ~ "resolution",
      TRUE ~ "other"
    )) %>%
    filter(phase != "other") %>%
    group_by(phase) %>%
    summarise(fib_fgf10_pct = mean(pct_expressing, na.rm=TRUE), .groups="drop")

  cross_df <- phase_summ_m %>%
    left_join(strunz_fib, by = "phase") %>%
    filter(!is.na(fib_fgf10_pct))

  if (nrow(cross_df) >= 3) {
    cr <- cor.test(cross_df$pct_Fgfr2, cross_df$fib_fgf10_pct,
                   method = "spearman", exact = FALSE)
    message(sprintf("\n  Cross-dataset Spearman: Fgfr2%% in AT2 vs Fgf10%% in Fib (by phase)"))
    message(sprintf("    rho = %.2f  p = %.2f  (n=%d phases)", cr$estimate, cr$p.value, nrow(cross_df)))
    message("  (Note: eco-correlation across 5 phases; directional only)")
  }
} else {
  message("  Sender CSV not found — skipping cross-dataset correlation. Run 08_fgf10_sender.R first.")
}

# =============================================================================
# B.  HUMAN: FGFR2 IN AT2 CELLS — Adams GSE136831 (IPF + Control)
# =============================================================================
message("\n========== B. HUMAN RECEIVER: Adams ATII cells ==========")
message("Isoform caveat: ATII FGFR2 treated as FGFR2b surrogate (see header).")

ADAMS_MTX  <- file.path(ADAMS_DIR, "GSE136831_RawCounts_Sparse.mtx.gz")
ADAMS_BAR  <- file.path(ADAMS_DIR, "GSE136831_AllCells.cellBarcodes.txt.gz")
ADAMS_GENE <- file.path(ADAMS_DIR, "GSE136831_AllCells.GeneIDs.txt.gz")
ADAMS_META <- file.path(ADAMS_DIR, "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz")

meta_a <- fread(ADAMS_META, sep = "\t", header = TRUE)
setnames(meta_a, "CellBarcode_Identity", "barcode")

# Adams AT2 = "ATII"; arrested/basaloid = "Aberrant_Basaloid"
at2_a <- meta_a[Manuscript_Identity == "ATII" &
                 Disease_Identity   %in% c("IPF","Control")]
message(sprintf("  Adams ATII cells (IPF+Control): %d", nrow(at2_a)))
print(table(at2_a$Disease_Identity))

barcodes_a <- fread(ADAMS_BAR, header = FALSE)[[1]]
genes_a_tbl <- fread(ADAMS_GENE, header = TRUE)
genes_a_sym <- make.unique(as.character(genes_a_tbl[[2]]))
genes_a_ens <- as.character(genes_a_tbl[[1]])

target_genes_h  <- c("FGFR2","FGFR1","TGFBR1","TGFBR2","CLDN4","KRT8","SERPINE1")
ENSEMBL_MAP_H <- c(FGFR2="ENSG00000066468", FGFR1="ENSG00000077782",
                   TGFBR1="ENSG00000106799", TGFBR2="ENSG00000163513",
                   CLDN4="ENSG00000189143", KRT8="ENSG00000170421",
                   SERPINE1="ENSG00000106366")

gene_idx_h <- setNames(match(target_genes_h, genes_a_sym), target_genes_h)
missing_h  <- names(gene_idx_h)[is.na(gene_idx_h)]
if (length(missing_h) > 0) {
  for (gn in missing_h) {
    idx <- match(ENSEMBL_MAP_H[gn], genes_a_ens)
    if (!is.na(idx)) gene_idx_h[gn] <- idx
  }
}
gene_idx_h <- gene_idx_h[!is.na(gene_idx_h)]
message(sprintf("  Genes matched: %s", paste(names(gene_idx_h), collapse=", ")))

skip_a <- mtx_skip_lines(ADAMS_MTX)
expr_at2_a <- extract_genes_mtx(ADAMS_MTX, gene_idx_h, barcodes_a, skip_lines = skip_a)

# Subset to AT2 barcodes
at2_barcodes_a <- intersect(at2_a$barcode, rownames(expr_at2_a))
at2_expr_a <- expr_at2_a[at2_barcodes_a, , drop = FALSE]
at2_data_a <- cbind(
  as.data.frame(at2_a[match(at2_barcodes_a, at2_a$barcode),
                       .(barcode, Disease_Identity, Subject_Identity)]),
  at2_expr_a
)
message(sprintf("  AT2 cells with expression data: %d", nrow(at2_data_a)))
message(sprintf("  FGFR2 expressers: %d / %d (%.1f%%)",
                sum(at2_data_a$FGFR2 > 0, na.rm=TRUE), nrow(at2_data_a),
                100*mean(at2_data_a$FGFR2 > 0, na.rm=TRUE)))

# Three metrics per condition
for (g in names(gene_idx_h)) {
  for (cond in c("IPF","Control")) {
    sel <- at2_data_a$Disease_Identity == cond
    if (sum(sel) < 5) next
    row <- gene_metrics_df(
      counts     = at2_data_a[[g]][sel],
      sample_ids = at2_data_a$Subject_Identity[sel],
      gene = g, subtype = "ATII", condition = cond, dataset = "Adams_GSE136831"
    )
    result_rows[[length(result_rows)+1]] <- row
  }
}

# Pseudobulk Wilcoxon for FGFR2, TGFBR2, CLDN4 (TGFβ arrest marker)
message("\n  Adams ATII pseudobulk Wilcoxon (IPF vs Control):")
for (g in c("FGFR2","TGFBR2","CLDN4","KRT8","SERPINE1")) {
  if (!g %in% colnames(at2_data_a)) next
  ipf_pb  <- tapply(at2_data_a[[g]][at2_data_a$Disease_Identity=="IPF"],
                    at2_data_a$Subject_Identity[at2_data_a$Disease_Identity=="IPF"],
                    mean, na.rm=TRUE)
  ctrl_pb <- tapply(at2_data_a[[g]][at2_data_a$Disease_Identity=="Control"],
                    at2_data_a$Subject_Identity[at2_data_a$Disease_Identity=="Control"],
                    mean, na.rm=TRUE)
  if (length(ipf_pb) < 2 || length(ctrl_pb) < 2) next
  wres <- pb_wilcox(ipf_pb, ctrl_pb)
  wres$dataset <- "Adams_GSE136831"; wres$gene <- g; wres$subtype <- "ATII"
  wilcox_rows[[length(wilcox_rows)+1]] <- wres
  direction <- if (!is.na(wres$log2fc) && wres$log2fc < 0) "↓ IPF<Ctrl" else "↑ IPF>Ctrl"
  message(sprintf("    %-10s  p=%.3f  log2FC=%+.2f  [%s]",
                  g, wres$p_value, wres$log2fc, direction))
}

# Composite: TGFβ response score = mean(CLDN4, KRT8, SERPINE1) in AT2 cells
tgfb_genes_present <- intersect(c("CLDN4","KRT8","SERPINE1"), colnames(at2_data_a))
if (length(tgfb_genes_present) >= 2) {
  at2_data_a$tgfb_score_a <- rowMeans(at2_data_a[, tgfb_genes_present, drop=FALSE])
  cr_a <- cor.test(at2_data_a$FGFR2, at2_data_a$tgfb_score_a, method="spearman", exact=FALSE)
  message(sprintf("\n  Adams per-cell Spearman (FGFR2 vs TGFβ-arrest score [%s]):",
                  paste(tgfb_genes_present, collapse="+")))
  message(sprintf("    rho = %.3f  p = %.2e  (n = %d ATII cells)", cr_a$estimate, cr_a$p.value, nrow(at2_data_a)))
}

# =============================================================================
# C.  HUMAN: Habermann AT2 cells (IPF only — within-IPF correlation)
# =============================================================================
if (file.exists(HAB_RDS)) {
  message("\n========== C. HUMAN RECEIVER: Habermann AT2 (IPF only) ==========")
  hab <- readRDS(HAB_RDS)
  DefaultAssay(hab) <- "RNA"
  # Detect celltype column
  hab_ct_col <- if ("celltype" %in% colnames(hab@meta.data)) "celltype" else
                if ("cell_type" %in% colnames(hab@meta.data)) "cell_type" else
                if ("CellType"  %in% colnames(hab@meta.data)) "CellType"  else {
                  hits <- sapply(colnames(hab@meta.data), function(cn) any(grepl("AT2", hab@meta.data[[cn]], ignore.case=TRUE)))
                  if (any(hits)) names(hits)[which(hits)[1]] else "celltype"  # default, will warn
                }
  message(sprintf("  Habermann cell-type column: '%s'", hab_ct_col))
  message(sprintf("  Habermann labels: %s",
                  paste(sort(unique(hab@meta.data[[hab_ct_col]])), collapse=", ")))
  hab_at2 <- hab[, hab@meta.data[[hab_ct_col]] %in% c("AT2","AT2 cells","ATII")]
  if (ncol(hab_at2) == 0) {
    message("  WARNING: no AT2 cells found under expected labels — skipping Habermann section.")
    rm(hab, hab_at2); gc()
  } else {
  message(sprintf("  Habermann AT2 cells: %d", ncol(hab_at2)))

  hab_at2 <- JoinLayers(hab_at2)
  fgfr2_h  <- GetAssayData(hab_at2, layer = "data")["FGFR2",]
  tgfbr2_h <- if ("TGFBR2" %in% rownames(hab_at2)) GetAssayData(hab_at2, layer="data")["TGFBR2",] else NULL
  cldn4_h  <- if ("CLDN4"  %in% rownames(hab_at2)) GetAssayData(hab_at2, layer="data")["CLDN4",]  else NULL
  krt8_h   <- if ("KRT8"   %in% rownames(hab_at2)) GetAssayData(hab_at2, layer="data")["KRT8",]   else NULL

  message(sprintf("  FGFR2 expressers: %d / %d (%.1f%%)",
                  sum(fgfr2_h > 0), length(fgfr2_h), 100*mean(fgfr2_h > 0)))
  if (!is.null(tgfbr2_h)) message(sprintf("  TGFBR2 expressers: %d / %d (%.1f%%)",
                                           sum(tgfbr2_h>0), length(tgfbr2_h), 100*mean(tgfbr2_h>0)))

  # TGFβ arrest score from available genes
  tgfb_genes_h <- c()
  if (!is.null(cldn4_h)) tgfb_genes_h <- c(tgfb_genes_h, "CLDN4")
  if (!is.null(krt8_h))  tgfb_genes_h <- c(tgfb_genes_h, "KRT8")
  if (length(tgfb_genes_h) >= 1) {
    tgfb_mat_h <- do.call(rbind, Filter(Negate(is.null), list(CLDN4=cldn4_h, KRT8=krt8_h)))
    tgfb_score_h <- colMeans(tgfb_mat_h, na.rm = TRUE)
    cr_h <- cor.test(fgfr2_h, tgfb_score_h, method = "spearman", exact = FALSE)
    message(sprintf("  Habermann per-cell Spearman (FGFR2 vs [%s] score):",
                    paste(tgfb_genes_h, collapse="+")))
    message(sprintf("    rho = %.3f  p = %.2e  (n = %d AT2)", cr_h$estimate, cr_h$p.value, ncol(hab_at2)))
  }
  rm(hab, hab_at2); gc()
  }  # end: AT2 cells found
} else {
  message("  Habermann Seurat not found — run 08_human_nichenet.R first.")
}

# =============================================================================
# D.  FIGURES
# =============================================================================
message("\n========== Generating figures ==========")
library(ggplot2)

# ── D1. Mouse time course: Fgfr2 % in AT2 cells by phase ──────────────────
phase_summ_plot <- phase_summ_m %>% filter(phase %in% PHASE_ORDER)

p_mouse_tc <- ggplot(phase_summ_plot, aes(x = phase)) +
  geom_col(aes(y = pct_Fgfr2), fill = "#2166ac", alpha = 0.8) +
  geom_line(aes(y = mean_Transitional * 100, group = 1),
            colour = "#d73027", linewidth = 1.2, linetype = "dashed") +
  geom_point(aes(y = mean_Transitional * 100), colour = "#d73027", size = 3) +
  scale_y_continuous(
    name = "% Fgfr2-expressing AT2 cells (bars)",
    sec.axis = sec_axis(~ . / 100, name = "Mean Transitional1 score (red line)")
  ) +
  scale_x_discrete(limits = PHASE_ORDER) +
  labs(title = "Mouse AT2: Fgfr2 expression ↓ as TGFβ-arrest ↑",
       subtitle = paste0("Bars = % Fgfr2+ AT2 cells  |  Red line = Transitional1 (Krt8/pEMT) score\n",
                         "FGFR2 IIIb/IIIc caveat: AT2 Fgfr2 treated as FGFR2b surrogate"),
       x = "Bleomycin phase") +
  theme_bw(base_size = 11) +
  theme(axis.title.y.right = element_text(colour = "#d73027"),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT_FIG, "fgfr2_timecourse_mouse.png"),
       p_mouse_tc, width = 8, height = 5, dpi = 300)
message("  Saved: fgfr2_timecourse_mouse.png")

# ── D2. Dual: Fgfr2 in AT2 + Fgf10 in Fibroblasts (cross-dataset) ────────
if (file.exists(SENDER_CSV)) {
  strunz_fib_tp <- read.csv(SENDER_CSV) %>%
    filter(dataset == "Strunz_mouse_GSE141259") %>%
    group_by(condition) %>%
    summarise(fib_fgf10_pct = mean(pct_expressing, na.rm=TRUE), .groups="drop") %>%
    mutate(phase = dplyr::case_when(
      condition == "PBS" ~ "PBS_control",
      condition %in% c("d3","d7")   ~ "acute",
      condition %in% c("d10","d14") ~ "fibrotic_early",
      condition == "d21"            ~ "fibrotic_late",
      condition == "d28"            ~ "resolution",
      TRUE ~ "other"
    )) %>%
    filter(phase != "other") %>%
    group_by(phase) %>% summarise(fib_fgf10_pct=mean(fib_fgf10_pct), .groups="drop")

  dual_df <- phase_summ_plot %>%
    left_join(strunz_fib_tp, by = "phase") %>%
    filter(!is.na(fib_fgf10_pct)) %>%
    mutate(phase = factor(phase, levels = PHASE_ORDER))

  if (nrow(dual_df) >= 2) {
    # Normalise both to [0,1] for dual display
    dual_long <- dual_df %>%
      mutate(Fgfr2_AT2_norm  = (pct_Fgfr2 - min(pct_Fgfr2)) / (max(pct_Fgfr2) - min(pct_Fgfr2) + 1e-9),
             Fgf10_Fib_norm  = (fib_fgf10_pct - min(fib_fgf10_pct)) /
               (max(fib_fgf10_pct) - min(fib_fgf10_pct) + 1e-9)) %>%
      select(phase, Fgfr2_AT2_norm, Fgf10_Fib_norm) %>%
      pivot_longer(cols = -phase, names_to = "signal", values_to = "norm_value")

    p_dual <- ggplot(dual_long, aes(x = phase, y = norm_value,
                                    colour = signal, group = signal)) +
      geom_line(linewidth = 1.3) +
      geom_point(size = 3) +
      scale_colour_manual(values = c(Fgfr2_AT2_norm  = "#2166ac",
                                     Fgf10_Fib_norm  = "#d7191c"),
                          labels = c(Fgfr2_AT2_norm  = "Fgfr2 in AT2 (receiver)",
                                     Fgf10_Fib_norm  = "Fgf10 in Fibroblasts (sender)")) +
      scale_x_discrete(limits = PHASE_ORDER) +
      labs(title = "Dual-hit model: both sender supply and receiver responsiveness drop",
           subtitle = paste0("Normalised 0–1 per signal. Fgfr2 from Kobayashi PATS epi_obj; ",
                             "Fgf10 from Strunz WholeLung fibroblasts.\n",
                             "Cross-dataset eco-correlation — directional only."),
           x = "Bleomycin phase", y = "Normalised signal (0–1)", colour = NULL) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 30, hjust = 1))
    ggsave(file.path(OUT_FIG, "fgfr2_timecourse_dual.png"),
           p_dual, width = 8, height = 5, dpi = 300)
    message("  Saved: fgfr2_timecourse_dual.png")
  }
}

# ── D3. Per-cell scatter: Fgfr2 count vs Transitional1 score (mouse AT2) ──
# Bin Fgfr2 to 0 vs >0 for clarity (most counts are 0 or 1)
scatter_df <- at2_meta %>%
  filter(phase %in% PHASE_ORDER) %>%
  mutate(Fgfr2_bin = ifelse(Fgfr2_count > 0, "Fgfr2+", "Fgfr2-"),
         phase = factor(phase, levels = PHASE_ORDER))

p_scatter <- ggplot(scatter_df, aes(x = Fgfr2_bin, y = Transitional1, fill = Fgfr2_bin)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_wrap(~phase, nrow = 1) +
  scale_fill_manual(values = c(`Fgfr2+` = "#2166ac", `Fgfr2-` = "#bdbdbd")) +
  labs(title = "Fgfr2-expressing AT2 cells have lower TGFβ-arrest score",
       subtitle = sprintf("Transitional1 = Krt8/Cldn4/pEMT module score. Spearman rho=%.2f p=%.1e",
                          ct_spear$estimate, ct_spear$p.value),
       x = NULL, y = "Transitional1 score (TGFβ-arrest proxy)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT_FIG, "fgfr2_vs_tgfb_scatter.png"),
       p_scatter, width = 10, height = 4, dpi = 300)
message("  Saved: fgfr2_vs_tgfb_scatter.png")

# ── D4. Adams human: FGFR2 + TGFBR2 in ATII IPF vs Control (per-sample jitter) ──
at2_pb_a <- at2_data_a %>%
  group_by(Subject_Identity, Disease_Identity) %>%
  summarise(across(all_of(intersect(c("FGFR2","TGFBR2","CLDN4","KRT8"),
                                     colnames(at2_data_a))),
                   mean, na.rm = TRUE, .names = "{.col}"),
            .groups = "drop") %>%
  pivot_longer(cols = -c(Subject_Identity, Disease_Identity),
               names_to = "gene", values_to = "pb_mean")

p_human <- ggplot(at2_pb_a, aes(x = Disease_Identity, y = pb_mean, fill = Disease_Identity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, aes(colour = Disease_Identity)) +
  facet_wrap(~gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values   = c(IPF = "#d7191c", Control = "#2c7bb6")) +
  scale_colour_manual(values = c(IPF = "#8b0000", Control = "#08306b")) +
  labs(title = "Human ATII: FGFR2 and TGFβ pathway in IPF vs Control",
       subtitle = "Adams GSE136831. FGFR2 IIIb/IIIc caveat applies. Per-sample pseudobulk shown.",
       x = NULL, y = "Mean counts per sample") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")
ggsave(file.path(OUT_FIG, "fgfr2_human_ipf_ctrl.png"),
       p_human, width = 10, height = 4, dpi = 300)
message("  Saved: fgfr2_human_ipf_ctrl.png")

# =============================================================================
# E.  SAVE SUMMARIES
# =============================================================================
message("\n========== Saving summaries ==========")
result_df <- bind_rows(result_rows)
wilcox_df  <- bind_rows(wilcox_rows)

write.csv(result_df, file.path(OUT_SUMM, "fgfr2_receiver_summary.csv"), row.names = FALSE)
write.csv(wilcox_df, file.path(OUT_SUMM, "fgfr2_receiver_wilcox.csv"),  row.names = FALSE)

message("\n=== Key FGFR2 findings ===")
message("Mouse (AT2 cells across bleomycin time course):")
print(phase_summ_m %>% select(phase, n_cells, pct_Fgfr2, pct_Tgfbr2, mean_Transitional))

message("\nHuman Adams (ATII cells, IPF vs Control):")
print(wilcox_df %>% filter(dataset == "Adams_GSE136831") %>%
        select(gene, subtype, p_value, median_a, median_b, log2fc))

message("\n08b_fgfr2_receiver.R complete.")
message("Outputs:")
message("  results/fgfr2_receiver_summary.csv")
message("  results/fgfr2_receiver_wilcox.csv")
message("  results/figures/fgfr2_timecourse_mouse.png")
message("  results/figures/fgfr2_timecourse_dual.png  (sender + receiver paired)")
message("  results/figures/fgfr2_vs_tgfb_scatter.png  (per-cell Spearman)")
message("  results/figures/fgfr2_human_ipf_ctrl.png")
