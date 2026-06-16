#!/usr/bin/env Rscript
# =============================================================================
# 08_fgf10_sender.R — FGF10 sender-side deficiency in IPF
# Run from repo root:  Rscript scripts/08_fgf10_sender.R
#
# Hypothesis: FGF10 from fibroblasts is reduced in IPF, supporting FGF10
# supplementation to promote AT2→AT1 repair. Tests the SUPPLY side (this script)
# vs the receiver side (FGFR2b on AT2, confirmed in 04_nichenet.R).
#
# Three datasets:
#   A. Adams   GSE136831  — human, IPF + Control, n=32 IPF / 28 Control subjects
#   B. Habermann GSE135893 — human, IPF + Control, n=12 IPF / 10 Control subjects
#   C. Strunz  GSE141259  — mouse, bleomycin time course (whole-lung)
#
# Three FGF10 metrics (required because FGF10 is low-abundance and dropout-prone;
# per-cell mean alone is unreliable for sparse genes):
#   (a) pct_expressing  = % fibroblasts with FGF10 > 0
#   (b) pseudobulk      = mean FGF10 per subject → Wilcoxon IPF vs Control
#   (c) mean_expressers = mean(FGF10 | FGF10 > 0) — per-cell, expressers only
#
# Compositional analysis (mechanistically distinct from per-cell downregulation):
#   % PLIN2+ lipofibroblasts in IPF vs Control
#   Fewer lipofibroblasts (fewer FGF10-producing cells) vs lower expression
#   per cell differ in therapeutic implication.
#
# Caveats honoured:
#   - FGF10 is low-abundance; all three metrics required for reliable inference
#   - Whole-lung dissociation (Strunz) distorts fibroblast composition; note it
#   - Habermann PLIN2+ annotation context: PLIN2+ may mark a pathological
#     fibroblast subset in human IPF rather than classical alveolar lipofibroblasts
#   - Association only — does not demonstrate rescue; supports supplementation rationale
#   - Adams has no pre-annotated PLIN2+ label; lipofibroblasts inferred from PLIN2 expr
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(tidyverse)
  library(data.table)
})

ROOT      <- here::here()
OUT_FIG   <- file.path(ROOT, "results/figures")
OUT_SUMM  <- file.path(ROOT, "results")
ADAMS_DIR <- file.path(ROOT, "data/human/adams_gse136831")
HAB_DIR   <- file.path(ROOT, "data/human/habermann_gse135893")
MOUSE_RDS <- file.path(ROOT, "results/qc/strunz_wholelun_qc.rds")

dir.create(OUT_FIG,  showWarnings = FALSE)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Three FGF10 metrics for a named count vector + sample IDs
fgf10_metrics_df <- function(counts, sample_ids, subtype, condition, dataset) {
  pb <- tapply(counts, sample_ids, mean, na.rm = TRUE)  # pseudobulk per sample
  data.frame(
    dataset             = dataset,
    condition           = condition,
    subtype             = subtype,
    n_cells             = length(counts),
    n_samples           = length(unique(sample_ids)),
    pct_expressing      = 100 * mean(counts > 0),
    mean_expressers     = if (any(counts > 0)) mean(counts[counts > 0]) else NA_real_,
    pseudobulk_mean     = mean(pb, na.rm = TRUE),
    pseudobulk_median   = median(pb, na.rm = TRUE),
    stringsAsFactors    = FALSE
  )
}

# Wilcoxon on pseudobulk values; returns p-value and direction summary
pb_wilcox <- function(pb_ipf, pb_ctrl) {
  if (length(pb_ipf) < 2 || length(pb_ctrl) < 2) {
    return(data.frame(p_value = NA_real_, median_ipf = median(pb_ipf),
                      median_ctrl = median(pb_ctrl), log2fc = NA_real_))
  }
  wt <- wilcox.test(pb_ipf, pb_ctrl, exact = FALSE)
  data.frame(
    p_value      = wt$p.value,
    median_ipf   = median(pb_ipf,  na.rm = TRUE),
    median_ctrl  = median(pb_ctrl, na.rm = TRUE),
    log2fc       = log2((median(pb_ipf, na.rm=TRUE) + 0.001) /
                        (median(pb_ctrl, na.rm=TRUE) + 0.001))
  )
}

# Extract multiple genes from a gzipped MTX in a single awk pass (memory-efficient)
# Returns a cells × genes data.frame of raw counts (0-filled for non-expressers)
# gene_indices: named integer vector (gene_name -> 1-based row index in MTX)
extract_genes_mtx <- function(mtx_gz, gene_indices, barcodes, skip_lines = 3) {
  if (length(gene_indices) == 0) return(NULL)
  idx_str <- paste(gene_indices, collapse = "|")
  # Single awk pass: match any of the gene indices in field 1
  cmd <- sprintf(
    "zcat '%s' | awk 'NR>%d && ($1 ~ /^(%s)$/) {print $1, $2, $3}'",
    mtx_gz, skip_lines, idx_str
  )
  message(sprintf("  Extracting %d gene(s) from %s (single awk pass — may take ~5 min)…",
                  length(gene_indices), basename(mtx_gz)))
  dt <- tryCatch(
    fread(cmd = cmd, header = FALSE, col.names = c("gene_idx","cell_idx","count")),
    error = function(e) {
      message("  awk extraction failed: ", e$message)
      data.table(gene_idx = integer(), cell_idx = integer(), count = numeric())
    }
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

# Determine how many header lines to skip in a gzipped MTX
mtx_skip_lines <- function(mtx_gz) {
  con <- gzcon(file(mtx_gz, "rb"))
  on.exit(close(con))
  lines <- readLines(con, n = 10)
  sum(startsWith(lines, "%")) + 1L  # comment lines + dimension line
}

# MTX dimensions from header
mtx_dims <- function(mtx_gz) {
  con <- gzcon(file(mtx_gz, "rb"))
  on.exit(close(con))
  lines <- readLines(con, n = 10)
  dim_line <- lines[!startsWith(lines, "%")][1]
  as.integer(strsplit(trimws(dim_line), "\\s+")[[1]])
}

# Collect results
summary_rows <- list()
wilcox_rows  <- list()

# =============================================================================
# A.  ADAMS GSE136831 (human IPF + Control)
# =============================================================================
message("\n========== A. ADAMS GSE136831 ==========")

ADAMS_MTX  <- file.path(ADAMS_DIR, "GSE136831_RawCounts_Sparse.mtx.gz")
ADAMS_BAR  <- file.path(ADAMS_DIR, "GSE136831_AllCells.cellBarcodes.txt.gz")
ADAMS_GENE <- file.path(ADAMS_DIR, "GSE136831_AllCells.GeneIDs.txt.gz")
ADAMS_META <- file.path(ADAMS_DIR, "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz")

# A1. Load metadata
meta_a <- fread(ADAMS_META, sep = "\t", header = TRUE)
setnames(meta_a, "CellBarcode_Identity", "barcode")
message(sprintf("  Adams metadata: %d cells, %d columns", nrow(meta_a), ncol(meta_a)))
message("  Disease breakdown: "); print(table(meta_a$Disease_Identity))

# Keep only IPF and Control (exclude COPD)
meta_a_filt <- meta_a[Disease_Identity %in% c("IPF","Control")]
message(sprintf("  IPF+Control cells: %d  (IPF=%d, Control=%d)",
                nrow(meta_a_filt),
                sum(meta_a_filt$Disease_Identity=="IPF"),
                sum(meta_a_filt$Disease_Identity=="Control")))

# Stromal cells only
stromal_a <- meta_a_filt[CellType_Category == "Stromal"]
# Focus on Fibroblast and Myofibroblast (exclude pericyte/SMC for FGF10 analysis)
fb_a <- stromal_a[Manuscript_Identity %in% c("Fibroblast","Myofibroblast")]
message(sprintf("  Fibroblast+Myofibroblast cells: %d", nrow(fb_a)))
print(table(fb_a$Manuscript_Identity, fb_a$Disease_Identity))

# A2. Load barcodes and genes to build index mapping
message("  Loading Adams barcodes and gene list…")
barcodes_a <- fread(ADAMS_BAR, header = FALSE)[[1]]

# Adams gene file: two-column TSV with header ("Ensembl_GeneID", "HGNC_EnsemblAlt_GeneID")
genes_a_tbl <- fread(ADAMS_GENE, header = TRUE)
message(sprintf("  Adams gene file: %d rows × %d cols; cols: %s",
                nrow(genes_a_tbl), ncol(genes_a_tbl),
                paste(colnames(genes_a_tbl), collapse=", ")))

# Use gene symbols (col2) as primary matching; fall back to Ensembl (col1) if needed
genes_a_sym <- if (ncol(genes_a_tbl) >= 2) make.unique(as.character(genes_a_tbl[[2]])) else character(0)
genes_a_ens <- as.character(genes_a_tbl[[1]])
message(sprintf("  Adams: %d barcodes, %d genes (symbol sample: %s)",
                length(barcodes_a), nrow(genes_a_tbl),
                paste(head(genes_a_sym, 3), collapse=", ")))

# Determine target genes
target_genes_a <- c("FGF10","PLIN2","ACTA2","HAS1","COL1A1")
ENSEMBL_MAP_A  <- c(FGF10 ="ENSG00000107480", PLIN2 ="ENSG00000028116",
                    ACTA2 ="ENSG00000107796", HAS1  ="ENSG00000105509",
                    COL1A1="ENSG00000108821")

gene_idx_a <- setNames(match(target_genes_a, genes_a_sym), target_genes_a)
missing_a  <- names(gene_idx_a)[is.na(gene_idx_a)]
if (length(missing_a) > 0) {
  message("  Symbols not found; trying Ensembl IDs for: ", paste(missing_a, collapse=", "))
  for (gn in missing_a) {
    ens_id <- ENSEMBL_MAP_A[gn]
    idx    <- match(ens_id, genes_a_ens)
    if (!is.na(idx)) gene_idx_a[gn] <- idx
  }
}
gene_idx_a <- gene_idx_a[!is.na(gene_idx_a)]
message(sprintf("  Genes matched for extraction: %s", paste(names(gene_idx_a), collapse=", ")))

# A3. MTX header check
skip_a <- mtx_skip_lines(ADAMS_MTX)
dims_a  <- mtx_dims(ADAMS_MTX)
message(sprintf("  Adams MTX dims: %d × %d (nnz=%d); skip_lines=%d",
                dims_a[1], dims_a[2], dims_a[3], skip_a))

# Determine orientation: rows=genes or rows=cells?
n_genes_expected <- nrow(genes_a_tbl)
n_cells_expected <- length(barcodes_a)
if (dims_a[1] == n_genes_expected && dims_a[2] == n_cells_expected) {
  message("  MTX orientation: rows=genes, cols=cells (standard)")
  mtx_row_is_gene <- TRUE
} else if (dims_a[1] == n_cells_expected && dims_a[2] == n_genes_expected) {
  message("  MTX orientation: rows=cells, cols=genes (transposed)")
  mtx_row_is_gene <- FALSE
  # For transposed MTX, swap: gene is in field 2, cell in field 1
  # Rebuild awk command accordingly below
} else {
  message(sprintf("  WARNING: MTX dims (%d × %d) don't match expected (%d genes × %d cells).",
                  dims_a[1], dims_a[2], n_genes_expected, n_cells_expected))
  message("  Assuming standard orientation (rows=genes).")
  mtx_row_is_gene <- TRUE
}

if (!mtx_row_is_gene) {
  # Transposed: field 2 = gene, field 1 = cell → swap in awk
  idx_str <- paste(gene_idx_a, collapse = "|")
  cmd_a <- sprintf(
    "zcat '%s' | awk 'NR>%d && ($2 ~ /^(%s)$/) {print $2, $1, $3}'",
    ADAMS_MTX, skip_a, idx_str
  )
  message("  Extracting genes (transposed orientation)…")
  dt_a <- tryCatch(
    fread(cmd = cmd_a, header = FALSE, col.names = c("gene_idx","cell_idx","count")),
    error = function(e) data.table(gene_idx=integer(), cell_idx=integer(), count=numeric())
  )
} else {
  dt_a <- NULL  # will use extract_genes_mtx below
}

if (mtx_row_is_gene) {
  expr_a <- extract_genes_mtx(ADAMS_MTX, gene_idx_a, barcodes_a, skip_lines = skip_a)
} else {
  # Reconstruct from dt_a
  rev_idx_a <- setNames(names(gene_idx_a), gene_idx_a)
  expr_a <- as.data.frame(
    matrix(0.0, nrow = length(barcodes_a), ncol = length(gene_idx_a),
           dimnames = list(barcodes_a, names(gene_idx_a)))
  )
  # Vectorized fill
  if (nrow(dt_a) > 0) {
    gnames_v <- rev_idx_a[as.character(dt_a$gene_idx)]
    valid_v  <- !is.na(gnames_v) & dt_a$cell_idx >= 1 & dt_a$cell_idx <= nrow(expr_a)
    for (gn in unique(na.omit(gnames_v[valid_v]))) {
      sel <- valid_v & gnames_v == gn
      expr_a[dt_a$cell_idx[sel], gn] <- dt_a$count[sel]
    }
  }
}

# A4. Merge expression with metadata for fibroblast cells
fb_barcodes_a <- intersect(fb_a$barcode, rownames(expr_a))
if (length(fb_barcodes_a) < nrow(fb_a)) {
  message(sprintf("  NOTE: %d fibroblast barcodes not found in MTX (format mismatch?)",
                  nrow(fb_a) - length(fb_barcodes_a)))
}
expr_fb_a <- expr_a[fb_barcodes_a, , drop = FALSE]
fb_data_a <- cbind(
  as.data.frame(fb_a[, .(barcode, Manuscript_Identity, Disease_Identity, Subject_Identity)]),
  expr_fb_a
)
message(sprintf("  Fibroblast cells with expression data: %d", nrow(fb_data_a)))
message(sprintf("  FGF10 expressers: %d / %d (%.1f%%)",
                sum(fb_data_a$FGF10 > 0, na.rm=TRUE), nrow(fb_data_a),
                100*mean(fb_data_a$FGF10 > 0, na.rm=TRUE)))

# A5. Subtype: within "Fibroblast", identify PLIN2+ (lipofibroblast) using PLIN2 expr
if ("PLIN2" %in% colnames(fb_data_a)) {
  fb_data_a$subtype <- dplyr::case_when(
    fb_data_a$Manuscript_Identity == "Myofibroblast" ~ "Myofibroblast",
    fb_data_a$Manuscript_Identity == "Fibroblast" & fb_data_a$PLIN2 > 0 ~ "PLIN2+ Fibroblast (lipofibroblast)",
    fb_data_a$Manuscript_Identity == "Fibroblast" ~ "Fibroblast (PLIN2-)"
  )
  message("\n  Adams fibroblast subtype breakdown (PLIN2-based):")
  print(table(fb_data_a$subtype, fb_data_a$Disease_Identity))
  message("  NOTE: PLIN2+ label inferred from PLIN2 expression, not pre-annotated.")
  message("  Whole-lung dissociation may distort PLIN2+ fraction — treat compositional change cautiously.")
} else {
  fb_data_a$subtype <- fb_data_a$Manuscript_Identity
  message("  PLIN2 not extracted — using Manuscript_Identity labels only.")
}

# A6. Three FGF10 metrics per subtype × condition
for (st in unique(fb_data_a$subtype)) {
  for (cond in c("IPF","Control")) {
    sel <- fb_data_a$subtype == st & fb_data_a$Disease_Identity == cond
    if (sum(sel) < 5) next
    row <- fgf10_metrics_df(
      counts     = fb_data_a$FGF10[sel],
      sample_ids = fb_data_a$Subject_Identity[sel],
      subtype    = st, condition = cond, dataset = "Adams_GSE136831"
    )
    summary_rows[[length(summary_rows)+1]] <- row
  }
}

# A7. Pseudobulk Wilcoxon: IPF vs Control within each subtype
message("\n  Adams pseudobulk Wilcoxon (IPF vs Control):")
for (st in unique(fb_data_a$subtype)) {
  ipf_pb  <- tapply(fb_data_a$FGF10[fb_data_a$subtype==st & fb_data_a$Disease_Identity=="IPF"],
                    fb_data_a$Subject_Identity[fb_data_a$subtype==st & fb_data_a$Disease_Identity=="IPF"],
                    mean, na.rm=TRUE)
  ctrl_pb <- tapply(fb_data_a$FGF10[fb_data_a$subtype==st & fb_data_a$Disease_Identity=="Control"],
                    fb_data_a$Subject_Identity[fb_data_a$subtype==st & fb_data_a$Disease_Identity=="Control"],
                    mean, na.rm=TRUE)
  wres <- pb_wilcox(ipf_pb, ctrl_pb)
  wres$dataset <- "Adams_GSE136831"; wres$subtype <- st
  wilcox_rows[[length(wilcox_rows)+1]] <- wres
  message(sprintf("    %-40s  p=%.3f  log2FC=%.2f  (IPF med=%.4f, Ctrl med=%.4f)",
                  st, wres$p_value, wres$log2fc, wres$median_ipf, wres$median_ctrl))
}

# A8. Compositional: PLIN2+ fraction per subject
if ("PLIN2" %in% colnames(fb_data_a)) {
  comp_a <- fb_data_a %>%
    filter(Manuscript_Identity == "Fibroblast") %>%
    group_by(Subject_Identity, Disease_Identity) %>%
    summarise(plin2_pct = 100 * mean(PLIN2 > 0), n_cells = n(), .groups = "drop")
  wcomp_a <- pb_wilcox(
    comp_a$plin2_pct[comp_a$Disease_Identity=="IPF"],
    comp_a$plin2_pct[comp_a$Disease_Identity=="Control"]
  )
  message(sprintf("\n  Adams compositional (PLIN2+ %% of Fibroblasts):"))
  message(sprintf("    IPF median=%.1f%%  Control median=%.1f%%  p=%.3f",
                  wcomp_a$median_ipf, wcomp_a$median_ctrl, wcomp_a$p_value))
  message("  NOTE: decrease in PLIN2+ fraction could reflect either fewer lipofibroblasts")
  message("  or dropout artifact from whole-lung dissociation (Strunz data shares this caveat).")
}

# =============================================================================
# B.  HABERMANN GSE135893 (human IPF + Control)
# =============================================================================
message("\n========== B. HABERMANN GSE135893 ==========")

HAB_META <- file.path(HAB_DIR, "GSE135893_IPF_metadata.csv.gz")
HAB_MTX  <- file.path(HAB_DIR, "GSE135893_matrix.mtx.gz")
HAB_GENE <- file.path(HAB_DIR, "GSE135893_genes.tsv.gz")
HAB_BAR  <- file.path(HAB_DIR, "GSE135893_barcodes.tsv.gz")

# B1. Load metadata
meta_h <- read.csv(HAB_META, row.names = 1)
message(sprintf("  Habermann metadata: %d cells", nrow(meta_h)))
message("  Diagnosis: "); print(table(meta_h$Diagnosis))

# Fibroblast subtypes — keep IPF and Control only
HAB_FB_TYPES <- c("Fibroblasts","PLIN2+ Fibroblasts","HAS1 High Fibroblasts","Myofibroblasts")
fb_meta_h <- meta_h %>%
  filter(celltype %in% HAB_FB_TYPES,
         Diagnosis %in% c("IPF","Control")) %>%
  mutate(condition = Diagnosis)

message(sprintf("  Fibroblast cells (IPF+Control): %d", nrow(fb_meta_h)))
print(table(fb_meta_h$celltype, fb_meta_h$condition))

# B2. Read gene list and barcode list
# Habermann gene file may have Ensembl IDs in col1 and symbols in col2 (10x format)
# or only gene symbols in col1 — detect and use symbols for matching
genes_h_raw <- read.table(gzfile(HAB_GENE), header=FALSE, stringsAsFactors=FALSE)
if (ncol(genes_h_raw) >= 2 && grepl("^ENSG", genes_h_raw[[1]][1])) {
  genes_h <- genes_h_raw[[2]]   # col2 = gene symbol
  message("  Habermann gene file: col1=EnsemblID, col2=symbol — using col2")
} else {
  genes_h <- genes_h_raw[[1]]   # col1 already gene symbols
}
bars_h <- read.table(gzfile(HAB_BAR), header=FALSE)[[1]]
message(sprintf("  Habermann: %d barcodes, %d genes (first 3: %s)",
                length(bars_h), length(genes_h), paste(head(genes_h,3), collapse=", ")))

# B3. Find FGF10 index; for Habermann we load only fibroblast barcodes (manageable size)
fb_barcodes_h <- rownames(fb_meta_h)
fb_barcodes_h <- intersect(fb_barcodes_h, bars_h)
fb_col_idx_h  <- match(fb_barcodes_h, bars_h)  # column indices for fibroblasts

fgf10_idx_h  <- match("FGF10",  genes_h)
plin2_idx_h  <- match("PLIN2",  genes_h)
acta2_idx_h  <- match("ACTA2",  genes_h)
message(sprintf("  FGF10 gene index: %s | PLIN2: %s | ACTA2: %s",
                fgf10_idx_h, plin2_idx_h, acta2_idx_h))

# B4. Load full MTX and subset to fibroblast columns (1.1 GB — manageable)
message("  Reading Habermann MTX (subset to fibroblast cells)…")
mat_h <- readMM(HAB_MTX)
# mat_h is genes × cells (standard orientation)
target_rows_h <- c(FGF10=fgf10_idx_h, PLIN2=plin2_idx_h, ACTA2=acta2_idx_h)
target_rows_h <- target_rows_h[!is.na(target_rows_h)]

expr_fb_h <- as.data.frame(
  t(as.matrix(mat_h[target_rows_h, fb_col_idx_h]))
)
colnames(expr_fb_h) <- names(target_rows_h)
rownames(expr_fb_h) <- fb_barcodes_h
rm(mat_h); gc()

fb_data_h <- cbind(
  fb_meta_h[fb_barcodes_h, c("celltype","condition","Sample_Name")],
  expr_fb_h
)
message(sprintf("  Habermann fibroblast cells with expression: %d", nrow(fb_data_h)))
message(sprintf("  FGF10 expressers: %d / %d (%.1f%%)",
                sum(fb_data_h$FGF10 > 0, na.rm=TRUE), nrow(fb_data_h),
                100*mean(fb_data_h$FGF10 > 0, na.rm=TRUE)))

# B5. Three metrics per subtype × condition
for (ct in HAB_FB_TYPES) {
  for (cond in c("IPF","Control")) {
    sel <- fb_data_h$celltype == ct & fb_data_h$condition == cond
    if (sum(sel) < 5) next
    row <- fgf10_metrics_df(
      counts     = fb_data_h$FGF10[sel],
      sample_ids = fb_data_h$Sample_Name[sel],
      subtype    = ct, condition = cond, dataset = "Habermann_GSE135893"
    )
    summary_rows[[length(summary_rows)+1]] <- row
  }
}

# B6. Pseudobulk Wilcoxon per subtype
message("\n  Habermann pseudobulk Wilcoxon (IPF vs Control):")
for (ct in HAB_FB_TYPES) {
  ipf_pb  <- tapply(fb_data_h$FGF10[fb_data_h$celltype==ct & fb_data_h$condition=="IPF"],
                    fb_data_h$Sample_Name[fb_data_h$celltype==ct & fb_data_h$condition=="IPF"],
                    mean, na.rm=TRUE)
  ctrl_pb <- tapply(fb_data_h$FGF10[fb_data_h$celltype==ct & fb_data_h$condition=="Control"],
                    fb_data_h$Sample_Name[fb_data_h$celltype==ct & fb_data_h$condition=="Control"],
                    mean, na.rm=TRUE)
  if (length(ipf_pb) < 1 || length(ctrl_pb) < 1) next
  wres <- pb_wilcox(ipf_pb, ctrl_pb)
  wres$dataset <- "Habermann_GSE135893"; wres$subtype <- ct
  wilcox_rows[[length(wilcox_rows)+1]] <- wres
  message(sprintf("    %-30s  p=%.3f  log2FC=%.2f  (IPF med=%.4f, Ctrl med=%.4f)",
                  ct, wres$p_value, wres$log2fc, wres$median_ipf, wres$median_ctrl))
}

# B7. Compositional: PLIN2+ fraction in IPF vs Control (pre-annotated in Habermann)
# PLIN2+ Fibroblasts: n_IPF=1167 vs n_Control=6 — a striking enrichment in IPF
# (see caveat: PLIN2+ in Habermann may be a pathological subtype, not classical
#  alveolar lipofibroblast — use expression data to double-check PLIN2 levels)
all_fb_meta_h <- meta_h %>%
  filter(celltype %in% HAB_FB_TYPES,
         Diagnosis %in% c("IPF","Control"))
comp_h <- all_fb_meta_h %>%
  group_by(Sample_Name, Diagnosis) %>%
  summarise(
    total_fb   = n(),
    plin2_plus = sum(celltype == "PLIN2+ Fibroblasts"),
    plin2_pct  = 100 * plin2_plus / total_fb,
    .groups = "drop"
  )
wcomp_h <- pb_wilcox(
  comp_h$plin2_pct[comp_h$Diagnosis=="IPF"],
  comp_h$plin2_pct[comp_h$Diagnosis=="Control"]
)
message(sprintf("\n  Habermann compositional (PLIN2+ %% of all fibroblasts):"))
message(sprintf("    IPF median=%.1f%%  Control median=%.1f%%  p=%.3f",
                wcomp_h$median_ipf, wcomp_h$median_ctrl, wcomp_h$p_value))
message("  NOTE: PLIN2+ Fibroblasts are ENRICHED in IPF in Habermann (n=1167 IPF vs 6 Control).")
message("  This is the OPPOSITE of the mouse alveolar lipofibroblast depletion hypothesis.")
message("  Human PLIN2+ may represent a distinct pathological fibroblast subset,")
message("  not the FGF10-secreting alveolar lipofibroblast. Validate with FGF10 expression levels.")

# =============================================================================
# C.  STRUNZ MOUSE — BLEOMYCIN TIME COURSE
# =============================================================================
# Strunz metadata uses:
#   cell.type  — cell type labels (Fibroblasts, Myofibroblasts, ...)
#   grouping   — time point collapsed: d3, d7, d10, d14, d21, d28, PBS
#   identifier — per-mouse sample ID: {mouse_id}_{Bleo|PBS}_{day}
#                → use as pseudobulk unit (3-5 mice per time point)
# =============================================================================
message("\n========== C. STRUNZ MOUSE (bleomycin time course) ==========")

library(Seurat)
mouse_obj <- readRDS(MOUSE_RDS)
message(sprintf("  Loaded: %d cells, %d genes", ncol(mouse_obj), nrow(mouse_obj)))

# Subset to fibroblasts
Idents(mouse_obj) <- "cell.type"
all_fb_types_m <- grep("[Ff]ibroblast", unique(mouse_obj$cell.type), value=TRUE)
message(sprintf("  Fibroblast types: %s", paste(all_fb_types_m, collapse=", ")))

# Score lipofibroblast markers before subsetting
lipo_markers_m <- intersect(c("Plin2","Plin3","Cd36","Pparg"), rownames(mouse_obj))
if (length(lipo_markers_m) >= 1) {
  mouse_obj <- AddModuleScore(mouse_obj, features=list(lipo_markers_m), name="LipoFb")
  message(sprintf("  Lipofibroblast module score (Plin2 axis): %s",
                  paste(lipo_markers_m, collapse=", ")))
}

fb_obj_m <- subset(mouse_obj, idents = all_fb_types_m)
rm(mouse_obj); gc()
fb_obj_m <- JoinLayers(fb_obj_m)
DefaultAssay(fb_obj_m) <- "RNA"

message(sprintf("  Mouse fibroblasts: %d cells", ncol(fb_obj_m)))
cat("  grouping distribution:\n"); print(sort(table(fb_obj_m$grouping)))

# Fgf10 counts
fgf10_m_counts <- GetAssayData(fb_obj_m, layer="counts")["Fgf10",]
message(sprintf("  Fgf10 expressers: %d / %d (%.2f%%)",
                sum(fgf10_m_counts > 0), length(fgf10_m_counts),
                100*mean(fgf10_m_counts > 0)))

fb_meta_m <- fb_obj_m@meta.data
fb_meta_m$fgf10_count <- as.numeric(fgf10_m_counts[rownames(fb_meta_m)])

# Parse condition from identifier (contains "_Bleo_" or "_PBS_")
fb_meta_m$condition_m <- ifelse(grepl("_PBS_", fb_meta_m$identifier), "PBS", "Bleo")

# Lipofibroblast subtype via PLIN2-axis module score
if ("LipoFb1" %in% colnames(fb_meta_m)) {
  lipo_thresh_m <- quantile(
    fb_meta_m$LipoFb1[fb_meta_m$cell.type == "Fibroblasts"], 0.75, na.rm=TRUE)
  fb_meta_m$subtype_m <- dplyr::case_when(
    fb_meta_m$cell.type == "Myofibroblasts" ~ "Myofibroblasts",
    fb_meta_m$cell.type == "Fibroblasts" & fb_meta_m$LipoFb1 >= lipo_thresh_m ~ "Fibroblasts (Plin2-high)",
    fb_meta_m$cell.type == "Fibroblasts" ~ "Fibroblasts (Plin2-low)",
    TRUE ~ fb_meta_m$cell.type
  )
} else {
  fb_meta_m$subtype_m <- fb_meta_m$cell.type
}

# Bleomycin phase grouping based on actual `grouping` column values
# (d3, d7, d10, d14, d21, d28, PBS)
# Bleomycin model: acute d1-7; early fibrotic d10-14; peak fibrosis d21; resolution d28
fb_meta_m$phase <- dplyr::case_when(
  fb_meta_m$grouping == "PBS"  ~ "PBS_control",
  fb_meta_m$grouping == "d3"   ~ "acute",
  fb_meta_m$grouping == "d7"   ~ "acute",
  fb_meta_m$grouping == "d10"  ~ "fibrotic_early",
  fb_meta_m$grouping == "d14"  ~ "fibrotic_early",
  fb_meta_m$grouping == "d21"  ~ "fibrotic_late",
  fb_meta_m$grouping == "d28"  ~ "resolution",
  TRUE                          ~ "other"
)

message("\n  Cell counts by subtype × grouping:")
print(table(fb_meta_m$subtype_m, fb_meta_m$grouping))

# Summary by phase
phase_summary_m <- fb_meta_m %>%
  group_by(phase, subtype_m) %>%
  summarise(
    n_cells         = n(),
    n_mice          = n_distinct(identifier),
    pct_fgf10_pos   = 100 * mean(fgf10_count > 0),
    mean_expressers = if (any(fgf10_count>0)) mean(fgf10_count[fgf10_count>0]) else NA_real_,
    .groups = "drop"
  )
message("\n  Fgf10+ fraction by phase:")
print(phase_summary_m, n=30)

# Metrics per grouping timepoint (using identifier as per-mouse unit)
for (tp in unique(fb_meta_m$grouping)) {
  for (st in unique(fb_meta_m$subtype_m)) {
    sel <- fb_meta_m$grouping == tp & fb_meta_m$subtype_m == st
    if (sum(sel) < 3) next
    row <- fgf10_metrics_df(
      counts     = fb_meta_m$fgf10_count[sel],
      sample_ids = fb_meta_m$identifier[sel],  # per-mouse pseudobulk unit
      subtype    = st,
      condition  = tp,
      dataset    = "Strunz_mouse_GSE141259"
    )
    summary_rows[[length(summary_rows)+1]] <- row
  }
}

# Pseudobulk Wilcoxon: bleomycin fibrotic phase vs PBS (all fibroblasts)
# identifier allows per-mouse pseudobulk (3-5 replicates per phase)
for (ph in c("fibrotic_early","fibrotic_late","resolution")) {
  bleo_pb <- tapply(
    fb_meta_m$fgf10_count[fb_meta_m$phase == ph],
    fb_meta_m$identifier[fb_meta_m$phase == ph],
    mean, na.rm=TRUE)
  pbs_pb <- tapply(
    fb_meta_m$fgf10_count[fb_meta_m$phase == "PBS_control"],
    fb_meta_m$identifier[fb_meta_m$phase == "PBS_control"],
    mean, na.rm=TRUE)
  if (length(bleo_pb) >= 2 && length(pbs_pb) >= 2) {
    wres <- pb_wilcox(bleo_pb, pbs_pb)
    wres$dataset <- "Strunz_mouse_GSE141259"; wres$subtype <- paste0("AllFibroblasts_", ph)
    wilcox_rows[[length(wilcox_rows)+1]] <- wres
    message(sprintf("  Mouse pseudobulk Wilcoxon (Bleo_%s vs PBS):  p=%.3f  log2FC=%.2f",
                    ph, wres$p_value, wres$log2fc))
  }
}

message("  NOTE: Strunz `identifier` encodes per-mouse samples — pseudobulk is valid.")
message("  Whole-lung dissociation may still distort fibroblast subtype composition (note it).")

# =============================================================================
# D.  FIGURES
# =============================================================================
message("\n========== Generating figures ==========")

library(ggplot2)

# --- D1: Human FGF10 pct expressing by subtype × condition ---
sum_df <- bind_rows(summary_rows) %>%
  filter(dataset %in% c("Adams_GSE136831","Habermann_GSE135893"))

if (nrow(sum_df) > 0) {
  p_pct <- ggplot(sum_df,
                  aes(x = reorder(subtype, pct_expressing), y = pct_expressing,
                      fill = condition)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    facet_wrap(~dataset, scales = "free_x") +
    scale_fill_manual(values = c(IPF="#d7191c", Control="#2c7bb6")) +
    coord_flip() +
    labs(title = "% FGF10-expressing fibroblasts (IPF vs Control)",
         subtitle = "Metric (a): fraction with FGF10 > 0; low-abundance dropouts mitigated",
         x = NULL, y = "% cells expressing FGF10", fill = "Condition") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT_FIG, "fgf10_pct_expressing.png"),
         p_pct, width = 11, height = 5, dpi = 300)
  message("  Saved: fgf10_pct_expressing.png")

  # --- D2: Pseudobulk mean bar with sample-level jitter ---
  # For Habermann (has enough samples to show spread)
  hab_fb_pb <- fb_data_h %>%
    group_by(celltype, condition, Sample_Name) %>%
    summarise(pb_fgf10 = mean(FGF10, na.rm=TRUE), .groups="drop")
  if (nrow(hab_fb_pb) > 0) {
    p_pb <- ggplot(hab_fb_pb, aes(x=celltype, y=pb_fgf10, fill=condition)) +
      geom_bar(data = hab_fb_pb %>% group_by(celltype, condition) %>%
                 summarise(pb_fgf10=mean(pb_fgf10), .groups="drop"),
               stat="identity", position=position_dodge(0.8), width=0.7, alpha=0.7) +
      geom_jitter(aes(colour=condition),
                  position=position_jitterdodge(dodge.width=0.8, jitter.width=0.15),
                  size=2, alpha=0.9) +
      scale_fill_manual(values = c(IPF="#d7191c", Control="#2c7bb6")) +
      scale_colour_manual(values = c(IPF="#8b0000", Control="#08306b")) +
      coord_flip() +
      labs(title = "Pseudobulk FGF10 per sample — Habermann",
           subtitle = "Metric (b): mean FGF10 per subject; dots = individual samples",
           x = NULL, y = "Mean FGF10 counts (per-sample)", fill = "Condition") +
      theme_bw(base_size = 11) +
      guides(colour = "none")
    ggsave(file.path(OUT_FIG, "fgf10_pseudobulk_habermann.png"),
           p_pb, width = 9, height = 5, dpi = 300)
    message("  Saved: fgf10_pseudobulk_habermann.png")
  }

  # --- D3: Mean-among-expressers bar ---
  p_expr <- ggplot(sum_df %>% filter(!is.na(mean_expressers)),
                   aes(x = reorder(subtype, mean_expressers), y = mean_expressers,
                       fill = condition)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    facet_wrap(~dataset, scales = "free_x") +
    scale_fill_manual(values = c(IPF="#d7191c", Control="#2c7bb6")) +
    coord_flip() +
    labs(title = "FGF10 per-cell expression (expressers only)",
         subtitle = "Metric (c): mean(FGF10 | FGF10 > 0); lower = per-cell downregulation",
         x = NULL, y = "Mean FGF10 counts among expressers", fill = "Condition") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT_FIG, "fgf10_mean_expressers.png"),
         p_expr, width = 11, height = 5, dpi = 300)
  message("  Saved: fgf10_mean_expressers.png")
}

# --- D4: Habermann compositional: PLIN2+ fraction per sample ---
if (nrow(comp_h) > 0) {
  p_comp <- ggplot(comp_h, aes(x=Diagnosis, y=plin2_pct, fill=Diagnosis)) +
    geom_boxplot(outlier.shape=NA, alpha=0.7) +
    geom_jitter(width=0.15, size=2.5, aes(colour=Diagnosis)) +
    scale_fill_manual(values = c(IPF="#d7191c", Control="#2c7bb6")) +
    scale_colour_manual(values = c(IPF="#8b0000", Control="#08306b")) +
    labs(title = "PLIN2+ fibroblast fraction — Habermann (compositional)",
         subtitle = sprintf("IPF med=%.1f%%, Control med=%.1f%%, p=%.3f (Wilcoxon pseudobulk)\nNOTE: PLIN2+ enriched in IPF here — may be pathological, not alveolar lipofibroblast",
                            wcomp_h$median_ipf, wcomp_h$median_ctrl, wcomp_h$p_value),
         x = "Diagnosis", y = "% PLIN2+ Fibroblasts of all fibroblasts") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(OUT_FIG, "fgf10_compositional_plin2_habermann.png"),
         p_comp, width = 5, height = 5, dpi = 300)
  message("  Saved: fgf10_compositional_plin2_habermann.png")
}

# --- D5: Mouse time course (using grouping column; per-mouse pseudobulk available) ---
GROUPING_ORDER <- c("PBS","d3","d7","d10","d14","d21","d28")
ts_data_m <- fb_meta_m %>%
  filter(grouping %in% GROUPING_ORDER) %>%
  group_by(grouping, condition_m, subtype_m) %>%
  summarise(pct_fgf10 = 100*mean(fgf10_count>0),
            n_cells   = n(),
            n_mice    = n_distinct(identifier),
            .groups   = "drop") %>%
  mutate(grouping = factor(grouping, levels=GROUPING_ORDER))

p_tc <- ggplot(ts_data_m, aes(x=grouping, y=pct_fgf10,
                               colour=subtype_m, group=interaction(subtype_m, condition_m),
                               linetype=condition_m)) +
  geom_line(linewidth=1.1) +
  geom_point(aes(size=n_mice), alpha=0.9) +
  scale_size_continuous(range=c(2,6), name="n mice") +
  scale_linetype_manual(values=c(Bleo="solid", PBS="dashed")) +
  labs(title = "Mouse Fgf10+ fibroblast fraction — Strunz bleomycin time course",
       subtitle = "Solid = bleomycin; dashed = PBS. Per-mouse pseudobulk Wilcoxon in output CSV.",
       x = "Time point", y = "% Fgf10-expressing fibroblasts",
       colour = "Subtype", linetype = "Condition") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle=30, hjust=1))
ggsave(file.path(OUT_FIG, "fgf10_mouse_timecourse.png"),
       p_tc, width = 9, height = 5, dpi = 300)
message("  Saved: fgf10_mouse_timecourse.png")

# =============================================================================
# E.  SUMMARY CSV
# =============================================================================
message("\n========== Saving summary ==========")

summary_df <- bind_rows(summary_rows)
wilcox_df  <- bind_rows(wilcox_rows) %>%
  select(dataset, subtype, p_value, median_ipf, median_ctrl, log2fc)

write.csv(summary_df, file.path(OUT_SUMM, "fgf10_sender_summary.csv"), row.names=FALSE)
write.csv(wilcox_df,  file.path(OUT_SUMM, "fgf10_sender_wilcox.csv"),  row.names=FALSE)

message("\n=== FGF10 sender summary (all datasets) ===")
print(head(summary_df %>% select(dataset, subtype, condition, n_cells,
                                  pct_expressing, mean_expressers, pseudobulk_mean), 40))

message("\n=== Wilcoxon pseudobulk results (IPF vs Control) ===")
print(head(wilcox_df, 20))

message("\n08_fgf10_sender.R complete.")
message("Outputs:")
message("  results/fgf10_sender_summary.csv")
message("  results/fgf10_sender_wilcox.csv")
message("  results/figures/fgf10_pct_expressing.png")
message("  results/figures/fgf10_pseudobulk_habermann.png")
message("  results/figures/fgf10_mean_expressers.png")
message("  results/figures/fgf10_compositional_plin2_habermann.png")
message("  results/figures/fgf10_mouse_timecourse.png")
