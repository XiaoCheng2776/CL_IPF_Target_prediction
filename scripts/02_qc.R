#!/usr/bin/env Rscript
# =============================================================================
# 02_qc.R — QC and preprocessing for all mouse datasets
# Run from repo root:  Rscript scripts/02_qc.R
# Output: results/qc/*.rds  (one filtered Seurat object per dataset)
#         results/figures/qc_*.png
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)       # v5
  library(harmony)
  library(Matrix)
  library(tidyverse)
  library(patchwork)
})

ROOT    <- here::here()
DATA    <- file.path(ROOT, "data/mouse")
OUT_RDS <- file.path(ROOT, "results/qc")
OUT_FIG <- file.path(ROOT, "results/figures")
dir.create(OUT_RDS, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

MITO_PATTERN  <- "^mt-"
NFEATURE_MIN  <- 200
NFEATURE_MAX  <- 7000
PERCENT_MT_MAX <- 20

# Load an MTX and orient to genes × cells regardless of deposit orientation.
# readMM returns dgTMatrix; Seurat 5 requires dgCMatrix — coerce before returning.
read_mtx_oriented <- function(mtx_path, barcodes_path, features_path) {
  barcodes <- read_lines(barcodes_path)
  features <- read_lines(features_path)
  m        <- readMM(gzfile(mtx_path))
  if (nrow(m) == length(barcodes) && ncol(m) == length(features)) {
    message("  [orient] cells×genes (", nrow(m), "×", ncol(m), ") — transposing")
    m <- t(m)
  } else if (nrow(m) == length(features) && ncol(m) == length(barcodes)) {
    message("  [orient] genes×cells (", nrow(m), "×", ncol(m), ") — ok")
  } else {
    stop("Dimension mismatch: matrix ", nrow(m), "×", ncol(m),
         "  barcodes=", length(barcodes), "  genes=", length(features))
  }
  rownames(m) <- features
  colnames(m) <- barcodes
  as(m, "CsparseMatrix")
}

# Generic QC + processing function -----------------------------------------
# normalize_method: "sctransform" (default) or "lognorm"
# shard_by: column name to split on before SCTransform (avoids dense-residuals OOM
#           on large objects); shards are merged and variable features selected
#           via SelectIntegrationFeatures before RunPCA.
process_one <- function(counts, dataset_name, meta = NULL,
                        batch_key = NULL, n_dims = 30,
                        normalize_method = c("sctransform", "lognorm"),
                        shard_by = NULL) {
  normalize_method <- match.arg(normalize_method)

  obj <- CreateSeuratObject(counts = counts, project = dataset_name,
                            min.cells = 3, min.features = 200)
  if (!is.null(meta)) {
    shared <- intersect(colnames(obj), rownames(meta))
    obj <- AddMetaData(obj, metadata = meta[shared, , drop = FALSE])
  }
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = MITO_PATTERN)

  p_pre <- VlnPlot(obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                   ncol = 3, pt.size = 0) +
    plot_annotation(title = paste(dataset_name, "— pre-filter"))
  ggsave(file.path(OUT_FIG, paste0("qc_prefilter_", dataset_name, ".png")),
         p_pre, width = 12, height = 4, dpi = 300)

  obj <- subset(obj,
                subset = nFeature_RNA > NFEATURE_MIN &
                         nFeature_RNA < NFEATURE_MAX &
                         percent.mt   < PERCENT_MT_MAX)
  message(dataset_name, ": ", ncol(obj), " cells after QC filter")

  if (normalize_method == "lognorm") {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
  } else if (!is.null(shard_by) && shard_by %in% colnames(obj@meta.data)) {
    message("  SCTransform: sharding by '", shard_by, "'")
    shards <- SplitObject(obj, split.by = shard_by)
    shards <- lapply(shards, function(s) {
      SCTransform(s, method = "glmGamPoi", vst.flavor = "v2",
                  ncells = min(5000L, ncol(s)), verbose = FALSE)
    })
    obj <- merge(shards[[1]], y = shards[-1])
    VariableFeatures(obj) <- SelectIntegrationFeatures(
      object.list = shards, nfeatures = 3000
    )
  } else {
    obj <- SCTransform(obj, method = "glmGamPoi", vst.flavor = "v2",
                       ncells = min(5000L, ncol(obj)), verbose = FALSE)
  }

  obj <- RunPCA(obj, verbose = FALSE)
  if (!is.null(batch_key) && batch_key %in% colnames(obj@meta.data)) {
    obj <- RunHarmony(obj, group.by.vars = batch_key, verbose = FALSE)
  }
  red <- if (!is.null(batch_key) && batch_key %in% colnames(obj@meta.data)) "harmony" else "pca"
  obj <- RunUMAP(obj, reduction = red, dims = 1:n_dims, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = red, dims = 1:n_dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  obj
}

# =============================================================================
# 1.  Strunz GSE141259  (primary — bleomycin time course, whole lung)
# =============================================================================
strunz_dir <- file.path(DATA, "strunz_gse141259")

if (file.exists(file.path(OUT_RDS, "strunz_wholelun_qc.rds")) &&
    file.exists(file.path(OUT_RDS, "strunz_highres_qc.rds"))) {
  message("=== Strunz GSE141259 — skipping (RDS exists) ===")
} else {
  message("=== Strunz GSE141259 ===")

  # WholeLung resolution (sender characterisation)
  wl_counts <- read_mtx_oriented(
    mtx_path      = file.path(strunz_dir, "GSE141259_WholeLung_rawcounts.mtx.gz"),
    barcodes_path = file.path(strunz_dir, "GSE141259_WholeLung_barcodes.txt.gz"),
    features_path = file.path(strunz_dir, "GSE141259_WholeLung_genes.txt.gz")
  )
  wl_meta <- read.csv(file.path(strunz_dir, "GSE141259_WholeLung_cellinfo.csv.gz"),
                      row.names = 1)
  strunz_wl <- process_one(wl_counts, "strunz_wholelun", meta = wl_meta,
                           batch_key = "identifier", normalize_method = "lognorm")

  # HighResolution subset (receiver trajectory — prefer for AT2/transitional/AT1)
  hr_counts <- read_mtx_oriented(
    mtx_path      = file.path(strunz_dir, "GSE141259_HighResolution_rawcounts.mtx.gz"),
    barcodes_path = file.path(strunz_dir, "GSE141259_HighResolution_barcodes.txt.gz"),
    features_path = file.path(strunz_dir, "GSE141259_HighResolution_genes.txt.gz")
  )
  hr_meta <- read_tsv(file.path(strunz_dir, "GSE141259_HighResolution_cellinfo.csv.gz"),
                      show_col_types = FALSE) %>%
    column_to_rownames("cell_barcode")
  strunz_hr <- process_one(hr_counts, "strunz_highres", meta = hr_meta,
                           batch_key = "sample_id", shard_by = "sample_id")

  saveRDS(strunz_wl, file.path(OUT_RDS, "strunz_wholelun_qc.rds"))
  saveRDS(strunz_hr, file.path(OUT_RDS, "strunz_highres_qc.rds"))
}

# =============================================================================
# 2.  Kobayashi GSE141634  (PATS sorted epithelium — MTECplus fraction)
# =============================================================================
if (file.exists(file.path(OUT_RDS, "kobayashi_qc.rds"))) {
  message("=== Kobayashi GSE141634 — skipping (RDS exists) ===")
} else {
  message("=== Kobayashi GSE141634 ===")
  kob_dir <- file.path(DATA, "kobayashi_gse141634")
  kob_mat <- read.table(
    gzfile(file.path(kob_dir, "GSM4210295_MTECplus.tsv.gz")),
    header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
  )
  kob_counts <- Matrix(as.matrix(kob_mat), sparse = TRUE)
  kobayashi <- process_one(kob_counts, "kobayashi_pats")
  saveRDS(kobayashi, file.path(OUT_RDS, "kobayashi_qc.rds"))
}

# =============================================================================
# 3.  Choi GSE145031  (DATPs — 6 samples: PBS/Day14/Day28 × Tomato/nonTomato)
# =============================================================================
if (file.exists(file.path(OUT_RDS, "choi_qc.rds"))) {
  message("=== Choi GSE145031 — skipping (RDS exists) ===")
} else {
  message("=== Choi GSE145031 ===")
  choi_dir <- file.path(DATA, "choi_gse145031")

  choi_samples <- list(
    PBS_Tomato    = "GSM4304609_PBS_AT2_Tomato",
    PBS_nonTomato = "GSM4304610_PBS_AT2_nonTomato",
    D14_Tomato    = "GSM4304611_Day14_AT2_Tomato",
    D14_nonTomato = "GSM4304612_Day14_AT2_nonTomato",
    D28_Tomato    = "GSM4304613_Day28_AT2_Tomato",
    D28_nonTomato = "GSM4304614_Day28_AT2_nonTomato"
  )

  choi_list <- lapply(names(choi_samples), function(nm) {
    pfx <- choi_samples[[nm]]
    counts <- ReadMtx(
      mtx      = file.path(choi_dir, paste0(pfx, ".mtx.gz")),
      cells    = file.path(choi_dir, paste0(pfx, "_barcodes.tsv.gz")),
      features = file.path(choi_dir, paste0(pfx, "_gene.tsv.gz")),
      feature.column = 2
    )
    obj <- CreateSeuratObject(counts = counts, project = nm,
                              min.cells = 3, min.features = 200)
    obj$sample <- nm
    obj$timepoint <- sub("_(Tomato|nonTomato)$", "", nm)
    obj$lineage   <- ifelse(grepl("Tomato$", nm) & !grepl("non", nm),
                            "Tomato_pos", "Tomato_neg")
    obj
  })

  choi_merged <- merge(choi_list[[1]], y = choi_list[-1],
                       add.cell.ids = names(choi_samples))
  choi_merged <- JoinLayers(choi_merged)   # v5: collapse per-sample count layers before extraction
  choi <- process_one(GetAssayData(choi_merged, layer = "counts"),
                      "choi_datp",
                      meta             = choi_merged@meta.data,
                      batch_key        = "sample",
                      normalize_method = "lognorm")
  saveRDS(choi, file.path(OUT_RDS, "choi_qc.rds"))
}

message("02_qc.R complete. Objects saved to results/qc/")
