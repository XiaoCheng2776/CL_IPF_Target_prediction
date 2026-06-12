#!/usr/bin/env python3
"""
08_adams_prefilter.py — extract only relevant cell types from Adams GSE136831
Produces data/human/adams_gse136831/adams_epi_sender.rds (via R) from a filtered MTX.
Run from repo root:  conda run -n gf python scripts/08_adams_prefilter.py

Adams MTX is genes × cells (rows=genes, cols=cells), standard GEO MTX format.
We stream-filter to keep only relevant barcodes, write a small filtered MTX,
then call an inline R snippet to convert to Seurat RDS.
"""
import gzip, re, subprocess, sys
from pathlib import Path

ROOT      = Path(__file__).resolve().parents[1]
ADAMS_DIR = ROOT / "data/human/adams_gse136831"
OUT_MTX   = ADAMS_DIR / "adams_filtered.mtx.gz"
OUT_BARS  = ADAMS_DIR / "adams_filtered_barcodes.txt"
OUT_RDS   = ADAMS_DIR / "adams_epi_sender.rds"

KEEP_TYPES = {
    "ATII", "ATI", "Aberrant_Basaloid",
    "Fibroblast", "Myofibroblast", "SMC",
    "Macrophage", "Macrophage_Alveolar",
}
KEEP_DIAG  = {"IPF", "Control"}

# --- 1. Read metadata and build keep-barcode set ---
meta_file = ADAMS_DIR / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"
bar_col, ct_col, diag_col = None, None, None
keep_barcodes = set()

print("Reading Adams metadata…")
with gzip.open(meta_file, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    bar_col  = header.index("CellBarcode_Identity")
    ct_col   = header.index("Manuscript_Identity")
    diag_col = header.index("Disease_Identity")
    for line in f:
        parts = line.rstrip("\n").split("\t")
        barcode  = parts[bar_col].strip('"')
        ct       = parts[ct_col].strip('"')
        diag     = parts[diag_col].strip('"')
        if ct in KEEP_TYPES and diag in KEEP_DIAG:
            keep_barcodes.add(barcode)

print(f"  Keeping {len(keep_barcodes)} barcodes from {len(KEEP_TYPES)} cell types")

# --- 2. Read all barcodes file to get column indices ---
bar_file = ADAMS_DIR / "GSE136831_AllCells.cellBarcodes.txt.gz"
barcodes = []
with gzip.open(bar_file, "rt") as f:
    for line in f:
        barcodes.append(line.strip())

keep_idx = {i + 1 for i, b in enumerate(barcodes) if b in keep_barcodes}  # 1-based
keep_list = sorted(keep_idx)
old2new = {old: new for new, old in enumerate(keep_list, 1)}
n_new_cols = len(keep_list)
new_barcodes = [barcodes[i - 1] for i in keep_list]
print(f"  Column indices to keep: {n_new_cols}")

# Save filtered barcodes
with open(OUT_BARS, "w") as f:
    f.write("\n".join(new_barcodes) + "\n")

# --- 3. Stream-filter the MTX ---
mtx_file = ADAMS_DIR / "GSE136831_RawCounts_Sparse.mtx.gz"
print("Streaming MTX file — this may take a few minutes…")

entries = []
n_rows = 0
with gzip.open(mtx_file, "rt") as f:
    for line in f:
        if line.startswith("%"):
            continue
        parts = line.split()
        if len(parts) == 3 and n_rows == 0:
            # First non-comment line is dimensions
            n_rows = int(parts[0])
            # n_cols = int(parts[1])  # ignored
            n_rows_orig = n_rows
            break

    for line in f:
        parts = line.split()
        row, col, val = int(parts[0]), int(parts[1]), parts[2]
        if col in keep_idx:
            entries.append((row, old2new[col], val))

print(f"  Retained {len(entries)} non-zero entries")

with gzip.open(OUT_MTX, "wt") as f:
    f.write("%%MatrixMarket matrix coordinate real general\n")
    f.write(f"{n_rows_orig} {n_new_cols} {len(entries)}\n")
    for r, c, v in entries:
        f.write(f"{r} {c} {v}\n")
print(f"  Written: {OUT_MTX}")

# --- 4. Call R to build Seurat RDS ---
rscript = subprocess.run(["which", "Rscript"], capture_output=True, text=True).stdout.strip()
if not rscript:
    rscript = "/home/shiqi/miniforge3/envs/ipf/bin/Rscript"

r_code = f"""
suppressPackageStartupMessages({{library(Seurat); library(Matrix); library(tidyverse)}})
root <- "{ROOT}"
adams_dir <- "{ADAMS_DIR}"

mat <- readMM("{OUT_MTX}")
bars  <- readLines("{OUT_BARS}")
genes <- read.table(
  gzfile(file.path(adams_dir, "GSE136831_AllCells.GeneIDs.txt.gz")),
  header = FALSE)[[1]]
rownames(mat) <- make.unique(as.character(genes))
colnames(mat) <- bars

meta <- read.table(
  gzfile(file.path(adams_dir, "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz")),
  header = TRUE, sep = "\\t")
rownames(meta) <- meta$CellBarcode_Identity
meta <- meta[bars, ]

obj <- CreateSeuratObject(counts = mat, meta.data = meta, project = "Adams")
obj$cell_type  <- obj$Manuscript_Identity
obj$diagnosis  <- obj$Disease_Identity
obj$dataset    <- "Adams"

obj$nichenet_group <- dplyr::case_when(
  obj$cell_type %in% c("ATII")              ~ "AT2",
  obj$cell_type %in% c("ATI")               ~ "AT1",
  obj$cell_type %in% c("Aberrant_Basaloid") ~ "Basaloid",
  obj$cell_type %in% c("Fibroblast","Myofibroblast","SMC",
                        "Macrophage","Macrophage_Alveolar") ~ "Sender",
  TRUE ~ "Other"
)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
message(sprintf("Adams filtered Seurat: %d cells", ncol(obj)))
print(table(obj$nichenet_group))
saveRDS(obj, "{OUT_RDS}")
message("Saved: {OUT_RDS}")
"""

result = subprocess.run([rscript, "--vanilla", "-e", r_code],
                        capture_output=False)
sys.exit(result.returncode)
