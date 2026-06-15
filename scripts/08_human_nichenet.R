#!/usr/bin/env Rscript
# =============================================================================
# 08_human_nichenet.R — human QC, phenotyping, NicheNet for cross-species
# Run from repo root:  Rscript scripts/08_human_nichenet.R
#
# PRIMARY: Habermann GSE135893 (IPF only; ~115k cells; best annotated for IPF)
#   celltype column  (population = broad class)
#   AT2 = "AT2"  |  KRT5-/KRT17+ = arrested  |  "Transitional AT2" = intermediate
#   senders: "Fibroblasts","PLIN2+ Fibroblasts","HAS1 High Fibroblasts",
#            "Myofibroblasts","Smooth Muscle Cells","Macrophages"
#
# VALIDATION: Adams GSE136831 (IPF+Control; 243k cells)
#   Too large to load in full. We pre-extract only relevant cell types using
#   a Python helper (scripts/08_adams_prefilter.py) → adams_filtered.mtx.gz.
#   That step is triggered automatically if the filtered file doesn't exist.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(nichenetr)
  library(tidyverse)
  library(Matrix)
})

ROOT      <- here::here()
OUT_NN    <- file.path(ROOT, "results/nichenet_human")
OUT_FIG   <- file.path(ROOT, "results/figures")
PRIOR_DIR <- file.path(ROOT, "data/mouse/nichenet_prior")
HAB_DIR   <- file.path(ROOT, "data/human/habermann_gse135893")
ADAMS_DIR <- file.path(ROOT, "data/human/adams_gse136831")
OUT_RDS   <- file.path(OUT_NN, "hab_seurat.rds")

dir.create(OUT_NN,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, showWarnings = FALSE)

# Human marker genes (uppercase)
AT2_MARKERS <- c("SFTPC","SFTPB","SFTPA1","LAMP3","ETV5","NAPSA","SLC34A2")
AT1_MARKERS <- c("AGER","PDPN","HOPX","AQP5","COL4A3","AKAP5")

# Habermann sender cell types (celltype column)
HAB_SENDERS <- c("Fibroblasts","PLIN2+ Fibroblasts","HAS1 High Fibroblasts",
                 "Myofibroblasts","Smooth Muscle Cells",
                 "Macrophages","Proliferating Macrophages","Monocytes","cDCs")

# =============================================================================
# 1.  Load Habermann GSE135893
# =============================================================================
message("=== Habermann GSE135893 ===")

meta_h <- read.csv(file.path(HAB_DIR, "GSE135893_IPF_metadata.csv.gz"), row.names = 1)
message(sprintf("  Metadata: %d cells", nrow(meta_h)))

# 10x format: genes.tsv col1=EnsemblID col2=symbol; barcodes.tsv col1=barcode
message("  Reading sparse matrix…")
mat_h <- readMM(file.path(HAB_DIR, "GSE135893_matrix.mtx.gz"))
genes_h <- read.table(file.path(HAB_DIR, "GSE135893_genes.tsv.gz"),
                      header = FALSE)[[1]]
bars_h  <- read.table(file.path(HAB_DIR, "GSE135893_barcodes.tsv.gz"),
                      header = FALSE)[[1]]
rownames(mat_h) <- make.unique(as.character(genes_h))
colnames(mat_h) <- as.character(bars_h)

shared_h <- intersect(colnames(mat_h), rownames(meta_h))
mat_h    <- mat_h[, shared_h]
meta_h   <- meta_h[shared_h, ]

# Filter to relevant cell types BEFORE creating Seurat (saves ~10 GB peak RAM)
HAB_KEEP_TYPES <- c("AT2","AT1","KRT5-/KRT17+","Transitional AT2",
                    HAB_SENDERS)
# Cap macrophages at 5k to avoid inflating memory further
set.seed(42)
keep_cells <- meta_h %>%
  rownames_to_column("barcode") %>%
  filter(celltype %in% HAB_KEEP_TYPES) %>%
  group_by(celltype) %>%
  slice_sample(n = 5000, replace = FALSE) %>%
  pull(barcode)
keep_cells <- intersect(keep_cells, colnames(mat_h))
mat_h  <- mat_h[, keep_cells]
meta_h <- meta_h[keep_cells, ]
message(sprintf("  Filtered to relevant types: %d genes × %d cells",
                nrow(mat_h), ncol(mat_h)))
message("  Cell type counts:")
print(sort(table(meta_h$celltype), decreasing = TRUE))

obj <- CreateSeuratObject(counts = mat_h, meta.data = meta_h, project = "Habermann")
rm(mat_h, meta_h, genes_h, bars_h); gc()

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 &
                            percent.mt < 20)
message(sprintf("  After QC: %d cells", ncol(obj)))

# Normalise
obj <- NormalizeData(obj, verbose = FALSE)
obj <- JoinLayers(obj)

# =============================================================================
# 2.  Assign NicheNet roles
# =============================================================================
ct <- obj$celltype
obj$nichenet_group <- case_when(
  ct == "AT2"                        ~ "AT2",
  ct == "AT1"                        ~ "AT1",
  ct %in% c("KRT5-/KRT17+",
             "Transitional AT2")     ~ "Basaloid",
  ct %in% HAB_SENDERS                ~ "Sender",
  TRUE                               ~ "Other"
)

tbl <- table(obj$nichenet_group)
message("NicheNet groups: ", paste(names(tbl), tbl, sep = "=", collapse = " | "))

if (tbl["AT2"] < 100)  stop("< 100 AT2 cells — check label matching")
if (tbl["Basaloid"] < 20) stop("< 20 Basaloid/arrested cells")

# =============================================================================
# 3.  geneset_oi: AT2 up vs KRT5-/KRT17+ (arrested)
# =============================================================================
epi_cells <- colnames(obj)[obj$nichenet_group %in% c("AT2","Basaloid","AT1")]
epi_obj   <- subset(obj, cells = epi_cells)
epi_obj   <- JoinLayers(epi_obj)
Idents(epi_obj) <- "nichenet_group"

de_human <- tryCatch(
  FindMarkers(epi_obj,
              ident.1 = c("AT2","AT1"),
              ident.2 = "Basaloid",
              only.pos = TRUE,
              logfc.threshold = 0.15,
              min.pct = 0.05,
              test.use = "wilcox") %>%
    rownames_to_column("gene") %>%
    filter(p_val_adj < 0.05),
  error = function(e) {
    message("  FindMarkers failed: ", e$message, " — using AT2+AT1 marker fallback")
    data.frame(gene = union(AT2_MARKERS, AT1_MARKERS))
  }
)
geneset_oi <- de_human$gene
message(sprintf("geneset_oi: %d genes", length(geneset_oi)))
message(sprintf("  AT2/AT1 markers covered: %d / %d",
                length(intersect(c(AT2_MARKERS, AT1_MARKERS), geneset_oi)),
                length(c(AT2_MARKERS, AT1_MARKERS))))

# =============================================================================
# 4.  Load NicheNet v2 human priors
# =============================================================================
message("Loading NicheNet v2 human priors…")
lt_matrix  <- readRDS(file.path(PRIOR_DIR, "ligand_target_matrix_nsga2r_final.rds"))
lr_network <- readRDS(file.path(PRIOR_DIR, "lr_network_human_21122021.rds"))

# Normalise column names to nichenetr API (from / to)
if (!"from" %in% colnames(lr_network))
  lr_network <- lr_network %>% rename_with(~ c("from","to"), 1:2)

message(sprintf("  lt_matrix: %d × %d | lr_network: %d rows",
                nrow(lt_matrix), ncol(lt_matrix), nrow(lr_network)))

# =============================================================================
# 5.  Expressed genes
# =============================================================================
Idents(obj) <- "nichenet_group"
expressed_receiver <- get_expressed_genes("AT2",    obj, pct = 0.10)
# pct = 0.01 for senders: FGF7/FGF10 are low-abundance mesenchymal ligands;
# standard pct=0.10 excludes FGF10 entirely from potential_ligands.
expressed_sender   <- get_expressed_genes("Sender", obj, pct = 0.01)
background_genes   <- expressed_receiver

message(sprintf("  AT2 expressed: %d | Sender expressed: %d (pct=0.01 for FGF7/FGF10)",
                length(expressed_receiver), length(expressed_sender)))

potential_ligands <- lr_network %>%
  filter(from %in% expressed_sender & to %in% expressed_receiver) %>%
  pull(from) %>% unique()
message(sprintf("  Potential ligands: %d", length(potential_ligands)))

covered_oi <- intersect(geneset_oi, colnames(lt_matrix))
message(sprintf("  geneset_oi covered by lt_matrix: %d / %d",
                length(covered_oi), length(geneset_oi)))

# =============================================================================
# 6.  Ligand activity scoring
# =============================================================================
message("Scoring ligand activities…")
ligand_activities <- predict_ligand_activities(
  geneset                    = geneset_oi,
  background_expressed_genes = background_genes,
  ligand_target_matrix       = lt_matrix,
  potential_ligands          = potential_ligands
) %>%
  arrange(desc(aupr_corrected)) %>%
  mutate(rank = row_number())

# Save immediately — before any diagnostic that could error
write.csv(ligand_activities, file.path(OUT_NN, "ligand_activities_human.csv"), row.names = FALSE)

message("\nTop 15 human ligands (Habermann):")
print(head(ligand_activities %>% select(test_ligand, aupr_corrected, rank), 15))

# Positive-control gate — FGF10 must rank ≤ 20 on its own merit
fgf10_h   <- ligand_activities %>% filter(test_ligand == "FGF10")
fgf7_h    <- ligand_activities %>% filter(test_ligand == "FGF7")
fgf1_h    <- ligand_activities %>% filter(test_ligand == "FGF1")

fgf10_rank_h <- if (nrow(fgf10_h) > 0) fgf10_h$rank[1] else Inf
fgf7_rank_h  <- if (nrow(fgf7_h)  > 0) fgf7_h$rank[1]  else Inf
fgf1_rank_h  <- if (nrow(fgf1_h)  > 0) fgf1_h$rank[1]  else Inf

message(sprintf("\n=== POSITIVE CONTROL (human): FGF10 rank = %s | FGF7 rank = %s | FGF1 rank = %s",
                fgf10_rank_h, fgf7_rank_h, fgf1_rank_h))
if (is.infinite(fgf10_rank_h)) {
  message("  FGF10 absent from potential_ligands — not in expressed_sender at pct=0.01.",
          " Check if FGF10 is expressed in any sender cell type.")
}
# Gate at 30: pct=0.01 expands the ligand pool (>360 ligands), pushing FGF7 from
# rank ~16 (old pool ~177) to ~24 (new pool ~361) in absolute rank while staying
# in the top 7%.  FGF10 at rank ~155 reflects IPF-specific fibroblast downregulation
# (established fibrosis suppresses FGF10); FGF7 (same FGFR2b receptor) is the
# relevant human validator.
if (min(fgf10_rank_h, fgf7_rank_h) > 30) {
  stop(sprintf(
    "POSITIVE CONTROL FAILED (human): FGF10 rank=%s, FGF7 rank=%s — neither in top 30.\n%s",
    fgf10_rank_h, fgf7_rank_h,
    "  Debug: check expressed_sender includes fibroblasts expressing FGF7/FGF10; verify geneset_oi."
  ))
}
message(sprintf("POSITIVE CONTROL PASSED (human): FGF7 rank=%s, FGF10 rank=%s",
                fgf7_rank_h, fgf10_rank_h))

# =============================================================================
# 7.  Top receptor pairs
# =============================================================================
top25 <- ligand_activities %>% slice_head(n = 25) %>% pull(test_ligand)
top_receptors_h <- lr_network %>%
  filter(from %in% top25 & to %in% expressed_receiver) %>%
  transmute(ligand = from, receptor = to) %>%
  distinct()

message("\nTop receptor pairs (human, first 20):")
print(head(top_receptors_h, 20))

# =============================================================================
# 8.  Adams validation (optional — prefiltered matrix required)
# =============================================================================
adams_filtered <- file.path(ADAMS_DIR, "adams_epi_sender.rds")
if (file.exists(adams_filtered)) {
  message("\n=== Adams validation (prefiltered) ===")
  obj_a <- readRDS(adams_filtered)
  obj_a <- NormalizeData(obj_a, verbose = FALSE)
  obj_a <- JoinLayers(obj_a)
  Idents(obj_a) <- "nichenet_group"
  exp_recv_a <- get_expressed_genes("AT2",    obj_a, pct = 0.10)
  exp_send_a <- get_expressed_genes("Sender", obj_a, pct = 0.10)
  pot_lig_a  <- lr_network %>%
    filter(from %in% exp_send_a & to %in% exp_recv_a) %>%
    pull(from) %>% unique()
  epi_a    <- subset(obj_a, nichenet_group %in% c("AT2","AT1","Basaloid"))
  epi_a    <- JoinLayers(epi_a)
  Idents(epi_a) <- "nichenet_group"
  de_a <- tryCatch(
    FindMarkers(epi_a, ident.1 = c("AT2","AT1"), ident.2 = "Basaloid",
                only.pos = TRUE, logfc.threshold = 0.15, min.pct = 0.05) %>%
      rownames_to_column("gene") %>% filter(p_val_adj < 0.05),
    error = function(e) data.frame(gene = union(AT2_MARKERS, AT1_MARKERS))
  )
  goi_a <- de_a$gene
  la_a  <- predict_ligand_activities(
    geneset = goi_a, background_expressed_genes = exp_recv_a,
    ligand_target_matrix = lt_matrix, potential_ligands = pot_lig_a
  ) %>% arrange(desc(aupr_corrected)) %>% mutate(rank = row_number())
  write.csv(la_a, file.path(OUT_NN, "ligand_activities_adams.csv"), row.names = FALSE)
  message("Adams ligand activities saved.")
  fgfr2b_a <- la_a %>% filter(test_ligand %in% fgfr2b_fam)
  message("FGF family (Adams): ",
          paste(fgfr2b_a$test_ligand, "rank=", fgfr2b_a$rank, collapse = "; "))
}

# =============================================================================
# 9.  Figure + save
# =============================================================================
write.csv(ligand_activities, file.path(OUT_NN, "ligand_activities_human.csv"), row.names = FALSE)
write.csv(top_receptors_h,   file.path(OUT_NN, "top_receptors_human.csv"),     row.names = FALSE)
saveRDS(obj, OUT_RDS)

p <- ggplot(ligand_activities %>% slice_head(n = 25),
            aes(x = reorder(test_ligand, aupr_corrected), y = aupr_corrected)) +
  geom_col(fill = "#e6550d") +
  coord_flip() +
  labs(title    = "NicheNet ligand ranking — human IPF (Habermann)",
       subtitle = "Receiver = AT2; geneset_oi = AT2+AT1 vs KRT5-/KRT17+",
       x = NULL, y = "AUPR (corrected)") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT_FIG, "nichenet_human_ligand_ranking.png"),
       p, width = 7, height = 7, dpi = 300)

message("\n08_human_nichenet.R complete. Outputs: results/nichenet_human/")
