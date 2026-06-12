#!/usr/bin/env Rscript
# =============================================================================
# 03_phenotype.R — epithelial identification, trajectory, and fate labelling
# Run from repo root:  Rscript scripts/03_phenotype.R
# Input:  results/qc/strunz_highres_qc.rds  (primary receiver dataset)
#         results/qc/kobayashi_qc.rds
#         results/qc/choi_qc.rds
# Output: results/markers/epi_obj.rds           (labelled epithelial subset)
#         results/markers/fate_labels.csv
#         results/figures/phenotype_*.png
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(slingshot)    # pseudotime
  library(tidyverse)
  library(patchwork)
})

ROOT    <- here::here()
OUT_MRK <- file.path(ROOT, "results/markers")
OUT_FIG <- file.path(ROOT, "results/figures")
dir.create(OUT_MRK, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# Marker gene sets (mouse, lowercase) ----------------------------------------
at2_markers        <- c("Sftpc","Sftpb","Sftpa1","Lamp3","Etv5","Napsa","Slc34a2")
at1_markers        <- c("Ager","Pdpn","Hopx","Aqp5","Col4a3","Akap5")
airway_club        <- c("Scgb1a1","Scgb3a2","Cyp2f2")
basal_markers      <- c("Krt5","Krt15","Trp63")
transitional_mouse <- c("Krt8","Cldn4","Krt18","Sfn","Cdkn1a","Lgals3","Tpm1")
senescence         <- c("Cdkn1a","Cdkn2a","Trp53","Glb1")
partial_emt        <- c("Vim","Fn1","Cdh2","Snai2","Zeb1")
regen_success      <- c(at2_markers, at1_markers)

# =============================================================================
# 1.  Load QC'd objects and merge for trajectory analysis
# =============================================================================
message("Loading QC'd objects…")
strunz_hr <- readRDS(file.path(ROOT, "results/qc/strunz_highres_qc.rds"))
kobayashi  <- readRDS(file.path(ROOT, "results/qc/kobayashi_qc.rds"))
choi       <- readRDS(file.path(ROOT, "results/qc/choi_qc.rds"))

# =============================================================================
# 2.  Epithelial identification
#     The Strunz HighResolution cellinfo already carries Strunz-defined labels in
#     the "cell_type" column: AT2, AT2 activated, AT1, Krt8+ ADI, Basal,
#     Ciliated, Club, etc.  Use those directly; add module scores for scoring only.
# =============================================================================
epi_types <- c("AT2", "AT2 activated", "AT1", "Krt8+ ADI")

if ("cell_type" %in% colnames(strunz_hr@meta.data)) {
  Idents(strunz_hr) <- "cell_type"
  epi_obj <- subset(strunz_hr, idents = epi_types)
  message("Epithelial cells from Strunz labels: ", ncol(epi_obj))
} else {
  # Fallback: infer from module scores if metadata join failed
  strunz_hr <- AddModuleScore(strunz_hr, features = list(at2_markers), name = "AT2_Score")
  strunz_hr <- AddModuleScore(strunz_hr, features = list(at1_markers), name = "AT1_Score")
  strunz_hr <- AddModuleScore(strunz_hr,
                               features = list(transitional_mouse), name = "Transitional")
  epi_obj <- subset(strunz_hr,
                    subset = AT2_Score1 > 0 | AT1_Score1 > 0 | Transitional1 > 0)
  message("Epithelial cells (module score fallback): ", ncol(epi_obj))
}

# Add module scores on the epithelial subset for fate labelling
epi_obj <- AddModuleScore(epi_obj, features = list(regen_success),      name = "RegenSuccess")
epi_obj <- AddModuleScore(epi_obj, features = list(transitional_mouse), name = "Transitional")
epi_obj <- AddModuleScore(epi_obj, features = list(senescence),         name = "Senescence")
epi_obj <- AddModuleScore(epi_obj, features = list(partial_emt),        name = "PartialEMT")

# =============================================================================
# 3.  Inject Strunz paper UMAP coordinates as a reference reduction
#     read_tsv confirms the file is tab-delimited; reuse umap_1/umap_2 from
#     the original embedding rather than relying solely on the re-computed one.
# =============================================================================
cellinfo <- read_tsv(
  file.path(ROOT, "data/mouse/strunz_gse141259/GSE141259_HighResolution_cellinfo.csv.gz"),
  col_select = c(cell_barcode, umap_1, umap_2, cell_type, time_point, louvain_cluster),
  show_col_types = FALSE
) %>% column_to_rownames("cell_barcode")

shared <- intersect(colnames(epi_obj), rownames(cellinfo))
if (length(shared) < ncol(epi_obj))
  warning(ncol(epi_obj) - length(shared), " epithelial cells absent from cellinfo; ",
          "Strunz UMAP will cover ", length(shared), "/", ncol(epi_obj), " cells only.")

strunz_coords <- as.matrix(cellinfo[shared, c("umap_1", "umap_2")])
colnames(strunz_coords) <- c("STUMAP_1", "STUMAP_2")
epi_obj[["strunz_umap"]] <- CreateDimReducObject(
  embeddings = strunz_coords[colnames(epi_obj)[colnames(epi_obj) %in% shared], ],
  key = "STUMAP_", assay = DefaultAssay(epi_obj)
)

# =============================================================================
# 4.  UMAP of epithelial subset
# =============================================================================
p_umap <- DimPlot(epi_obj, label = TRUE) +
  FeaturePlot(epi_obj, features = c("Sftpc","Ager","Krt8"), ncol = 3)
ggsave(file.path(OUT_FIG, "phenotype_epi_umap.png"),
       p_umap, width = 14, height = 5, dpi = 300)

# FGFR2b isoform flag: standard Fgfr2 quantification does NOT separate IIIb
# (epithelial) from IIIc (mesenchymal). Use with caution; annotate in outputs.
p_fgfr2 <- FeaturePlot(epi_obj, features = "Fgfr2") +
  labs(title = "Fgfr2 (NOTE: IIIb + IIIc isoforms combined)")
ggsave(file.path(OUT_FIG, "phenotype_fgfr2_flag.png"),
       p_fgfr2, width = 6, height = 5, dpi = 300)

# =============================================================================
# 5.  Pseudotime — slingshot (AT2 start, AT1 end)
# =============================================================================
# Slingshot requires a reduced-dim embedding and cluster labels.
# Adjust start.clus / end.clus after inspecting cluster IDs above.
rd  <- Embeddings(epi_obj, "umap")
cl  <- Idents(epi_obj)

sds <- slingshot::slingshot(rd, clusterLabels = cl,
                             start.clus = "AT2",
                             end.clus   = "AT1")
epi_obj$pseudotime <- slingshot::slingPseudotime(sds)[, 1]

# =============================================================================
# 6.  Fate labelling — transcriptional score-based; timepoint = sanity check only
#
#     healthy_baseline: d14_PBS cells — reference anchor, excluded from the
#       completing-vs-arrested contrast AND from Geneformer training.
#       Rationale: resting AT2 in an uninjured lung did not complete a
#       regenerative transit; merging them into "completing" would conflate
#       homeostatic AT2 identity with productive post-injury regeneration.
#
#     completing / arrested: defined by RegenSuccess vs Transitional module
#       scores on bleomycin cells only. Timepoint cross-tab is a post-hoc
#       sanity check — it must not define the labels.
# =============================================================================
bleo_cells  <- epi_obj$time_point != "d14_PBS"
q_regen_hi  <- quantile(epi_obj$RegenSuccess1[bleo_cells],  0.75, na.rm = TRUE)
q_regen_lo  <- quantile(epi_obj$RegenSuccess1[bleo_cells],  0.25, na.rm = TRUE)
q_trans_hi  <- quantile(epi_obj$Transitional1[bleo_cells],  0.75, na.rm = TRUE)
q_trans_lo  <- quantile(epi_obj$Transitional1[bleo_cells],  0.25, na.rm = TRUE)

epi_obj$fate_label <- dplyr::case_when(
  epi_obj$time_point == "d14_PBS"                                              ~ "healthy_baseline",
  epi_obj$cell_type  == "AT1"                                                  ~ "completing",
  epi_obj$RegenSuccess1 > q_regen_hi & epi_obj$Transitional1 < q_trans_lo     ~ "completing",
  epi_obj$Transitional1 > q_trans_hi & epi_obj$RegenSuccess1 < q_regen_lo     ~ "arrested",
  TRUE                                                                          ~ "intermediate"
)
message("\n=== Fate label distribution ===")
print(table(epi_obj$fate_label))

# Sanity check: arrested fraction should rise at late timepoints.
# Day 28 excluded (n=48; too sparse for reliable percentages).
message("\n=== Arrested % by timepoint (sanity check; day 28 excluded) ===")
ct <- epi_obj@meta.data %>%
  filter(!time_point %in% c("d14_PBS", "day 28")) %>%
  count(time_point, fate_label) %>%
  group_by(time_point) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  filter(fate_label == "arrested") %>%
  arrange(time_point)
print(ct, n = Inf)
# Expect pct to trend upward from day 2 → day 21/36/54.
# If it does not, debug score thresholds or epithelial subset composition.

p_fate <- DimPlot(epi_obj, group.by = "fate_label",
                  cols = c(completing       = "#2ca25f",
                           arrested         = "#de2d26",
                           intermediate     = "#bdbdbd",
                           healthy_baseline = "#fdae6b")) +
  labs(title = "Fate labels (score-based; d14_PBS = healthy_baseline)")
ggsave(file.path(OUT_FIG, "phenotype_fate_labels.png"),
       p_fate, width = 7, height = 5, dpi = 300)

# =============================================================================
# 7.  Save
#     healthy_baseline cells are retained in epi_obj for reference but are
#     flagged in fate_labels.csv; downstream scripts (04, 06) must filter them
#     out before running the completing-vs-arrested contrast.
# =============================================================================
write.csv(
  epi_obj@meta.data[, intersect(
    c("fate_label","pseudotime","cell_type","time_point","louvain_cluster",
      "RegenSuccess1","Transitional1","Senescence1","PartialEMT1"),
    colnames(epi_obj@meta.data)
  )],
  file.path(OUT_MRK, "fate_labels.csv")
)
saveRDS(epi_obj, file.path(OUT_MRK, "epi_obj.rds"))
message("03_phenotype.R complete.")
