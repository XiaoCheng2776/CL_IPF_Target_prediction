#!/usr/bin/env Rscript
# =============================================================================
# 09_manuscript_figures.R — publication-ready figures for the IPF FGF10 study
#
# Styling mirrors MSC VS lsc.R (cardiac GO-plot script):
#   COL_FGF   = "#E07B39"  warm orange  — FGF/regeneration highlight
#   COL_REGEN = "#B83232"  brick red    — agonist candidates / positive control
#   COL_CONT  = "#4472B8"  steel blue   — contrast / control condition
#   theme_classic, black text, str_wrap GO labels, 300 dpi PNG + combined PDF
#
# Outputs (results/figures/):
#   Fig1_overview.png          — pipeline schematic + UMAP + NicheNet bar
#   Fig2_candidates.png        — agonist vs block dot plot + composite score
#   Fig3_sender.png            — FGF10 fibroblast time course + human IPF/Ctrl
#   Fig4_receiver.png          — AT2 FGFR2 expresser % + Krt8/pEMT co-expression
#   manuscript_figures.pdf     — all 4 combined (one figure per page)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

ROOT    <- here::here()
OUT_FIG <- file.path(ROOT, "results/figures")
dir.create(OUT_FIG, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# GLOBAL STYLE (verbatim palette from cardiac script)
# ─────────────────────────────────────────────────────────────────────────────
COL_FGF   <- "#E07B39"   # warm orange  — FGF/regeneration positive control
COL_REGEN <- "#B83232"   # brick red    — agonist candidates
COL_CONT  <- "#4472B8"   # steel blue   — contrast / PBS / Control condition
COL_BLOCK <- "#6A5B9F"   # muted violet — TGFβ block candidates
COL_UNCRT <- "#888888"   # mid grey     — uncertain direction
COL_NS    <- "#D3D3D3"   # light grey   — non-significant / neutral

AX_SZ     <- 14
TTL_SZ    <- 15
STRIP_SZ  <- 13
DPI       <- 300

base_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text        = element_text(size = AX_SZ,   color = "black"),
    axis.title       = element_text(size = AX_SZ,   color = "black"),
    plot.title       = element_text(size = TTL_SZ,  color = "black", face = "bold"),
    plot.subtitle    = element_text(size = AX_SZ - 2, color = "black"),
    legend.text      = element_text(size = AX_SZ - 2, color = "black"),
    legend.title     = element_text(size = AX_SZ - 2, color = "black", face = "bold"),
    strip.text       = element_text(size = STRIP_SZ, color = "black", face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA)
  )

fmt_p <- function(p) {
  if (is.na(p)) return("n.s.")
  if (p >= 0.05) return("n.s.")
  if (p < 0.001) return(sprintf("p = %.2e", p))
  return(sprintf("p = %.3f", p))
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT_FIG, paste0(name, ".png")), p, width = w, height = h, dpi = DPI)
  message(sprintf("  Saved: %s.png (%g × %g in)", name, w, h))
  invisible(p)
}

# ─────────────────────────────────────────────────────────────────────────────
# Load shared objects (once)
# ─────────────────────────────────────────────────────────────────────────────
message("Loading epi_obj.rds …")
epi <- readRDS(file.path(ROOT, "results/markers/epi_obj.rds"))
ct_col <- if ("cell_type" %in% colnames(epi@meta.data)) "cell_type" else "cell.type"
epi_meta <- epi@meta.data
umap_coords <- as.data.frame(Embeddings(epi, "umap"))
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_df <- cbind(umap_coords, epi_meta[rownames(umap_coords), ])
umap_df$cell_label <- umap_df[[ct_col]]

# Phase for receiver analyses
tp_col <- if ("time_point" %in% colnames(epi_meta)) "time_point" else "timepoint"
umap_df$day_num <- suppressWarnings(
  as.numeric(gsub(".*?([0-9]+).*", "\\1", umap_df[[tp_col]]))
)
umap_df$phase <- dplyr::case_when(
  grepl("PBS|NC", umap_df[[tp_col]], ignore.case = TRUE) ~ "PBS",
  !is.na(umap_df$day_num) & umap_df$day_num <= 10       ~ "Acute",
  !is.na(umap_df$day_num) & umap_df$day_num <= 14       ~ "Fibrotic\nearly",
  !is.na(umap_df$day_num) & umap_df$day_num <= 21       ~ "Fibrotic\nlate",
  !is.na(umap_df$day_num) & umap_df$day_num > 21        ~ "Resolution",
  TRUE ~ "other"
)
PHASE_ORDER <- c("PBS", "Acute", "Fibrotic\nearly", "Fibrotic\nlate", "Resolution")

# AT2 Fgfr2 counts from SCT
at2_cells <- rownames(epi_meta)[epi_meta[[ct_col]] %in% c("AT2", "AT2 activated")]
fgfr2_cnt <- tryCatch(
  as.numeric(GetAssayData(epi[, at2_cells], assay = "SCT", layer = "counts")["Fgfr2", ]),
  error = function(e) rep(NA_real_, length(at2_cells))
)
names(fgfr2_cnt) <- at2_cells

trans_col <- if ("Transitional1" %in% colnames(epi_meta)) "Transitional1" else
             if ("Transitional11" %in% colnames(epi_meta)) "Transitional11" else NULL

message("epi_obj loaded. Freeing Seurat object …")
rm(epi); gc()

# Load CSVs
nn_df       <- read.csv(file.path(ROOT, "results/nichenet/ligand_activities.csv"))
cand_df     <- read.csv(file.path(ROOT, "results/candidates_ranked.csv"))
sender_summ <- read.csv(file.path(ROOT, "results/fgf10_sender_summary.csv"))
sender_wilx <- read.csv(file.path(ROOT, "results/fgf10_sender_wilcox.csv"))
recv_summ   <- read.csv(file.path(ROOT, "results/fgfr2_receiver_summary.csv"))
recv_wilx   <- read.csv(file.path(ROOT, "results/fgfr2_receiver_wilcox.csv"))

# =============================================================================
# FIGURE 1 — Pipeline schematic + Epithelial UMAP + NicheNet ligand bar
# =============================================================================
message("\n=== Figure 1 ===")

# ── 1A: Pipeline schematic ───────────────────────────────────────────────────
schem_boxes <- tribble(
  ~x0, ~x1, ~y0, ~y1, ~label,                            ~fill,
  0.0,  1.4, 3.6, 4.4, "scRNA-seq\nData",                 "#DDDDDD",
  1.8,  3.2, 3.6, 4.4, "QC + Integration\n(Seurat/Harmony)", "#DDDDDD",
  3.6,  5.0, 3.6, 4.4, "Epithelial Fate\nLabels (Slingshot)", COL_FGF,
  0.0,  1.4, 1.6, 2.6, "Mouse NicheNet\n(ligand AUPR)",   COL_CONT,
  1.8,  3.2, 1.6, 2.6, "CellChat\n(comm. prob.)",          COL_CONT,
  3.6,  5.0, 1.6, 2.6, "Geneformer\n(in silico OE/DEL)",  COL_CONT,
  0.0,  1.4, 0.0, 1.0, "Human NicheNet\n(Habermann)",      COL_CONT,
  1.8,  5.0, 0.0, 1.0, "Composite Ranking\n→ Candidates",  COL_REGEN
)

schem_arrows <- tribble(
  ~x, ~xend, ~y,  ~yend,
  1.4, 1.8,   4.0, 4.0,   # Data → QC
  3.2, 3.6,   4.0, 4.0,   # QC → Fate
  0.7, 0.7,   3.6, 2.6,   # Fate → NicheNet
  2.5, 2.5,   3.6, 2.6,   # Fate → CellChat
  4.3, 4.3,   3.6, 2.6,   # Fate → Geneformer
  0.7, 0.7,   1.6, 1.0,   # NicheNet → Human
  0.7, 1.8,   0.5, 0.5,   # Human → Composite
  2.5, 2.5,   1.6, 1.0,   # CellChat → Composite
  4.3, 4.3,   1.6, 1.0    # Geneformer → Composite
)

p_schem <- ggplot() +
  geom_rect(data = schem_boxes,
            aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill),
            color = "black", linewidth = 0.6, alpha = 0.85) +
  scale_fill_identity() +
  geom_text(data = schem_boxes,
            aes(x = (x0 + x1) / 2, y = (y0 + y1) / 2, label = label),
            size = 3.5, fontface = "bold", color = "black") +
  geom_segment(data = schem_arrows,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               color = "grey30", linewidth = 0.7) +
  labs(title = "A   Analysis pipeline") +
  coord_cartesian(xlim = c(-0.2, 5.3), ylim = c(-0.3, 4.8)) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = TTL_SZ, color = "black",
                                  margin = margin(b = 4)))

# ── 1B: Epithelial fate UMAP ─────────────────────────────────────────────────
CELL_COLS <- c(
  "AT2"          = COL_FGF,    # warm orange — starting population
  "AT2 activated"= "#F5B97A",  # lighter orange — injury-activated
  "Krt8+ ADI"   = COL_REGEN,  # brick red — arrested/disease state
  "AT1"          = COL_CONT    # steel blue — regeneration destination
)
# Randomise plot order so no type covers another
set.seed(42)
umap_plot <- umap_df[sample(nrow(umap_df)), ]

p_umap <- ggplot(umap_plot, aes(UMAP1, UMAP2, color = cell_label)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(
    values = CELL_COLS,
    labels = c("AT2" = "AT2", "AT2 activated" = "AT2 activated",
               "Krt8+ ADI" = "Krt8⁺ ADI (arrested)", "AT1" = "AT1"),
    guide  = guide_legend(override.aes = list(size = 3, alpha = 1))
  ) +
  labs(title = "B   Epithelial cell states", x = "UMAP 1", y = "UMAP 2",
       color = NULL) +
  base_theme +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "right")

# ── 1C: NicheNet top-30 ligand bar, Fgf7 + Fgf10 highlighted ────────────────
nn_top30 <- nn_df %>%
  arrange(rank) %>%
  slice_head(n = 30) %>%
  mutate(
    highlight = case_when(
      test_ligand == "Fgf7"  ~ "Fgf7 (rank 3, FGFR2b)",
      test_ligand == "Fgf10" ~ "Fgf10 (rank 27)",
      TRUE                   ~ "other"
    ),
    bar_col = case_when(
      test_ligand == "Fgf7"  ~ COL_FGF,
      test_ligand == "Fgf10" ~ COL_REGEN,
      TRUE                   ~ COL_NS
    ),
    ligand_label = factor(test_ligand, levels = rev(test_ligand))
  )

p_nn <- ggplot(nn_top30, aes(x = aupr_corrected, y = ligand_label, fill = bar_col)) +
  geom_col(width = 0.75, color = NA) +
  scale_fill_identity() +
  geom_text(
    data = filter(nn_top30, highlight != "other"),
    aes(label = highlight, x = aupr_corrected + 0.003),
    hjust = 0, size = 3.8, fontface = "bold",
    color = ifelse(filter(nn_top30, highlight != "other")$test_ligand == "Fgf7",
                   COL_FGF, COL_REGEN)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35)),
                     name = "NicheNet AUPR (corrected)") +
  labs(title = "C   Top 30 sender ligands (mouse NicheNet)",
       subtitle = "Completing vs arrested fate (geneset_oi). Highlighted: FGFR2b ligands.",
       y = "Ligand") +
  base_theme +
  theme(axis.text.y = element_text(size = 11))

fig1 <- (p_schem / p_umap) | p_nn +
  plot_layout(widths = c(1, 1.1)) +
  plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))

save_fig(fig1, "Fig1_overview", w = 14, h = 9)

# =============================================================================
# FIGURE 2 — Agonist vs block receptor dot plot + composite score
# =============================================================================
message("\n=== Figure 2 ===")

CAND_SHOW <- c("Sdc1","Erbb2","Egfr","Fzd2","Sdc4","Fgfr2","Cd9",
               "Nrp1","Erbb3","Cd63","Bmpr2","App","Tgfbr2","Tgfbr1",
               "Itgb1","Axl","Ramp1","Notch2","Cd74")

dir_cols <- c(agonist        = COL_REGEN,
              agonist_uncertain = COL_FGF,
              block_TGFb     = COL_CONT,
              uncertain      = COL_UNCRT)

dir_labels <- c(agonist        = "Agonist",
                agonist_uncertain = "Agonist (uncertain)",
                block_TGFb     = "Block (TGFβ)",
                uncertain      = "Uncertain")

# ── 2A: Composite score lollipop ─────────────────────────────────────────────
score_cols <- c("score_nn","score_nn_h","score_cc","score_gf","score_del")

cand_plot <- cand_df %>%
  filter(receptor %in% CAND_SHOW) %>%
  mutate(
    n_methods    = rowSums(!is.na(across(all_of(score_cols)))),
    direction    = factor(therapeutic_direction, levels = names(dir_cols)),
    is_posctrl   = receptor == "Fgfr2",
    receptor_lbl = if_else(is_posctrl,
                           paste0(receptor, "*"),
                           receptor),
    receptor_lbl = factor(receptor_lbl,
                          levels = rev(receptor_lbl[order(composite_score)]))
  ) %>%
  arrange(composite_score)

p_lollipop <- ggplot(cand_plot,
                     aes(x = composite_score, y = receptor_lbl,
                         color = direction, size = n_methods)) +
  geom_segment(aes(xend = 0, yend = receptor_lbl), linewidth = 0.5,
               color = "grey80") +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = dir_cols, labels = dir_labels,
                     name = "Therapeutic\ndirection") +
  scale_size_continuous(name = "Supporting\nmethods",
                        range = c(3, 9), breaks = 1:5) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     name = "Composite score (0–1)") +
  annotate("text", x = 0.01, y = "Fgfr2*", hjust = 0,
           label = "* IIIb/IIIc isoform caveat", size = 3.2,
           color = COL_REGEN, fontface = "italic") +
  labs(title = "A   Receptor candidates: composite score",
       y = "Receptor") +
  base_theme +
  theme(axis.text.y = element_text(size = 12),
        legend.box  = "vertical")

# ── 2B: Method-by-method heatmap tile ────────────────────────────────────────
score_long <- cand_plot %>%
  select(receptor_lbl, direction, all_of(score_cols)) %>%
  pivot_longer(all_of(score_cols), names_to = "method", values_to = "score") %>%
  mutate(
    method = recode(method,
                    score_nn    = "Mouse\nNicheNet",
                    score_nn_h  = "Human\nNicheNet",
                    score_cc    = "CellChat",
                    score_gf    = "Geneformer\nOE",
                    score_del   = "Geneformer\nDEL"),
    method = factor(method, levels = c("Mouse\nNicheNet","Human\nNicheNet",
                                       "CellChat","Geneformer\nOE","Geneformer\nDEL"))
  )

p_heat <- ggplot(score_long, aes(x = method, y = receptor_lbl)) +
  geom_tile(aes(fill = score), color = "white", linewidth = 0.4) +
  geom_point(data = filter(score_long, is.na(score)),
             shape = 4, size = 2, color = "grey60") +
  scale_fill_gradient2(low = "white", mid = "#F5B97A", high = COL_REGEN,
                       midpoint = 0.5, na.value = NA,
                       name = "Score\n(0–1)", limits = c(0, 1)) +
  labs(title = "B   Per-method scores", x = NULL, y = NULL) +
  base_theme +
  theme(axis.text.y  = element_text(size = 11),
        axis.text.x  = element_text(size = 10, angle = 0),
        legend.position = "right",
        panel.grid    = element_blank())

fig2 <- p_lollipop | p_heat +
  plot_layout(widths = c(1.4, 1)) +
  plot_annotation(
    title = "Figure 2 — Ranked receptor candidates",
    theme = theme(plot.title = element_text(face = "bold", size = TTL_SZ + 1,
                                            color = "black"))
  )

save_fig(fig2, "Fig2_candidates", w = 14, h = 8)

# =============================================================================
# FIGURE 3 — FGF10 sender-side: mouse time course + human IPF vs Control
# =============================================================================
message("\n=== Figure 3 ===")

# ── 3A: Mouse Fgf10+ fibroblast fraction over bleomycin time course ──────────
strunz_tc <- sender_summ %>%
  filter(dataset == "Strunz_mouse_GSE141259",
         subtype  == "Fibroblasts (Plin2-low)") %>%
  mutate(
    phase = recode(condition,
                   PBS = "PBS", d3 = "d3", d7 = "d7",
                   d10 = "d10", d14 = "d14", d21 = "d21", d28 = "d28"),
    phase = factor(phase, levels = c("PBS","d3","d7","d10","d14","d21","d28")),
    bar_col = if_else(condition == "PBS", COL_CONT, COL_REGEN)
  )

# Wilcoxon annotation: fibrotic_early vs PBS p=0.029
# Show bracket over d10 + d14 bars vs PBS
p_mouse_tc <- ggplot(strunz_tc, aes(x = phase, y = pct_expressing, fill = bar_col)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.4, alpha = 0.85) +
  scale_fill_identity() +
  # p=0.029 bracket for fibrotic early (d10+d14)
  annotate("segment", x = 1, xend = 5.5, y = 11, yend = 11,
           linewidth = 0.7, color = "black") +
  annotate("segment", x = 1,   xend = 1,   y = 10.5, yend = 11,
           linewidth = 0.7, color = "black") +
  annotate("segment", x = 5.5, xend = 5.5, y = 10.5, yend = 11,
           linewidth = 0.7, color = "black") +
  annotate("text", x = 3.25, y = 11.6, label = "p = 0.029",
           size = 4, color = "black", fontface = "bold") +
  annotate("text", x = 1, y = 0.8, label = "0%", size = 3.5, color = "white",
           fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     name = "Fgf10+ fibroblasts (%)") +
  labs(title = "A   Mouse bleomycin time course: Fgf10⁺ Plin2-low fibroblasts",
       subtitle = "Strunz GSE141259. Pseudobulk Wilcoxon vs PBS. p = 0.029 for fibrotic early (d10+d14 pooled).",
       x = "Time point") +
  base_theme

# ── 3B: Human FGF10 pseudobulk — IPF vs Control (Adams + Habermann) ─────────
human_fgf10 <- sender_summ %>%
  filter(dataset %in% c("Adams_GSE136831", "Habermann_GSE135893"),
         condition %in% c("IPF", "Control"),
         subtype %in% c("Fibroblast (PLIN2-)", "Myofibroblast",
                        "Fibroblasts", "Myofibroblasts",
                        "PLIN2+ Fibroblast (lipofibroblast)",
                        "PLIN2+ Fibroblasts", "HAS1 High Fibroblasts")) %>%
  mutate(
    ds_label = recode(dataset,
                      Adams_GSE136831     = "Adams",
                      Habermann_GSE135893 = "Habermann"),
    subtype_short = subtype %>%
      gsub("Fibroblast \\(PLIN2-\\)",              "PLIN2− Fib",    .) %>%
      gsub("Fibroblasts$",                         "Fibroblasts",   .) %>%
      gsub("Myofibroblast$",                       "Myofib",        .) %>%
      gsub("Myofibroblasts$",                      "Myofib",        .) %>%
      gsub("PLIN2\\+ Fibroblast \\(lipofibroblast\\)", "PLIN2+ Fib",.  ) %>%
      gsub("PLIN2\\+ Fibroblasts",                 "PLIN2+ Fib",    .) %>%
      gsub("HAS1 High Fibroblasts",                "HAS1-high Fib", .),
    panel_label  = paste0(ds_label, "\n", subtype_short),
    cond_col     = if_else(condition == "IPF", COL_REGEN, COL_CONT),
    condition    = factor(condition, levels = c("Control","IPF"))
  )

# Attach p-values from wilcox CSV
wlx_map <- sender_wilx %>%
  filter(dataset %in% c("Adams_GSE136831","Habermann_GSE135893")) %>%
  mutate(p_label = sapply(p_value, fmt_p))

# Add significance labels to human_fgf10
human_fgf10_ipf <- human_fgf10 %>% filter(condition == "IPF") %>%
  left_join(
    wlx_map %>% select(dataset, subtype, p_label),
    by = c("dataset", "subtype")
  )

p_human_fgf10 <- ggplot(human_fgf10,
                         aes(x = condition, y = pct_expressing,
                             fill = condition)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.4, alpha = 0.85) +
  geom_text(
    data = human_fgf10_ipf,
    aes(x = 2, y = pct_expressing + 1.5, label = p_label),
    size = 3.5, color = "black", fontface = "bold", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c(Control = COL_CONT, IPF = COL_REGEN),
                    guide   = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22)),
                     name = "FGF10+ cells (%)") +
  facet_wrap(~ panel_label, nrow = 2, scales = "free_y") +
  labs(title = "B   Human fibroblasts: FGF10 expression",
       subtitle = "Adams GSE136831 + Habermann GSE135893. Pseudobulk Wilcoxon. n.s. = p ≥ 0.05.",
       x = NULL) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 11))

fig3 <- p_mouse_tc / p_human_fgf10 +
  plot_layout(heights = c(1, 1.4)) +
  plot_annotation(
    title = "Figure 3 — FGF10 sender-side: supply is episodic in mouse and low in human IPF",
    theme = theme(plot.title = element_text(face = "bold", size = TTL_SZ + 1,
                                            color = "black"))
  )

save_fig(fig3, "Fig3_sender", w = 12, h = 11)

# =============================================================================
# FIGURE 4 — FGFR2 receiver: AT2 expresser % + co-expression with Krt8/pEMT
# =============================================================================
message("\n=== Figure 4 ===")

# ── 4A: AT2 FGFR2 expresser % ────────────────────────────────────────────────
# Left panel: mouse bleomycin phases (epi_obj)
mouse_fgfr2_pct <- recv_summ %>%
  filter(dataset == "Mouse_epi_obj", gene == "Fgfr2") %>%
  mutate(
    phase = factor(condition, levels = c("PBS_control","acute",
                                         "fibrotic_early","fibrotic_late","resolution")),
    phase_lbl = recode(condition,
                       PBS_control    = "PBS",
                       acute          = "Acute",
                       fibrotic_early = "Fibrotic\nearly",
                       fibrotic_late  = "Fibrotic\nlate",
                       resolution     = "Resolution")
  )

# All mouse Wilcoxon comparisons are n.s. (p > 0.38)
p_mouse_fgfr2 <- ggplot(mouse_fgfr2_pct,
                          aes(x = factor(phase_lbl,
                                         levels = c("PBS","Acute","Fibrotic\nearly",
                                                    "Fibrotic\nlate","Resolution")),
                              y = pct_expressing)) +
  geom_col(aes(fill = condition == "PBS_control"),
           width = 0.72, color = "black", linewidth = 0.4, alpha = 0.85) +
  scale_fill_manual(values = c(`TRUE` = COL_CONT, `FALSE` = COL_FGF),
                    guide = "none") +
  geom_text(aes(y = pct_expressing + 0.5, label = "n.s."),
            size = 3.5, color = "grey40") +
  # Remove n.s. from PBS itself
  scale_y_continuous(limits = c(0, 20),
                     expand = expansion(mult = c(0, 0.1)),
                     name = "Fgfr2⁺ AT2 cells (%)") +
  labs(title = "A   Mouse AT2: Fgfr2⁺ fraction over bleomycin time course",
       subtitle = "epi_obj (Kobayashi PATS). All comparisons vs PBS: n.s. (p > 0.38).",
       x = "Phase") +
  base_theme

# Right panel: human Adams ATII IPF vs Control
human_fgfr2 <- recv_summ %>%
  filter(dataset == "Adams_GSE136831", gene == "FGFR2",
         subtype == "ATII", condition %in% c("IPF","Control"))

fgfr2_wilx_h <- recv_wilx %>%
  filter(dataset == "Adams_GSE136831", gene == "FGFR2", subtype == "ATII")

p_val_fgfr2_h <- if (nrow(fgfr2_wilx_h) > 0) {
  p_raw <- fgfr2_wilx_h$p_value[1]
  if (is.na(p_raw) || p_raw >= 0.05) "n.s." else fmt_p(p_raw)
} else "n.s."

p_human_fgfr2 <- ggplot(human_fgfr2,
                          aes(x = factor(condition, levels = c("Control","IPF")),
                              y = pct_expressing,
                              fill = condition)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.4, alpha = 0.85) +
  scale_fill_manual(values = c(Control = COL_CONT, IPF = COL_REGEN),
                    guide = "none") +
  # Significance bracket
  annotate("segment", x = 1, xend = 2, y = 50, yend = 50,
           linewidth = 0.7, color = "black") +
  annotate("segment", x = 1, xend = 1, y = 48, yend = 50,
           linewidth = 0.7, color = "black") +
  annotate("segment", x = 2, xend = 2, y = 48, yend = 50,
           linewidth = 0.7, color = "black") +
  annotate("text", x = 1.5, y = 52,
           label = p_val_fgfr2_h,
           size = 4, color = "black", fontface = "bold") +
  scale_y_continuous(limits = c(0, 56),
                     expand = expansion(mult = c(0, 0.05)),
                     name = "FGFR2⁺ ATII cells (%)") +
  labs(title = "Human Adams ATII:\nFGFR2 IPF vs Control",
       subtitle = "Adams GSE136831.",
       x = NULL) +
  base_theme

# ── 4B: FGFR2 / Krt8-pEMT co-expression (epi_obj AT2 cells) ─────────────────
# Rebuild from epi_obj metadata (loaded above into umap_df)
# Keep only AT2 + AT2 activated cells in bleomycin phases
coexpr_df <- umap_df %>%
  filter(cell_label %in% c("AT2","AT2 activated"),
         phase %in% PHASE_ORDER) %>%
  mutate(
    Fgfr2_count  = fgfr2_cnt[rownames(.)],
    Fgfr2_status = if_else(!is.na(Fgfr2_count) & Fgfr2_count > 0,
                            "Fgfr2+", "Fgfr2-"),
    phase = factor(phase, levels = PHASE_ORDER)
  )

if (!is.null(trans_col) && trans_col %in% colnames(coexpr_df)) {
  coexpr_df$arrest_score <- coexpr_df[[trans_col]]
} else {
  coexpr_df$arrest_score <- NA_real_
}
coexpr_df <- filter(coexpr_df, !is.na(arrest_score))

p_coexpr <- ggplot(coexpr_df,
                    aes(x = Fgfr2_status, y = arrest_score, fill = Fgfr2_status)) +
  geom_violin(scale = "width", alpha = 0.75, trim = TRUE, linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.shape = NA,
               fill = "white", alpha = 0.85, linewidth = 0.5) +
  scale_fill_manual(values = c(`Fgfr2+` = COL_FGF, `Fgfr2-` = COL_NS),
                    guide = "none") +
  facet_wrap(~ phase, nrow = 1) +
  scale_y_continuous(name = "Krt8/pEMT arrest score\n(Transitional1 module)") +
  labs(
    title    = "B   Co-expression: Fgfr2 and Krt8/pEMT arrest score in AT2 cells",
    subtitle = paste0(
      "epi_obj AT2 + AT2-activated cells. Violin = score distribution in Fgfr2+ vs Fgfr2− cells.\n",
      "FGFR2 IIIb/IIIc isoform caveat: total Fgfr2 used as FGFR2b surrogate (AT2 cells predominantly IIIb)."
    ),
    x = NULL
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig4 <- (p_mouse_fgfr2 | p_human_fgfr2) / p_coexpr +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    title = "Figure 4 — FGFR2 receiver: expression maintained; co-expressed with arrest markers",
    theme = theme(plot.title = element_text(face = "bold", size = TTL_SZ + 1,
                                            color = "black"))
  )

save_fig(fig4, "Fig4_receiver", w = 14, h = 10)

# =============================================================================
# COMBINED PDF (one figure per page)
# =============================================================================
message("\nSaving combined PDF …")
pdf_path <- file.path(OUT_FIG, "manuscript_figures.pdf")
cairo_pdf(pdf_path, width = 14, height = 10, onefile = TRUE)
print(fig1)
print(fig2)
print(fig3)
print(fig4)
dev.off()
message(sprintf("  Saved: manuscript_figures.pdf"))

message("\n09_manuscript_figures.R complete.")
message("Outputs:")
message("  results/figures/Fig1_overview.png")
message("  results/figures/Fig2_candidates.png")
message("  results/figures/Fig3_sender.png")
message("  results/figures/Fig4_receiver.png")
message("  results/figures/manuscript_figures.pdf")
