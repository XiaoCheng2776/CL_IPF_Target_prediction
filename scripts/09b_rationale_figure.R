#!/usr/bin/env Rscript
# =============================================================================
# 09b_rationale_figure.R
# Fig 1 — Motivating the in-silico screen: AT2 bifurcation, arrest accumulation,
#          and Fgfr2 retention on stalled cells.
#
# Three panels (all from results/markers/epi_obj.rds):
#   A. UMAP: epithelial cell states + bifurcation arrows
#      completing → AT1 (Ager/Pdpn/Hopx)  vs  arrested → Krt8+/Cldn4+
#   B. Line: fraction of Krt8+ ADI cells across bleomycin time course (d2–d21)
#      day28 excluded (n=9); day36/54 shown as grey open points
#   C. Three UMAP feature plots: Krt8, Cldn4, Fgfr2
#      visual message: arrested cells express Krt8/Cldn4 AND retain Fgfr2
#
# Style: #E07B39/#B83232 (FGF/regen), #4472B8 (contrast),
#        theme_classic, black text, 300 dpi PNG + appended to manuscript PDF.
#
# Output:
#   results/figures/Fig1_rationale.png
#   results/figures/manuscript_figures.pdf  (Fig 1 prepended)
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

COL_FGF   <- "#E07B39"
COL_REGEN <- "#B83232"
COL_CONT  <- "#4472B8"
COL_NS    <- "#D3D3D3"
COL_DARK  <- "#222222"

AX_SZ  <- 13
TTL_SZ <- 14
DPI    <- 300

base_theme <- theme_classic(base_size = 13) +
  theme(
    axis.text    = element_text(size = AX_SZ, color = "black"),
    axis.title   = element_text(size = AX_SZ, color = "black"),
    plot.title   = element_text(size = TTL_SZ, color = "black", face = "bold"),
    plot.subtitle = element_text(size = AX_SZ - 1, color = "black"),
    legend.text  = element_text(size = AX_SZ - 1, color = "black"),
    legend.title = element_text(size = AX_SZ - 1, color = "black", face = "bold"),
    strip.text   = element_text(size = AX_SZ, color = "black", face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA)
  )

feat_theme <- base_theme +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    axis.title = element_text(size = 10, color = "black")
  )

# ─────────────────────────────────────────────────────────────────────────────
# Load epi_obj and extract everything needed up front
# ─────────────────────────────────────────────────────────────────────────────
message("Loading epi_obj.rds …")
epi  <- readRDS(file.path(ROOT, "results/markers/epi_obj.rds"))
meta <- epi@meta.data

# Detect column names defensively
ct_col    <- if ("cell_type"   %in% colnames(meta)) "cell_type"   else "cell.type"
tp_col    <- if ("time_point"  %in% colnames(meta)) "time_point"  else "timepoint"
pt_col    <- if ("pseudotime"  %in% colnames(meta)) "pseudotime"  else NULL
fl_col    <- if ("fate_label"  %in% colnames(meta)) "fate_label"  else NULL

# UMAP coords (stored as umap_1/umap_2 in metadata for this object)
if ("umap_1" %in% colnames(meta)) {
  umap_df <- meta[, c("umap_1","umap_2"), drop=FALSE]
  colnames(umap_df) <- c("UMAP1","UMAP2")
} else {
  umap_df <- as.data.frame(Embeddings(epi, "umap"))
  colnames(umap_df) <- c("UMAP1","UMAP2")
}

# SCT gene expression for feature panels
feat_genes <- intersect(c("Krt8","Cldn4","Fgfr2","Ager","Pdpn","Hopx"),
                        rownames(epi))
message(sprintf("  Extracting genes: %s", paste(feat_genes, collapse=", ")))
expr_mat <- GetAssayData(epi, assay = "SCT", layer = "counts")[feat_genes, ]

# Assemble working data frame
df <- cbind(
  umap_df,
  meta[rownames(umap_df), c(ct_col, tp_col,
                             if (!is.null(pt_col)) pt_col else character(0),
                             if (!is.null(fl_col)) fl_col else character(0),
                             "RegenSuccess1","Transitional1")]
)
colnames(df)[3] <- "cell_type"
colnames(df)[4] <- "time_point"
if (!is.null(pt_col))  colnames(df)[5] <- "pseudotime"
if (!is.null(fl_col))  colnames(df)[if (!is.null(pt_col)) 6 else 5] <- "fate_label"

for (g in feat_genes) {
  df[[g]] <- as.numeric(expr_mat[g, rownames(df)])
}

# Winsorise gene expression at 99th percentile for feature plots
for (g in feat_genes) {
  hi <- quantile(df[[g]], 0.99, na.rm = TRUE)
  df[[paste0(g,"_w")]] <- pmin(df[[g]], hi)
}

message("  epi_obj loaded. Freeing Seurat object …")
rm(epi, expr_mat); gc()

# ─────────────────────────────────────────────────────────────────────────────
# PANEL A — UMAP: epithelial cell states + bifurcation annotation
# ─────────────────────────────────────────────────────────────────────────────
CELL_COLS <- c(
  "AT2"           = COL_FGF,   # warm orange  — progenitor (starting)
  "AT2 activated" = "#F5B97A", # light orange — injury-activated
  "Krt8+ ADI"    = COL_REGEN, # brick red    — arrested (therapeutic target)
  "AT1"           = COL_CONT   # steel blue   — regeneration destination
)
CELL_LABELS <- c(
  "AT2"           = "AT2",
  "AT2 activated" = "AT2 activated",
  "Krt8+ ADI"    = "Krt8⁺ ADI\n(arrested)",
  "AT1"           = "AT1\n(destination)"
)

# Cluster centroids for label placement
centroids <- df %>%
  group_by(cell_type) %>%
  summarise(u1 = median(UMAP1), u2 = median(UMAP2), .groups = "drop")

# Shuffle point order
set.seed(42)
df_shuf <- df[sample(nrow(df)), ]

p_a <- ggplot(df_shuf, aes(UMAP1, UMAP2)) +
  # Background: all cells grey first, then coloured on top
  geom_point(data = filter(df_shuf, cell_type == "AT2 activated"),
             aes(color = cell_type), size = 0.4, alpha = 0.55) +
  geom_point(data = filter(df_shuf, cell_type == "AT2"),
             aes(color = cell_type), size = 0.4, alpha = 0.55) +
  geom_point(data = filter(df_shuf, cell_type %in% c("Krt8+ ADI","AT1")),
             aes(color = cell_type), size = 0.5, alpha = 0.75) +
  scale_color_manual(values = CELL_COLS, labels = CELL_LABELS,
                     guide = guide_legend(override.aes = list(size = 3, alpha = 1),
                                          ncol = 1)) +
  # Cluster labels (white box behind text for readability)
  geom_label(data = centroids,
             aes(x = u1, y = u2, label = cell_type),
             fill = "white", color = "black", alpha = 0.85,
             size = 3.2, linewidth = 0.3, fontface = "bold",
             label.padding = unit(0.18, "lines")) +
  # Trajectory arrows: AT2 → AT1 (completing) and AT2 → Krt8+ADI (arrested)
  annotate("curve",
           x    = centroids$u1[centroids$cell_type == "AT2"],
           xend = centroids$u1[centroids$cell_type == "AT1"],
           y    = centroids$u2[centroids$cell_type == "AT2"],
           yend = centroids$u2[centroids$cell_type == "AT1"],
           curvature = -0.30, linewidth = 1.1, color = COL_CONT,
           arrow = arrow(length = unit(0.22, "cm"), type = "closed")) +
  annotate("curve",
           x    = centroids$u1[centroids$cell_type == "AT2"],
           xend = centroids$u1[centroids$cell_type == "Krt8+ ADI"],
           y    = centroids$u2[centroids$cell_type == "AT2"],
           yend = centroids$u2[centroids$cell_type == "Krt8+ ADI"],
           curvature = 0.35, linewidth = 1.1, color = COL_REGEN,
           arrow = arrow(length = unit(0.22, "cm"), type = "closed")) +
  # Fate labels on the trajectories
  annotate("text",
           x = mean(c(centroids$u1[centroids$cell_type=="AT2"],
                      centroids$u1[centroids$cell_type=="AT1"])) - 0.5,
           y = mean(c(centroids$u2[centroids$cell_type=="AT2"],
                      centroids$u2[centroids$cell_type=="AT1"])) + 1.2,
           label = "Completing\n→ AT1",
           color = COL_CONT, size = 3.6, fontface = "bold", hjust = 0.5) +
  annotate("text",
           x = mean(c(centroids$u1[centroids$cell_type=="AT2"],
                      centroids$u1[centroids$cell_type=="Krt8+ ADI"])) + 0.3,
           y = mean(c(centroids$u2[centroids$cell_type=="AT2"],
                      centroids$u2[centroids$cell_type=="Krt8+ ADI"])) - 1.1,
           label = "Arrested\n(therapeutic target)",
           color = COL_REGEN, size = 3.6, fontface = "bold", hjust = 0.5) +
  labs(title = "A   Epithelial cell-state landscape",
       subtitle = "Mouse bleomycin model (Kobayashi PATS). AT2 progenitors bifurcate: completing or arrested.",
       x = "UMAP 1", y = "UMAP 2", color = "Cell type") +
  base_theme +
  feat_theme +
  theme(axis.title = element_text(size = 10),
        legend.position = "right",
        plot.subtitle = element_text(size = 10))

# ─────────────────────────────────────────────────────────────────────────────
# PANEL B — Arrested-cell fraction over bleomycin time course
# ─────────────────────────────────────────────────────────────────────────────
time_df <- df %>%
  mutate(day = suppressWarnings(as.numeric(gsub(".*?([0-9]+).*","\\1", time_point)))) %>%
  filter(!grepl("PBS", time_point, ignore.case = TRUE),
         !is.na(day)) %>%
  group_by(time_point, day) %>%
  summarise(
    n_total    = n(),
    n_arrested = sum(cell_type == "Krt8+ ADI"),
    frac_pct   = 100 * n_arrested / n_total,
    .groups    = "drop"
  ) %>%
  arrange(day) %>%
  mutate(
    low_n  = n_total < 50,   # day28: n=9
    include = day <= 21      # main series: d2-d21
  )

main_df  <- filter(time_df, include)
lown_df  <- filter(time_df, !include, !low_n)   # d36, d54
excl_df  <- filter(time_df, low_n)              # d28 (n=9)

# Peak annotation
peak_row <- main_df[which.max(main_df$frac_pct), ]

p_b <- ggplot(main_df, aes(x = day, y = frac_pct)) +
  # Low-n excluded days as open grey points (d36, d54)
  geom_point(data = lown_df, aes(x = day, y = frac_pct),
             shape = 1, size = 3, color = "grey60", stroke = 1) +
  geom_text(data = lown_df, aes(x = day, y = frac_pct + 2.5,
            label = paste0("d", day, "\n(n=", n_total, ")")),
            size = 2.8, color = "grey50", hjust = 0.5) +
  # Main line + area
  geom_area(fill = COL_REGEN, alpha = 0.15) +
  geom_line(color = COL_REGEN, linewidth = 1.1) +
  geom_point(aes(size = n_total), color = COL_REGEN, fill = "white",
             shape = 21, stroke = 1.5) +
  scale_size_continuous(name = "n cells", range = c(2, 5),
                        breaks = c(500, 1000, 1500)) +
  # Peak annotation
  annotate("segment",
           x = peak_row$day, xend = peak_row$day,
           y = peak_row$frac_pct + 1.5, yend = peak_row$frac_pct + 7,
           color = COL_REGEN, linewidth = 0.7, linetype = "dashed") +
  annotate("text", x = peak_row$day, y = peak_row$frac_pct + 9,
           label = sprintf("Peak: %.0f%%\n(day %d)", peak_row$frac_pct, peak_row$day),
           size = 3.4, color = COL_REGEN, fontface = "bold", hjust = 0.5) +
  # Excluded day28 footnote
  annotate("text", x = 21.5, y = 2,
           label = "day 28 excluded\n(n=9)", size = 2.6,
           color = "grey50", hjust = 0, fontface = "italic") +
  scale_x_continuous(breaks = c(2,4,6,8,10,12,15,21),
                     name = "Day after bleomycin") +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0.02, 0.18)),
                     name = "Krt8⁺ ADI cells (%)") +
  labs(title = "B   Arrest accumulates during fibrotic consolidation",
       subtitle = "Fraction of Krt8⁺ ADI cells per time point. n=9 at day 28 excluded.") +
  base_theme +
  theme(legend.position = c(0.88, 0.75),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.5,"cm"))

# ─────────────────────────────────────────────────────────────────────────────
# PANEL C — Feature plots: Krt8, Cldn4, Fgfr2 on UMAP
# ─────────────────────────────────────────────────────────────────────────────
# Custom colour scale: grey → orange → brick red (matches palette)
feat_scale <- function(gene, title_str = NULL) {
  hi_pct <- quantile(df[[paste0(gene,"_w")]], 0.95, na.rm = TRUE)
  scale_color_gradientn(
    colors   = c("grey92", "#F5B97A", COL_FGF, COL_REGEN),
    values   = rescale(c(0, hi_pct * 0.15, hi_pct * 0.5, hi_pct)),
    limits   = c(0, max(df[[paste0(gene,"_w")]], na.rm=TRUE)),
    na.value = "grey92",
    name     = if (!is.null(title_str)) title_str else gene,
    guide    = guide_colorbar(barwidth = 0.5, barheight = 3.5,
                              title.position = "top", title.hjust = 0.5)
  )
}

make_feat_plot <- function(gene, panel_letter, subtitle_txt) {
  # Sort: low expressers plotted first (high expressers on top)
  expr_col <- paste0(gene,"_w")
  d <- df[order(df[[expr_col]]), ]
  n_expr <- sum(df[[gene]] > 0, na.rm = TRUE)
  pct_expr <- 100 * n_expr / nrow(df)

  ggplot(d, aes(x = UMAP1, y = UMAP2, color = .data[[expr_col]])) +
    geom_point(size = 0.3, alpha = 0.7) +
    feat_scale(gene) +
    labs(
      title    = sprintf("%s   %s", panel_letter, subtitle_txt),
      subtitle = sprintf("%d / %d cells express (%.0f%%)", n_expr, nrow(df), pct_expr),
      x = "UMAP 1", y = "UMAP 2"
    ) +
    base_theme +
    feat_theme +
    theme(
      plot.title    = element_text(size = TTL_SZ, face = "bold"),
      plot.subtitle = element_text(size = 9),
      legend.position = "right"
    )
}

p_krt8  <- make_feat_plot("Krt8",  "C", "Krt8 (arrest marker)")
p_cldn4 <- make_feat_plot("Cldn4", " ", "Cldn4 (arrest marker)")
p_fgfr2 <- make_feat_plot("Fgfr2", " ", "Fgfr2 (target receptor)")

# Annotate Krt8+ADI region on the Fgfr2 plot with a rectangle or text
# Use centroid of Krt8+ADI cells
krt8adi_c <- df %>%
  filter(cell_type == "Krt8+ ADI") %>%
  summarise(u1 = median(UMAP1), u2 = median(UMAP2))

p_fgfr2 <- p_fgfr2 +
  annotate("text",
           x = krt8adi_c$u1, y = krt8adi_c$u2 + 1.5,
           label = "Receptor retained\nin arrested cells",
           color = COL_REGEN, size = 3.0, fontface = "bold", hjust = 0.5) +
  annotate("segment",
           x = krt8adi_c$u1, xend = krt8adi_c$u1,
           y = krt8adi_c$u2 + 0.9, yend = krt8adi_c$u2 + 0.2,
           color = COL_REGEN, linewidth = 0.7,
           arrow = arrow(length = unit(0.15,"cm"), type="closed"))

# ─────────────────────────────────────────────────────────────────────────────
# ASSEMBLE & SAVE
# ─────────────────────────────────────────────────────────────────────────────
fig_rationale <- (p_a | p_b) / (p_krt8 | p_cldn4 | p_fgfr2) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    title = paste0(
      "Figure 1 — AT2 regeneration arrests during bleomycin-induced fibrosis; ",
      "Fgfr2 is retained on stalled cells"
    ),
    caption = paste0(
      "Kobayashi PATS (GSE141634) + Choi DATPs (GSE145031). ",
      "Fgfr2 IIIb/IIIc isoform caveat: AT2 cells express predominantly IIIb; ",
      "total Fgfr2 used as surrogate."
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", size = TTL_SZ + 1, color = "black"),
      plot.caption = element_text(size = 9, color = "grey40", hjust = 0)
    )
  )

out_png <- file.path(OUT_FIG, "Fig1_rationale.png")
ggsave(out_png, fig_rationale, width = 15, height = 10, dpi = DPI)
message(sprintf("Saved: %s", out_png))

# Prepend to combined PDF (write as new file with this figure first)
pdf_path <- file.path(OUT_FIG, "manuscript_figures.pdf")
cairo_pdf(pdf_path, width = 15, height = 10, onefile = TRUE)
print(fig_rationale)
dev.off()
message(sprintf("Updated: %s (rationale figure only — rerun 09_manuscript_figures.R to append full set)", pdf_path))

message("\n09b_rationale_figure.R complete.")
message(sprintf("  %s", out_png))
