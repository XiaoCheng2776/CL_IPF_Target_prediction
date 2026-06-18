#!/usr/bin/env Rscript
# 09c_synthesis_figure.R
# Conserved-target synthesis figure: lollipop chart of IPF therapeutic candidates
# Main panel: conserved targets (human NicheNet evidence)
# Grey section: mouse-only candidates (pending human validation)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

COL_AGONIST <- "#E07B39"   # warm orange — predicted stimulate
COL_BLOCK   <- "#4472B8"   # steel blue — inhibit
COL_POSCTRL <- "#B83232"   # brick red — positive control
COL_MOUSE   <- "#9E9E9E"   # grey — mouse-only

out_dir <- "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ─── 1. Load data ────────────────────────────────────────────────────────────
df <- read.csv("results/candidates_ranked.csv", stringsAsFactors = FALSE)

fmt_rank <- function(r) {
  ifelse(is.na(r) | r == "", "—", paste0("#", as.integer(r)))
}

get_row <- function(t, df) {
  row <- df[tolower(df$receptor) == tolower(t), ]
  if (nrow(row) == 0) stop(paste("Target not found:", t))
  row[1, ]
}

# ─── 2. Define target sets ───────────────────────────────────────────────────
conserved_agonist <- c("Sdc1", "Erbb2", "Egfr", "Sdc4", "Fgfr2", "Cd9")
conserved_block   <- c("App", "Tgfbr2")
mouse_only_show   <- c("Cd74", "Ramp1", "Fzd2", "Bmpr2", "Nrp1", "Tgfbr1")

d_ag <- bind_rows(lapply(conserved_agonist, get_row, df)) |> arrange(composite_score)
d_bl <- bind_rows(lapply(conserved_block,   get_row, df)) |> arrange(composite_score)
d_mo <- bind_rows(lapply(mouse_only_show,   get_row, df)) |> arrange(composite_score)

# ─── 3. Y positions (ascending = visually bottom-to-top) ─────────────────────
# Layout (visual top → bottom): agonist arm, divider, block arm, divider, mouse section
n_ag <- nrow(d_ag); n_bl <- nrow(d_bl); n_mo <- nrow(d_mo)

y_mo_seq <- seq(1, by = 1, length.out = n_mo)                        # 1–6
y_bl_seq <- seq(max(y_mo_seq) + 3, by = 1, length.out = n_bl)        # 9–10
y_ag_seq <- seq(max(y_bl_seq) + 3, by = 1, length.out = n_ag)        # 13–18

d_ag$y <- y_ag_seq
d_bl$y <- y_bl_seq
d_mo$y <- y_mo_seq

# ─── 4. Colours, labels, display names ───────────────────────────────────────
d_ag$col      <- ifelse(tolower(d_ag$receptor) == "fgfr2", COL_POSCTRL, COL_AGONIST)
d_bl$col      <- COL_BLOCK
d_mo$col      <- COL_MOUSE

d_ag$fontface_col <- ifelse(tolower(d_ag$receptor) == "fgfr2", "bold", "plain")
d_bl$fontface_col <- "plain"
d_mo$fontface_col <- "plain"

d_ag$display  <- d_ag$receptor
d_ag$display[tolower(d_ag$receptor) == "fgfr2"] <- "Fgfr2 ★ (PC)"

d_bl$display  <- d_bl$receptor

d_mo$display  <- paste0(
  d_mo$receptor,
  ifelse(d_mo$therapeutic_direction == "block_TGFb", " (↓)",
  ifelse(grepl("agonist", d_mo$therapeutic_direction), " (↑?)", " (?)"))
)

d_ag$pt_size  <- 5.0;  d_bl$pt_size <- 5.0;  d_mo$pt_size <- 3.5

all_df <- bind_rows(d_ag, d_bl, d_mo)
all_df$rank_text <- paste0(
  fmt_rank(all_df$nn_best_ligand_rank), " / ",
  fmt_rank(all_df$nn_human_best_rank)
)

# ─── 5. Reference coordinates ────────────────────────────────────────────────
div_bl_mo  <- (max(y_mo_seq) + min(y_bl_seq)) / 2      # between mouse + block sections
div_ag_bl  <- (max(y_bl_seq) + min(y_ag_seq)) / 2      # between block + agonist sections
hdr_y_mo   <- div_bl_mo - 0.55   # label just below divider (reads: divider → label → targets)
hdr_y_bl   <- div_ag_bl - 0.55
hdr_y_ag   <- max(y_ag_seq) + 0.85  # above agonist arm (topmost label in figure)

y_min_plot <- min(all_df$y) - 0.9
y_max_plot <- max(all_df$y) + 2.5

# Mouse-only background rectangle
bg_mo <- data.frame(
  xmin = -0.01, xmax = 1.02,
  ymin = min(y_mo_seq) - 0.65,
  ymax = max(y_mo_seq) + 0.65
)

x_lo <- -0.62   # leave room for receptor labels
x_hi <-  1.16

# ─── 6. Lollipop panel ───────────────────────────────────────────────────────
p_lollipop <- ggplot(all_df) +
  # Mouse-only shaded background
  geom_rect(data = bg_mo,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#F0F0F0", color = "#BBBBBB", linewidth = 0.5) +
  # Vertical reference at x = 0
  geom_vline(xintercept = 0, color = "#CCCCCC", linewidth = 0.4) +
  # Section dividers
  geom_hline(yintercept = c(div_bl_mo, div_ag_bl),
             linetype = "dashed", color = "#BBBBBB", linewidth = 0.5) +
  # Lollipop stems
  geom_segment(aes(x = 0, xend = composite_score, y = y, yend = y, color = col),
               linewidth = 1.1, show.legend = FALSE) +
  # Lollipop points
  geom_point(aes(x = composite_score, y = y, color = col, size = pt_size),
             show.legend = FALSE) +
  scale_size_identity() +
  # Composite score value (right of point)
  geom_text(aes(x = composite_score + 0.03, y = y,
                label = sprintf("%.3f", composite_score)),
            hjust = 0, size = 3.7, color = "black") +
  # Receptor name labels (left margin)
  geom_text(aes(x = -0.04, y = y, label = display,
                color = col, fontface = fontface_col),
            hjust = 1, size = 4.5) +
  # Section headers (above each arm / section)
  annotate("text", x = 0.5, y = hdr_y_ag,
           label = "PREDICTED STIMULATE  (Agonist / Pro-regenerative)",
           color = COL_AGONIST, fontface = "bold", size = 4.8, hjust = 0.5) +
  annotate("text", x = 0.5, y = hdr_y_bl,
           label = paste0("INHIBIT  (Block TGFβ → Krt8⁺ arrest)"),
           color = COL_BLOCK, fontface = "bold", size = 4.8, hjust = 0.5) +
  annotate("text", x = 0.5, y = hdr_y_mo,
           label = "Mouse-only — human validation needed",
           color = "#666666", fontface = "bold.italic", size = 4.1, hjust = 0.5) +
  # Colour scale
  scale_color_identity() +
  # Axes
  scale_x_continuous(
    name   = "Composite score",
    limits = c(x_lo, x_hi),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_y_continuous(
    limits = c(y_min_plot, y_max_plot),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.title.y       = element_blank(),
    axis.line.y        = element_blank(),
    axis.text.x        = element_text(color = "black", size = 13),
    axis.title.x       = element_text(color = "black", size = 13),
    panel.grid.major.x = element_line(color = "#EBEBEB", linewidth = 0.3),
    legend.position    = "none",
    plot.margin        = margin(t = 10, r = 5, b = 10, l = 5)
  )

# ─── 7. Rank annotation panel (right) ────────────────────────────────────────
p_rank <- ggplot(all_df, aes(y = y)) +
  geom_rect(data = bg_mo,
            aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#F0F0F0", color = "#BBBBBB", linewidth = 0.5) +
  geom_hline(yintercept = c(div_bl_mo, div_ag_bl),
             linetype = "dashed", color = "#BBBBBB", linewidth = 0.5) +
  geom_text(aes(x = 0.5, label = rank_text, color = col),
            size = 4.0, hjust = 0.5) +
  annotate("text", x = 0.5, y = max(all_df$y) + 1.6,
           label = "Mouse NN /\nHuman NN",
           fontface = "bold", size = 3.9, color = "black",
           hjust = 0.5, lineheight = 1.2) +
  scale_color_identity() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(y_min_plot, y_max_plot), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 0))

# ─── 8. Combine and annotate ─────────────────────────────────────────────────
caption_txt <- paste0(
  "★ PC = positive control: FGF10→FGFR2b axis validated in vivo (Fgf10-deficient mice develop fibrosis more readily). ",
  "Large points = conserved (mouse + human evidence); small = mouse-only.\n",
  "Agonist arm: direction predicted from NicheNet ligand-activity scores + receptor biology; ",
  "Geneformer in silico perturbation pending GPU rerun.\n",
  "Block arm: TGFβ-driven Krt8⁺ arrest pathway; direction supported by mouse + human NicheNet and Geneformer deletion scores.\n",
  "Composite score: weighted average of NicheNet-mouse, NicheNet-human, CellChat, Geneformer-OE, and Geneformer-DEL subscores (0–1 scale).\n",
  "Mouse NN rank: top ligand-activity rank from Strunz NicheNet. Human NN rank: top rank from Adams/Habermann NicheNet."
)

fig <- (p_lollipop | p_rank) +
  plot_layout(widths = c(4.5, 1)) +
  plot_annotation(
    title   = "IPF In Silico Target Screen: Conserved Therapeutic Candidates",
    caption = caption_txt,
    theme = theme(
      plot.title   = element_text(size = 17, face = "bold", color = "black",
                                  margin = margin(b = 10)),
      plot.caption = element_text(size = 8.5, color = "grey35", hjust = 0,
                                  lineheight = 1.4, margin = margin(t = 12))
    )
  )

# ─── 9. Save ─────────────────────────────────────────────────────────────────
png_path <- file.path(out_dir, "Fig_synthesis.png")
pdf_path <- file.path(out_dir, "Fig_synthesis.pdf")

ggsave(png_path, fig, width = 14, height = 10, dpi = 300, bg = "white")

cairo_pdf(pdf_path, width = 14, height = 10)
print(fig)
dev.off()

message("Saved: ", png_path)
message("Saved: ", pdf_path)
message("09c_synthesis_figure.R complete.")
