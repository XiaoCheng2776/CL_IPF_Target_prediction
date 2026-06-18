#!/usr/bin/env Rscript
# 09d_crossspecies_scatter.R
# Four-quadrant cross-species scatter: mouse vs human NicheNet evidence
#
# Evidence metric: -log10(best_ligand_rank / N)   [higher = better]
#   Mouse: N = 154 ligands (NicheNet v2 mouse expressed-ligand pool, Strunz)
#   Human: N = 361 ligands (Habermann NicheNet; source: "top 7% of 361 ligands" §3)
# Hit threshold: rank ≤ 10  →  mouse score ≥ 1.19, human score ≥ 1.56
#
# Point types:
#   Filled circle  — scored in both networks (main quadrant)
#   Open triangle  — mouse rank only; plotted in bottom strip (y < 0, "not assessable in human")
#   Open diamond   — human rank only; plotted in left strip  (x < 0, "not assessable in mouse")
#   No NicheNet rank at all → excluded from scatter (CellChat/Geneformer-only candidates)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
  library(ggrepel)
})

COL_AGONIST <- "#E07B39"   # warm orange — agonist / agonist_uncertain
COL_BLOCK   <- "#4472B8"   # steel blue  — block
COL_POSCTRL <- "#B83232"   # brick red   — Fgfr2 positive control
COL_UNCERT  <- "#AAAAAA"   # grey        — uncertain direction

out_dir <- "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ─── 1. Parameters ───────────────────────────────────────────────────────────
N_MOUSE  <- 154     # mouse ligand pool (NicheNet v2 expressed ligands, Strunz)
N_HUMAN  <- 361     # human ligand pool (Habermann; stated in NicheNet log §3)
HIT_RANK <- 10      # hit threshold (rank ≤ 10)

X_THRESH <- -log10(HIT_RANK / N_MOUSE)   # 1.188
Y_THRESH <- -log10(HIT_RANK / N_HUMAN)   # 1.558

STRIP_Y  <- -0.55   # y-position for mouse-only strip (points with no human rank)
STRIP_X  <- -0.55   # x-position for human-only strip (points with no mouse rank)

# ─── 2. Load + annotate ──────────────────────────────────────────────────────
df <- read.csv("results/candidates_ranked.csv", stringsAsFactors = FALSE)

df <- df |>
  mutate(
    x_score = ifelse(!is.na(nn_best_ligand_rank),
                     -log10(nn_best_ligand_rank / N_MOUSE), NA_real_),
    y_score = ifelse(!is.na(nn_human_best_rank),
                     -log10(nn_human_best_rank / N_HUMAN), NA_real_),
    col = case_when(
      tolower(receptor) == "fgfr2"               ~ COL_POSCTRL,
      grepl("agonist", therapeutic_direction)    ~ COL_AGONIST,
      grepl("block",   therapeutic_direction)    ~ COL_BLOCK,
      TRUE                                       ~ COL_UNCERT
    ),
    is_pc = tolower(receptor) == "fgfr2",
    grp   = case_when(
      !is.na(x_score) & !is.na(y_score) ~ "both",
      !is.na(x_score) & is.na(y_score)  ~ "mouse_only",
      is.na(x_score)  & !is.na(y_score) ~ "human_only",
      TRUE                               ~ "neither"
    )
  )

d_both  <- filter(df, grp == "both")
d_mo    <- filter(df, grp == "mouse_only")   # no human rank → strip
d_ho    <- filter(df, grp == "human_only")   # no mouse rank → strip

cat(sprintf("Mouse N=%d  |  Human N=%d\n", N_MOUSE, N_HUMAN))
cat(sprintf("Mouse hit threshold: rank ≤ %d  -> score ≥ %.3f\n", HIT_RANK, X_THRESH))
cat(sprintf("Human hit threshold: rank ≤ %d  -> score ≥ %.3f\n", HIT_RANK, Y_THRESH))
cat("\nBoth scored:", paste(sort(d_both$receptor), collapse=", "), "\n")
cat("Mouse-only rank:", paste(sort(d_mo$receptor), collapse=", "), "\n")
cat("Human-only rank:", paste(sort(d_ho$receptor), collapse=", "), "\n")
cat("No NicheNet rank (excluded):", paste(sort(filter(df, grp=="neither")$receptor), collapse=", "), "\n")

# ─── 3. Assign plot coordinates ──────────────────────────────────────────────
d_both$x_plot <- d_both$x_score
d_both$y_plot <- d_both$y_score

d_mo$x_plot <- d_mo$x_score
d_mo$y_plot <- STRIP_Y   # placed in bottom "not assessable in human" strip

d_ho$x_plot <- STRIP_X   # placed in left "not assessable in mouse" strip
d_ho$y_plot <- d_ho$y_score

plot_df <- bind_rows(d_both, d_mo, d_ho)

# Quadrant membership for both-scored receptors
d_both <- d_both |>
  mutate(quadrant = case_when(
    x_score >= X_THRESH & y_score >= Y_THRESH ~ "conserved",
    x_score <  X_THRESH & y_score >= Y_THRESH ~ "human_hit",
    x_score >= X_THRESH & y_score <  Y_THRESH ~ "mouse_hit",
    TRUE                                       ~ "nonhit"
  ))

cat("\nConserved quadrant:", paste(d_both$receptor[d_both$quadrant=="conserved"], collapse=", "), "\n")
cat("Human-only hit (both scored):", paste(d_both$receptor[d_both$quadrant=="human_hit"], collapse=", "), "\n")
cat("Mouse-only hit (both scored):", paste(d_both$receptor[d_both$quadrant=="mouse_hit"], collapse=", "), "\n")

# ─── 4. Label logic ──────────────────────────────────────────────────────────
# Label: (a) all conserved hits, (b) Fgfr2, (c) user's key candidates,
#         (d) strip hits (above threshold on their respective axis), (e) others smaller
key_targets <- tolower(c("Sdc1","Sdc4","Erbb2","Egfr","Cd9","Fgfr2",
                          "App","Tgfbr2","Ramp1","Nrp1","Sort1","Cd63",
                          "Osmr","Il6st","Erbb3","Bmpr2"))

plot_df <- plot_df |>
  mutate(
    in_conserved = grp == "both" & x_plot >= X_THRESH & y_plot >= Y_THRESH,
    label_text   = receptor,
    label_size   = ifelse(in_conserved | is_pc, 4.3, 3.4),
    label_face   = ifelse(in_conserved | is_pc, "bold", "plain"),
    label_show   = in_conserved | is_pc |
                   tolower(receptor) %in% key_targets |
                   (grp == "mouse_only" & x_score >= X_THRESH) |
                   (grp == "human_only" & y_score >= Y_THRESH),
    label_text   = ifelse(label_show, label_text, ""),
    # Strip hit indicator (above threshold in respective axis)
    is_strip_hit = (grp == "mouse_only" & x_score >= X_THRESH) |
                   (grp == "human_only" & y_score >= Y_THRESH)
  )

# ─── 5. Plot geometry ────────────────────────────────────────────────────────
x_max <- max(plot_df$x_plot, na.rm=TRUE) + 0.30
y_max <- max(plot_df$y_plot, na.rm=TRUE) + 0.30
x_min <- STRIP_X - 0.45
y_min <- STRIP_Y - 0.45

# Quadrant label positions (main area)
qx_r  <- (X_THRESH + x_max) / 2       # right half
qx_l  <- (0 + X_THRESH) / 2           # left half (main area)
qy_t  <- (Y_THRESH + y_max) / 2       # top half
qy_b  <- (0 + Y_THRESH) / 2           # bottom half (main area)

rule_text <- sprintf(
  "Hit rule: rank ≤ %d\nMouse: score ≥ %.2f  (top %d of %d ligands, ~%.0f%%)\nHuman: score ≥ %.2f  (top %d of %d ligands, ~%.0f%%)",
  HIT_RANK,
  X_THRESH, HIT_RANK, N_MOUSE, 100*HIT_RANK/N_MOUSE,
  Y_THRESH, HIT_RANK, N_HUMAN, 100*HIT_RANK/N_HUMAN
)

p <- ggplot(plot_df, aes(x = x_plot, y = y_plot)) +

  # ── Not-assessable strips (grey backgrounds) ──────────────────────────────
  annotate("rect", xmin=x_min, xmax=0, ymin=y_min, ymax=y_max+0.2,
           fill="#EDEDED", color=NA) +
  annotate("rect", xmin=0, xmax=x_max+0.2, ymin=y_min, ymax=0,
           fill="#EDEDED", color=NA) +
  # Strip boundary lines
  annotate("segment", x=0, xend=0, y=0, yend=y_max+0.2,
           color="#BBBBBB", linewidth=0.5) +
  annotate("segment", x=0, xend=x_max+0.2, y=0, yend=0,
           color="#BBBBBB", linewidth=0.5) +
  # Strip labels
  annotate("text", x=(x_min+0)/2, y=(Y_THRESH+y_max)/2,
           label="No mouse\nNicheNet rank\n(not assessable)",
           color="#999999", size=3.2, fontface="italic", angle=90, hjust=0.5,
           lineheight=1.2) +
  annotate("text", x=(X_THRESH+x_max)/2, y=(y_min+0)/2,
           label="No human NicheNet rank  (not assessable in human)",
           color="#999999", size=3.2, fontface="italic", hjust=0.5,
           lineheight=1.2) +

  # ── Threshold lines ───────────────────────────────────────────────────────
  geom_vline(xintercept=X_THRESH, linetype="dashed",
             color="#444444", linewidth=0.65) +
  geom_hline(yintercept=Y_THRESH, linetype="dashed",
             color="#444444", linewidth=0.65) +

  # ── Quadrant labels (main area) ───────────────────────────────────────────
  annotate("text", x=qx_r, y=y_max+0.1,
           label="Conserved hit",
           color="#222222", fontface="bold", size=5.0, hjust=0.5) +
  annotate("text", x=qx_l, y=y_max+0.1,
           label="Human-only hit",
           color="#555555", fontface="bold", size=4.2, hjust=0.5) +
  annotate("text", x=qx_r, y=qy_b - 0.05,
           label="Mouse-only hit\n(both scored)",
           color="#555555", fontface="bold", size=3.8, hjust=0.5,
           lineheight=1.2) +
  annotate("text", x=qx_l, y=qy_b - 0.05,
           label="Non-hit",
           color="#BBBBBB", fontface="italic", size=3.8, hjust=0.5) +

  # ── Hit-rule box (bottom-right of main area) ──────────────────────────────
  annotate("text", x=x_max-0.03, y=0.04,
           label=rule_text, color="#555555", size=2.85,
           hjust=1, vjust=0, lineheight=1.3, fontface="plain") +

  # ── Points ────────────────────────────────────────────────────────────────
  # Both-scored: filled circles
  geom_point(data=filter(plot_df, grp=="both"),
             aes(color=col), size=4.5, shape=16) +
  # Mouse-only rank: open triangle (pointing up) in bottom strip
  geom_point(data=filter(plot_df, grp=="mouse_only"),
             aes(color=col), size=4.2, shape=24, fill="white", stroke=1.4) +
  # Human-only rank: open diamond in left strip
  geom_point(data=filter(plot_df, grp=="human_only"),
             aes(color=col), size=4.2, shape=23, fill="white", stroke=1.4) +
  # Positive-control star overlay on Fgfr2
  geom_point(data=filter(plot_df, is_pc),
             aes(color=col), size=9, shape=8, stroke=1.6) +

  # ── Labels (ggrepel) ─────────────────────────────────────────────────────
  geom_text_repel(
    aes(label=label_text, color=col, fontface=label_face, size=label_size),
    box.padding=0.5, point.padding=0.3,
    segment.color="#BBBBBB", segment.size=0.4, segment.alpha=0.8,
    max.overlaps=40, seed=42, na.rm=TRUE,
    min.segment.length=0.2
  ) +
  scale_size_identity() +

  # ── Colour and axis scales ─────────────────────────────────────────────────
  scale_color_identity() +
  scale_x_continuous(
    name   = paste0("–log₁₀(mouse NicheNet rank / ", N_MOUSE, ")"),
    limits = c(x_min, x_max+0.2),
    breaks = seq(0, floor(x_max+0.2), 0.5)
  ) +
  scale_y_continuous(
    name   = paste0("–log₁₀(human NicheNet rank / ", N_HUMAN, ")"),
    limits = c(y_min, y_max+0.25),
    breaks = seq(0, floor(y_max+0.25), 0.5)
  ) +

  # ── Theme ─────────────────────────────────────────────────────────────────
  theme_classic(base_size=14) +
  theme(
    axis.text    = element_text(color="black", size=12),
    axis.title   = element_text(color="black", size=13),
    plot.title   = element_text(size=15, face="bold", color="black",
                                margin=margin(b=8)),
    plot.caption = element_text(size=8.5, color="grey40", hjust=0,
                                lineheight=1.4, margin=margin(t=10)),
    plot.margin  = margin(t=20, r=20, b=10, l=15)
  ) +
  labs(
    title   = "Cross-Species NicheNet Evidence per Receptor",
    caption = paste0(
      "Filled circles: scored in both networks. Open △ = mouse rank only; open ◇ = human rank only ",
      "(not assessable in the other network; plotted in grey strips).\n",
      "Receptors with no NicheNet rank in either species are excluded (CellChat/Geneformer-only candidates: ",
      paste(sort(filter(df, grp=="neither")$receptor), collapse=", "), ").\n",
      "Color: orange = agonist; blue = block TGFβ; grey = uncertain direction; red ★ = positive control (FGFR2b).\n",
      "Itgb1: uncertain direction (both ANGPTL4 agonist-plausible and TGFB1 pro-fibrotic ligands). ",
      "Sdc1/Sdc4 overlap at identical coordinates (both: mouse rank 8, human rank 2)."
    )
  )

# ─── 6. Save ─────────────────────────────────────────────────────────────────
png_path <- file.path(out_dir, "Fig_crossspecies_scatter.png")
pdf_path <- file.path(out_dir, "Fig_crossspecies_scatter.pdf")

ggsave(png_path, p, width=10, height=9, dpi=300, bg="white")

cairo_pdf(pdf_path, width=10, height=9)
print(p)
dev.off()

message("Saved: ", png_path)
message("Saved: ", pdf_path)
message("09d_crossspecies_scatter.R complete.")
