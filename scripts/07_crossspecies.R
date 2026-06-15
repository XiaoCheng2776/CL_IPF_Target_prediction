#!/usr/bin/env Rscript
# =============================================================================
# 07_crossspecies.R — cross-species validation and final candidate prioritisation
# Run from repo root:  Rscript scripts/07_crossspecies.R
#
# NOTE: Human arm (Adams GSE136831, Habermann GSE135893) is intentionally deferred.
#       Currently runs mouse-only integration and produces a ranked candidate table.
#       When human data is added, re-run phases 3–4 on human objects and intersect
#       the ranked lists via orthology mapping in section 4 below.
#
# Input:  results/nichenet/ligand_activities.csv
#         results/nichenet/top_receptors.csv
#         results/cellchat/incoming_to_AT2cells.csv
#         results/geneformer/oe/receptors_ranked_oe.csv
#         results/geneformer/del/receptors_ranked_del.csv
#         results/nichenet_human/ligand_activities_human.csv  (Habermann)
#         results/nichenet_human/top_receptors_human.csv
# Output: results/candidates_ranked.csv
#         results/figures/candidates_dotplot.png
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

ROOT    <- here::here()
OUT_FIG <- file.path(ROOT, "results/figures")
dir.create(OUT_FIG, showWarnings = FALSE)

# =============================================================================
# 1.  Load per-method results
# =============================================================================
nn_ligands   <- read.csv(file.path(ROOT, "results/nichenet/ligand_activities.csv"))
nn_receptors <- read.csv(file.path(ROOT, "results/nichenet/top_receptors.csv"))
cc_at2       <- read.csv(file.path(ROOT, "results/cellchat/incoming_to_AT2cells.csv"))
gf_oe        <- read.csv(file.path(ROOT, "results/geneformer/oe/receptors_ranked_oe.csv"))
gf_del       <- read.csv(file.path(ROOT, "results/geneformer/del/receptors_ranked_del.csv"))

# Human NicheNet (Habermann GSE135893)
human_NN_file  <- file.path(ROOT, "results/nichenet_human/ligand_activities_human.csv")
human_rec_file <- file.path(ROOT, "results/nichenet_human/top_receptors_human.csv")
has_human_nn   <- file.exists(human_NN_file) && file.exists(human_rec_file)
if (has_human_nn) {
  nn_human_lig <- read.csv(human_NN_file)
  nn_human_rec <- read.csv(human_rec_file)
  message("Human NicheNet loaded: ", nrow(nn_human_rec), " receptor pairs")
} else {
  message("Human NicheNet not found — run 08_human_nichenet.R first. Mouse-only mode.")
}

# =============================================================================
# 2.  Build per-method receptor summaries
# =============================================================================

# NicheNet: best ligand rank + AUPR for each receptor
nn_rec <- nn_receptors %>%
  left_join(nn_ligands %>% select(test_ligand, aupr_corrected, rank),
            by = c("ligand" = "test_ligand")) %>%
  group_by(receptor) %>%
  summarise(nn_best_ligand_rank = min(rank,           na.rm = TRUE),
            nn_best_aupr        = max(aupr_corrected,  na.rm = TRUE),
            nn_ligands          = paste(sort(ligand),  collapse = ";"),
            .groups = "drop")

# CellChat: total incoming probability to AT2 cells, by receptor
cc_rec <- cc_at2 %>%
  group_by(receptor) %>%
  summarise(cc_total_prob = sum(prob, na.rm = TRUE),
            cc_pathways   = paste(unique(pathway_name), collapse = ";"),
            .groups = "drop")

# Geneformer: OE gain + DEL loss (columns: symbol, ensembl_id, regen_gain/loss)
gf_rec <- gf_oe %>%
  select(symbol, regen_gain) %>%
  rename(gf_regen_gain = regen_gain)

gf_del_rec <- gf_del %>%
  select(symbol, regen_loss) %>%
  rename(gf_regen_loss = regen_loss)

# =============================================================================
# 3.  Union candidate list and join all scores
# =============================================================================
# Union of receptors seen in any method (capitalise for join key)
all_receptors <- union(
  union(nn_rec$receptor, cc_rec$receptor),
  union(gf_rec$symbol,   gf_del_rec$symbol)
)

candidates <- tibble(receptor = all_receptors) %>%
  left_join(nn_rec,                                 by = "receptor") %>%
  left_join(cc_rec,                                 by = "receptor") %>%
  left_join(gf_rec,      by = c("receptor" = "symbol")) %>%
  left_join(gf_del_rec,  by = c("receptor" = "symbol"))

# =============================================================================
# 4.  Cross-species conservation (Habermann human NicheNet)
# =============================================================================
if (has_human_nn) {
  human_rec_ranked <- nn_human_rec %>%
    left_join(nn_human_lig %>% select(test_ligand, aupr_corrected, rank),
              by = c("ligand" = "test_ligand")) %>%
    group_by(receptor) %>%
    summarise(nn_human_best_rank = min(rank,          na.rm = TRUE),
              nn_human_best_aupr = max(aupr_corrected, na.rm = TRUE),
              nn_human_ligands   = paste(sort(ligand), collapse = ";"),
              .groups = "drop")
  candidates <- candidates %>%
    mutate(receptor_upper = toupper(receptor)) %>%
    left_join(human_rec_ranked, by = c("receptor_upper" = "receptor"))
  candidates$species_support <- ifelse(
    !is.na(candidates$nn_human_best_rank), "conserved", "mouse_only"
  )
  message(sprintf("Conserved (mouse+human NicheNet): %d",
                  sum(candidates$species_support == "conserved")))
} else {
  candidates$species_support    <- "mouse_only"
  candidates$nn_human_best_rank <- NA_real_
  candidates$nn_human_best_aupr <- NA_real_
  candidates$nn_human_ligands   <- NA_character_
}

# =============================================================================
# 5.  Composite scoring
# =============================================================================
# rank_norm: maps any numeric vector to 0-1 where 1 = best.
# higher_is_better=TRUE  → largest value → 1
# higher_is_better=FALSE → smallest value → 1
rank_norm <- function(x, higher_is_better = TRUE) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  r <- rank(if (higher_is_better) x else -x, na.last = "keep", ties.method = "average")
  (r - 1) / (max(r, na.rm = TRUE) - 1)
}

candidates <- candidates %>%
  mutate(
    # NicheNet mouse: low ligand rank = better
    score_nn      = rank_norm(nn_best_ligand_rank,  higher_is_better = FALSE),
    # NicheNet human: low ligand rank = better (bonus dimension)
    score_nn_h    = rank_norm(nn_human_best_rank,   higher_is_better = FALSE),
    # CellChat: higher incoming prob = better
    score_cc      = rank_norm(cc_total_prob,         higher_is_better = TRUE),
    # Geneformer OE: higher regen_gain = better
    score_gf      = rank_norm(gf_regen_gain,          higher_is_better = TRUE),
    # Geneformer DEL: more negative regen_loss = deletion hurts regeneration = better target
    score_del     = rank_norm(gf_regen_loss,           higher_is_better = FALSE),
    # Composite: equal weight; NA ignored; conserved gets extra human NicheNet column
    composite_score = rowMeans(
      cbind(score_nn, score_nn_h, score_cc, score_gf, score_del), na.rm = TRUE)
  ) %>%
  arrange(desc(composite_score))

# =============================================================================
# 6.  Annotations
# =============================================================================
# 6b. Therapeutic direction annotation
# =============================================================================
# NicheNet/CellChat identify WHICH ligands are signaling to AT2 cells in fibrosis.
# They do NOT indicate therapeutic direction.  TGFB1 is the most active sender
# signal but it DRIVES the arrested/Krt8+ state — its receptors are block candidates
# (antagonize to let AT2 cells escape arrest), not agonist targets.
# This column separates pro-regenerative agonist candidates from pro-fibrotic
# block candidates for downstream prioritisation.

# TGFβ-pathway annotation: receptors whose canonical ligand is TGFB1
# and/or which transduce TGFβ signaling (SMAD2/3 → pEMT, Krt8 upregulation,
# senescence) in the context of bleomycin/IPF epithelium.
TGFB_RECEPTORS <- c(
  "Tgfbr1", "Tgfbr2",   # direct TGFβ serine-kinase receptors (SMAD2/3)
  "App",                  # TGFB1→APP: top ligand both species; promotes AT2 senescence
  "Itgav",               # αv integrins activate latent TGFβ (αvβ6 canonical activator)
  "Itgb6", "Itgb8"       # αvβ6 / αvβ8 — epithelial latent-TGFβ activators
)

# Receptors with a known pro-fibrotic role independent of TGFβ
PROFIBROTIC_OTHER <- c(
  "Axl"   # AXL promotes pEMT/fibrosis; negative control in Geneformer
)

candidates <- candidates %>%
  mutate(
    therapeutic_direction = case_when(
      # ----- block candidates: TGFβ pathway -----
      receptor %in% TGFB_RECEPTORS ~ "block_TGFb",
      receptor %in% PROFIBROTIC_OTHER ~ "block_other",

      # ----- agonist candidates: validated pro-regenerative axes -----
      # FGFR2b positive control
      receptor %in% c("Fgfr2", "Fgfr1") ~ "agonist",
      # Syndecans: heparan-sulfate co-receptors that present FGF10 to FGFR2b;
      #   also downstream of pro-survival ANGPTL4
      receptor %in% c("Sdc1", "Sdc4") ~ "agonist",
      # EGF-family: AREG/HBEGF→EGFR/ERBB2 promotes AT2 proliferation & AT2→AT1
      receptor %in% c("Egfr", "Erbb2", "Erbb3") ~ "agonist",
      # HBEGF/AREG co-receptors (tetraspanin web)
      receptor %in% c("Cd9", "Cd63") ~ "agonist",
      # WNT: Fzd2 = top Geneformer OE hit (+0.008); WNT promotes AT2 self-renewal
      receptor == "Fzd2" ~ "agonist",
      # FGFR2b co-receptor / VEGF axis
      receptor == "Nrp1" ~ "agonist",
      # BMP axis: BMP4 maintains AT2 identity in alveolar organoids
      receptor %in% c("Bmpr1a", "Bmpr2") ~ "agonist_uncertain",
      # Adm/CGRP receptor: vasculoprotective; direct AT2 role uncertain
      receptor == "Ramp1" ~ "agonist_uncertain",

      # ----- uncertain: Itgb1 is driven by BOTH ANGPTL4 and TGFB1 -----
      receptor == "Itgb1" ~ "uncertain",
      # Notch promotes basaloid fate in fibrotic lung; likely block but not confirmed
      receptor %in% c("Notch1", "Notch2") ~ "uncertain",

      TRUE ~ "uncertain"
    ),

    notes = case_when(
      receptor == "Fgfr2" ~
        "POSITIVE CONTROL — FGF10→FGFR2b; IIIb (epithelial) vs IIIc (mesenchymal) isoform ambiguity",
      receptor %in% c("Tgfbr1","Tgfbr2") ~
        "TGFβ receptor (SMAD2/3) — drives Krt8+/arrested/fibrotic state; BLOCK candidate",
      receptor == "App" ~
        "Top ligand = TGFB1 (mouse rank 1, human rank 5); TGFβ→APP promotes AT2 senescence; BLOCK candidate",
      receptor == "Itgb1" ~
        "Mixed: ANGPTL4→ITGB1 (pro-survival) AND TGFB1→ITGB1 (pro-fibrotic); resolve with isoform/context data",
      receptor %in% c("Sdc1","Sdc4") ~
        "Heparan-sulfate co-receptor; presents FGF10 to FGFR2b (ANGPTL4/LPL ligands); agonist candidate",
      receptor %in% c("Egfr","Erbb2","Erbb3") ~
        "AREG/HBEGF→EGFR family; promotes AT2 proliferation and AT2→AT1 transition; agonist candidate",
      receptor %in% c("Cd9","Cd63") ~
        "Tetraspanin co-receptor for HBEGF/TIMP1; scaffolds EGFR signaling; agonist candidate",
      receptor == "Fzd2" ~
        "WNT receptor; highest Geneformer OE regen_gain (+0.008); WNT promotes AT2 self-renewal",
      receptor == "Nrp1" ~
        "FGFR2b co-receptor and VEGFA-axis; mouse NicheNet rank 8; agonist candidate",
      receptor %in% c("Bmpr1a","Bmpr2") ~
        "BMP receptor; BMP4 maintains AT2 identity; direction uncertain in fibrotic context",
      receptor %in% c("Notch1","Notch2") ~
        "Notch may promote basaloid/arrested fate in IPF epithelium; uncertain/block",
      receptor == "Axl" ~
        "Pro-fibrotic (AXL promotes pEMT/mesenchymal transition); BLOCK candidate; GF negative control",
      receptor == "Ramp1" ~
        "ADM/CGRP receptor; vasculoprotective; direct AT2 regeneration role uncertain",
      TRUE ~ NA_character_
    )
  )

# =============================================================================
# 7.  Print and save
# =============================================================================
label <- if (has_human_nn) "mouse+human" else "mouse-only"

message(sprintf("\n=== Final candidate ranking (%s) ===", label))
print(candidates %>%
        select(receptor, composite_score, therapeutic_direction, species_support,
               nn_best_ligand_rank, nn_human_best_rank,
               gf_regen_gain, gf_regen_loss) %>%
        head(20), n = 20)

message("\n=== AGONIST candidates (pro-regenerative; agonism promotes AT2→AT1) ===")
agonists <- candidates %>%
  filter(grepl("^agonist", therapeutic_direction)) %>%
  select(receptor, composite_score, therapeutic_direction, species_support,
         nn_ligands, nn_human_ligands, gf_regen_gain, gf_regen_loss, notes)
print(agonists, n = 30)

message("\n=== BLOCK candidates (pro-fibrotic; antagonism may release arrest) ===")
blocks <- candidates %>%
  filter(grepl("^block", therapeutic_direction)) %>%
  select(receptor, composite_score, therapeutic_direction, species_support,
         nn_ligands, nn_human_ligands, notes)
print(blocks, n = 20)

message("\n=== UNCERTAIN direction ===")
print(candidates %>%
        filter(therapeutic_direction == "uncertain") %>%
        select(receptor, composite_score, species_support, nn_ligands, notes), n = 20)

write.csv(candidates, file.path(ROOT, "results/candidates_ranked.csv"),
          row.names = FALSE)

# Dot-plot: colour by therapeutic_direction; top 20 by composite score
dir_palette <- c(
  agonist          = "#2ca25f",
  agonist_uncertain = "#99d8c9",
  block_TGFb       = "#de2d26",
  block_other      = "#fc9272",
  uncertain        = "#bdbdbd"
)

top20 <- candidates %>%
  filter(!is.na(composite_score)) %>%
  slice_head(n = 20) %>%
  mutate(receptor = fct_reorder(receptor, composite_score))

p <- ggplot(top20, aes(x = composite_score, y = receptor)) +
  geom_point(aes(colour = therapeutic_direction,
                 shape  = species_support,
                 size   = ifelse(is.na(gf_regen_gain), 2,
                                 pmax(abs(gf_regen_gain) * 200 + 2, 2)))) +
  scale_colour_manual(values = dir_palette) +
  scale_shape_manual(values = c(mouse_only = 16, conserved = 18)) +
  scale_size_continuous(range = c(2, 8)) +
  labs(title    = "Top 20 candidate receptors — mouse IPF",
       subtitle = paste0("Green = agonist (promote AT2→AT1)  |  Red = block (antagonize to release arrest)\n",
                         "Shape: circle = mouse-only, diamond = cross-species conserved"),
       x = "Composite score (0–1)", y = NULL,
       colour = "Therapeutic direction", shape = "Species support") +
  theme_bw(base_size = 11) +
  guides(size = "none")

ggsave(file.path(OUT_FIG, "candidates_dotplot.png"),
       p, width = 10, height = 7, dpi = 300)

message("07_crossspecies.R complete.")
message("Output: results/candidates_ranked.csv")
message("Figure: results/figures/candidates_dotplot.png")
