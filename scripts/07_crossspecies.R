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
# FGFR2 isoform flag (see CLAUDE.md guardrail)
candidates <- candidates %>%
  mutate(notes = case_when(
    receptor == "Fgfr2"  ~ "POSITIVE CONTROL — IIIb (epithelial) vs IIIc (mesenchymal) isoform ambiguity; OE regen_gain negative (ligand-limited, not receptor-limited)",
    receptor %in% c("Tgfbr1","Tgfbr2") ~ "Expected negative: TGF-β promotes arrested/fibrotic state",
    receptor == "Axl"    ~ "Pro-fibrotic; expected low OE regen_gain (internal negative control)",
    TRUE                 ~ NA_character_
  ))

# =============================================================================
# 7.  Print and save
# =============================================================================
label <- if (has_human_nn) "mouse+human" else "mouse-only"
message(sprintf("\n=== Final candidate ranking (%s) ===", label))
print(candidates %>%
        select(receptor, composite_score, species_support,
               score_nn, score_nn_h, score_cc, score_gf, score_del,
               nn_best_ligand_rank, nn_human_best_rank,
               gf_regen_gain, gf_regen_loss) %>%
        head(20), n = 20)

message("\n=== Top 10 with full evidence ===")
print(candidates %>%
        select(receptor, composite_score, species_support,
               nn_ligands, nn_human_ligands, cc_pathways,
               gf_regen_gain, gf_regen_loss, notes) %>%
        head(10), n = 10)

if (has_human_nn) {
  message("\n=== Conserved candidates (mouse + human NicheNet) ===")
  print(candidates %>%
          filter(species_support == "conserved") %>%
          select(receptor, composite_score, nn_best_ligand_rank,
                 nn_human_best_rank, nn_human_ligands, gf_regen_gain), n = 30)
}

write.csv(candidates, file.path(ROOT, "results/candidates_ranked.csv"),
          row.names = FALSE)

# Dot-plot: top 20 with composite score
top20 <- candidates %>%
  filter(!is.na(composite_score)) %>%
  slice_head(n = 20) %>%
  mutate(receptor = fct_reorder(receptor, composite_score),
         gf_gain_clipped = pmax(gf_regen_gain, -0.02, na.rm = FALSE))

p <- ggplot(top20, aes(x = composite_score, y = receptor)) +
  geom_point(aes(colour = species_support,
                 size   = ifelse(is.na(gf_regen_gain), 1, pmax(gf_regen_gain + 0.025, 0.5)))) +
  scale_colour_manual(values = c(mouse_only = "#2b8cbe", conserved = "#f03b20")) +
  scale_size_continuous(range = c(2, 8)) +
  labs(title    = "Top 20 candidate receptors — mouse IPF",
       subtitle = "Composite score: NicheNet rank + CellChat prob + Geneformer OE gain + DEL consistency",
       x = "Composite score (0–1)", y = NULL,
       colour = "Species support", size = "GF regen gain (size)") +
  theme_bw(base_size = 11)

ggsave(file.path(OUT_FIG, "candidates_dotplot.png"),
       p, width = 9, height = 7, dpi = 300)

message("07_crossspecies.R complete.")
message("Output: results/candidates_ranked.csv")
message("Figure: results/figures/candidates_dotplot.png")
