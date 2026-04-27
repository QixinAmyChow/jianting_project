#!/usr/bin/env Rscript
# fix_dedup_gsea_figs.R
# Regenerates 05c (slide 14) with:
#   - Semantic deduplication (keep smallest padj within a pathway family)
#   - Top-5 selection by smallest padj (not by |NES|)
# Reads from cached CSVs — does NOT re-run GSEA

library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(cowplot)

OUT_DIR     <- "bm_analysis_out/figures_update_3"
GSEA_CSV    <- file.path(OUT_DIR, "gsea_cancer_results.csv")
BROADER_CSV <- file.path(OUT_DIR, "gsea_broader_results.csv")

BOLD_LEG <- theme(
  legend.title = element_text(face = "bold", size = 11),
  legend.text  = element_text(face = "bold", size = 10)
)

# ─── Semantic dedup helper ───────────────────────────────────────────────────
# families: named list of regex patterns (checked in order; first match wins)
# Within each matched family, keep the row with smallest padj.
# Unmatched pathways each form their own singleton family → never deduplicated.
sem_dedup <- function(df, families) {
  assign_fam <- function(pw) {
    for (nm in names(families)) {
      if (grepl(families[[nm]], pw, ignore.case = TRUE, perl = TRUE))
        return(nm)
    }
    return(paste0("UNIQUE__", pw))   # singleton → never merged
  }
  df$sem_family <- vapply(df$pathway, assign_fam, character(1))
  df <- df %>%
    group_by(sem_family) %>%
    slice_min(padj, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-sem_family)
  df
}


# ═══════════════════════════════════════════════════════════════════════════════
# PART 2 — Slide 14: 05c_cluster1_db_gsea + 05c_cluster1_hypo_gsea (broader)
# ═══════════════════════════════════════════════════════════════════════════════
cat("=== SLIDE 14: 05c_cluster1_db_gsea + 05c_cluster1_hypo_gsea ===\n")
gsea_broader <- read.csv(BROADER_CSV)
gsea_broader$cluster <- as.character(gsea_broader$cluster)
cat("  Broader cache rows:", nrow(gsea_broader), "\n")

shorten_broad <- function(x, n = 55) {
  x <- gsub("^GOBP_|^REACTOME_|^GSE[0-9]+_", "", x)
  x <- gsub("_", " ", x)
  x <- gsub("  +", " ", trimws(x))
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

# Hypothesis-keyword filter (identical to gsea_broader_update3.R)
hypo_keywords <- list(
  "ECM remodeling"                       = "COLLAGEN|EXTRACELLULAR_MATRIX|ECM|MMP|MATRIX_METAL|FIBRONECT|LAMININ|INTEGRIN|ECM_REMODEL|BASEMENT_MEMBRANE",
  "TGF-b / EMT"                          = "TGF.?B|TRANSFORMING_GROWTH|EPITHELIAL_MESENCH|EMT|SMAD|TGFB",
  "Wnt / osteogenic niche"               = "WNT|OSTEOBLAST|OSTEOCLAST|BONE_REMODEL|NOTCH|HEDGEHOG|BMP",
  "Immunosuppression / T cell crosstalk" = "T_CELL|CD8|CD4|TREG|EXHAUST|CHECKPOINT|PD.?1|PDL|CTLA|LAG3|TIM3|IMMUNE_EVASI|IMMUNO_SUPPRESS|IL10|TGFB.*IMMUNE",
  "Macrophage polarization"              = "MACRO|M1|M2|TAM|POLARIZ|MYELOID|MONOCYTE|ALTERNATIVE_ACTIV|CLASSIC_ACTIV|INFLAMMATORY_MAC|ANTI.INFLAM",
  "Cytokine / NFkB / JAK-STAT"          = "CYTOKINE|NFKB|NF.KB|JAK.?STAT|STAT3|IL6|IL4|IL13|INTERFERON|CHEMOKINE|TNF"
)

gsea_hypo <- bind_rows(lapply(names(hypo_keywords), function(h) {
  kw <- hypo_keywords[[h]]
  gsea_broader %>%
    filter(grepl(kw, pathway, ignore.case = TRUE, perl = TRUE), !is.na(NES)) %>%
    mutate(hypothesis = h)
}))

cat("Pathways per hypothesis (before dedup):\n")
print(gsea_hypo %>% filter(cluster == "1") %>%
      group_by(hypothesis) %>% summarise(n = n_distinct(pathway)), row.names = FALSE)

# Semantic families for broader GSEA pathways
broad_families <- list(
  # ECM remodeling
  ECM_ORG       = "ECM_ORGANIZATION|EXTRACELLULAR_MATRIX_ORGAN|EXTRACELLULAR_MATRIX_ASSEMBLY",
  COLLAGEN_FORM = "COLLAGEN_FIBRIL|COLLAGEN_BIOSYN",
  COLLAGEN_GEN  = "COLLAGEN",
  MMP           = "MMP|MATRIX_METALLOPROTEIN",
  INTEGRIN      = "INTEGRIN",
  FIBRONECT     = "FIBRONECT",
  LAMININ       = "LAMININ",
  BASEMENT_MEM  = "BASEMENT_MEMBRANE",
  # TGF-b / EMT
  TGFB_SIG      = "TGFB_SIGNALING|TGF_BETA_SIGNAL|SIGNALING_BY_TGF|TGF_BETA_RECEPT",
  TGFB_GEN      = "TGFB|TGF.?B",
  EMT_TRANS     = "EPITHELIAL_MESENCH|EMT_TRANS",
  EMT_GEN       = "EMT",
  SMAD          = "SMAD",
  # Wnt / osteogenic
  WNT_CANON     = "CANONICAL_WNT|WNT_CANONICAL|BETA_CATENIN|WNT_LIGAND",
  WNT_GEN       = "WNT",
  OSTEOBLAST    = "OSTEOBLAST|OSTEOGENIC_DIFF",
  OSTEOCLAST    = "OSTEOCLAST|BONE_RESORPT",
  BONE_REMODEL  = "BONE_REMODEL|BONE_MINERAL|OSSIFICATION",
  BMP           = "BMP",
  NOTCH         = "NOTCH",
  HEDGEHOG      = "HEDGEHOG",
  # Immunosuppression
  CD8_EXHAUST   = "CD8.*EXHAUST|EXHAUST.*CD8|T_CELL_EXHAUST",
  CD8_GEN       = "CD8",
  TREG          = "TREG|REGULATORY_T_CELL|T_REG",
  CHECKPOINT    = "CHECKPOINT|PD.?1|PDL1|CTLA4|LAG3|TIM3",
  IL10          = "IL10|IL_10",
  T_CELL_ACTIV  = "T_CELL_ACTIV|T_CELL_STIMUL|T_CELL_PROLIF",
  T_CELL_GEN    = "T_CELL",
  # Macrophage polarization
  M2_POLARIZ    = "M2_POLARIZ|ALTERNATIVE_ACTIV|ANTI.INFLAM.*MAC",
  M1_POLARIZ    = "M1_POLARIZ|CLASSIC_ACTIV|INFLAM.*MAC",
  TAM           = "TAM|TUMOR.*MACRO|MACRO.*TUMOR",
  MYELOID_DIFF  = "MYELOID_DIFF|MYELOID_CELL_DIFF",
  MYELOID_GEN   = "MYELOID",
  MONOCYTE      = "MONOCYTE",
  # Cytokine / NFkB / JAK-STAT
  NFKB          = "NFKB|NF.KB",
  JAK_STAT      = "JAK.?STAT|STAT3_TARGET",
  STAT_GEN      = "STAT",
  IFN_GAMMA     = "INTERFERON_GAMMA|IFN_GAMMA",
  IFN_ALPHA     = "INTERFERON_ALPHA|IFN_ALPHA|TYPE_I_IFN",
  IFN_GEN       = "INTERFERON",
  CHEMOKINE     = "CHEMOKINE",
  TNF           = "TNF",
  IL6           = "IL6|IL_6",
  IL4_IL13      = "IL4|IL13|IL_4|IL_13",
  CYTOKINE_PROD = "CYTOKINE_PRODUCT|CYTOKINE_SECRET",
  CYTOKINE_GEN  = "CYTOKINE"
)

# Cluster 1: dedup per hypothesis, top 5 by smallest padj
cl1 <- gsea_hypo %>%
  filter(cluster == "1") %>%
  group_by(hypothesis) %>%
  group_modify(~ {
    pool    <- .x %>% slice_max(abs(NES), n = 15, with_ties = FALSE)
    pool_dd <- sem_dedup(pool, broad_families)
    pool_dd %>% slice_min(padj, n = 5, with_ties = FALSE)
  }) %>%
  ungroup() %>%
  arrange(hypothesis, NES) %>%
  mutate(
    path_short = shorten_broad(pathway),
    path_short = factor(path_short, levels = unique(path_short)),
    dir        = ifelse(NES > 0, "Up-reg", "Down-reg")
  )

cat("\nCluster 1 terms per hypothesis after dedup:\n")
print(cl1 %>% group_by(hypothesis) %>% summarise(n = n()), row.names = FALSE)

cat("\nDetailed terms:\n")
print(cl1 %>% select(hypothesis, path_short, NES, padj, database) %>%
        arrange(hypothesis, NES), row.names = FALSE)

# ── Fig: colored by database ─────────────────────────────────────────────────
db_colors <- c("GOBP" = "#2ca02c", "REACTOME" = "#ff7f0e", "IMMUNE" = "#9467bd")

fig_cl1_db <- ggplot(cl1, aes(x = NES, y = path_short, fill = database)) +
  geom_bar(stat = "identity", width = 0.75) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
  scale_fill_manual(values = db_colors, name = "Database") +
  facet_wrap(~hypothesis, scales = "free_y", ncol = 2) +
  labs(
    title    = "Cluster 1 (94.8% BC_ER+) — Hypothesis-Relevant GSEA Pathways",
    subtitle = "GO:BP + Reactome + C7-IMMUNESIGDB | top 5 by padj (semantically deduplicated)",
    x = "NES", y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text      = element_text(face = "bold", size = 9),
    axis.text.y     = element_text(size = 7.5),
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 9),
    legend.position = "bottom"
  ) + BOLD_LEG

ggsave(file.path(OUT_DIR, "05c_cluster1_db_gsea.pdf"),  fig_cl1_db, width = 16, height = 12)
ggsave(file.path(OUT_DIR, "05c_cluster1_db_gsea.png"),  fig_cl1_db, width = 16, height = 12, dpi = 150)
cat("  Saved: 05c_cluster1_db_gsea\n")

# ── Fig: colored by direction (slide 15 — 05c_cluster1_hypo_gsea) ───────────
fig_cl1_dir <- ggplot(cl1, aes(x = NES, y = path_short, fill = dir)) +
  geom_bar(stat = "identity", width = 0.75) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
  scale_fill_manual(values = c("Up-reg" = "#d62728", "Down-reg" = "#2171b5"), name = NULL) +
  facet_wrap(~hypothesis, scales = "free_y", ncol = 2) +
  labs(
    title    = "Cluster 1 (94.8% BC_ER+) — Hypothesis-Relevant GSEA Pathways",
    subtitle = "GO:BP + Reactome + C7-IMMUNESIGDB | top 5 by padj (semantically deduplicated) | red = up, blue = down",
    x = "NES", y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text      = element_text(face = "bold", size = 9),
    axis.text.y     = element_text(size = 7.5),
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 9),
    legend.position = "bottom"
  ) + BOLD_LEG

ggsave(file.path(OUT_DIR, "05c_cluster1_hypo_gsea.pdf"), fig_cl1_dir, width = 16, height = 12)
ggsave(file.path(OUT_DIR, "05c_cluster1_hypo_gsea.png"), fig_cl1_dir, width = 16, height = 12, dpi = 150)
cat("  Saved: 05c_cluster1_hypo_gsea\n")

cat("\n===== fix_dedup_gsea_figs.R complete =====\n")
cat("Updated: 05c_cluster1_db_gsea | 05c_cluster1_hypo_gsea\n")
