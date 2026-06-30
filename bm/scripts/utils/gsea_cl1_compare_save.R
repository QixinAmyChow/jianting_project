#!/usr/bin/env Rscript
# gsea_cl1_compare_save.R
# Persist the exact GSEA selections used by the current slide 12 and slide 14
# figures (after Option-4 reconciliation) to RDS + CSV.

suppressPackageStartupMessages({
  library(dplyr)
})

source("gsea_families_shared.R")

OUT_DIR     <- "bm/figures/update3"
CANCER_CSV  <- file.path(OUT_DIR, "gsea_cancer_results.csv")
BROADER_CSV <- file.path(OUT_DIR, "gsea_broader_results.csv")

cancer <- read.csv(CANCER_CSV)
cancer$cluster  <- as.character(cancer$cluster)
cancer$database <- ifelse(grepl("^HALLMARK_", cancer$pathway), "HALLMARK", "C6")
broader <- read.csv(BROADER_CSV)
broader$cluster <- as.character(broader$cluster)
gsea_all <- bind_rows(cancer, broader) %>% filter(!is.na(NES))

families <- shared_families

sem_dedup <- function(df, fams) {
  assign_fam <- function(pw) {
    for (nm in names(fams)) {
      if (grepl(fams[[nm]], pw, ignore.case = TRUE, perl = TRUE)) return(nm)
    }
    paste0("UNIQUE__", pw)
  }
  df$sem_family <- vapply(df$pathway, assign_fam, character(1))
  df %>% group_by(sem_family) %>%
    slice_min(padj, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(-sem_family)
}

# ── Slide 12 (05e): all clusters, top 5 up + top 5 down per cluster ─────────
clusters <- as.character(sort(as.integer(unique(gsea_all$cluster))))
slide12 <- bind_rows(lapply(clusters, function(cl) {
  d <- gsea_all %>% filter(cluster == cl)
  up_pool <- d %>% filter(NES > 0) %>% slice_max(abs(NES), n = 20, with_ties = FALSE)
  dn_pool <- d %>% filter(NES < 0) %>% slice_max(abs(NES), n = 20, with_ties = FALSE)
  up5 <- sem_dedup(up_pool, families) %>% slice_min(padj, n = 5, with_ties = FALSE)
  dn5 <- sem_dedup(dn_pool, families) %>% slice_min(padj, n = 5, with_ties = FALSE)
  bind_rows(dn5, up5) %>% mutate(cluster = cl)
})) %>% mutate(direction = ifelse(NES > 0, "Up", "Down"))

# ── Slide 14 (05f): cluster 1, top 5 per hypothesis group ───────────────────
cl1 <- gsea_all %>% filter(cluster == "1")
hypo_keywords <- list(
  "ECM remodeling"             = "COLLAGEN|EXTRACELLULAR_MATRIX|ECM|MMP|MATRIX_METAL|FIBRONECT|LAMININ|INTEGRIN|ECM_REMODEL|BASEMENT_MEMBRANE",
  "TGF-b / EMT"                = "TGF.?B|TRANSFORMING_GROWTH|EPITHELIAL_MESENCH|EMT|SMAD|TGFB",
  "Wnt / osteogenic niche"     = "WNT|OSTEOBLAST|OSTEOCLAST|BONE_REMODEL|NOTCH|HEDGEHOG|BMP|OSSIF",
  "Immunosuppression / T cell" = "T_CELL|TCELL|CD8|CD4|TREG|EXHAUST|CHECKPOINT|PD.?1|PDL|CTLA|LAG3|TIM3|IMMUNE_EVASI|IMMUNO_SUPPRESS|IL10|MDSC",
  "Macrophage polarization"    = "MACRO|MACROPHAGE|\\bM1\\b|\\bM2\\b|\\bTAM\\b|POLARIZ|MYELOID|MONOCYTE|ALTERNATIVE_ACTIV|CLASSIC_ACTIV|PPARG|WBP7",
  "Cytokine / NFkB / JAK-STAT" = "CYTOKINE|NFKB|NF.KB|JAK.?STAT|STAT3|IL6|IL4|IL13|INTERFERON|CHEMOKINE|TNF"
)
slide14 <- bind_rows(lapply(names(hypo_keywords), function(h) {
  d <- cl1 %>% filter(grepl(hypo_keywords[[h]], pathway, ignore.case = TRUE, perl = TRUE))
  if (nrow(d) == 0) return(NULL)
  pool <- d %>% slice_max(abs(NES), n = 30, with_ties = FALSE)
  sem_dedup(pool, families) %>%
    slice_min(padj, n = 5, with_ties = FALSE) %>%
    mutate(hypothesis = h)
})) %>% mutate(direction = ifelse(NES > 0, "Up", "Down"))

# ── Save ────────────────────────────────────────────────────────────────────
rds_path <- file.path(OUT_DIR, "gsea_slide12_14_selections.rds")
csv12    <- file.path(OUT_DIR, "gsea_slide12_selection.csv")
csv14    <- file.path(OUT_DIR, "gsea_slide14_cl1_selection.csv")

saveRDS(
  list(
    slide12       = slide12,
    slide14       = slide14,
    hypo_keywords = hypo_keywords,
    families      = families
  ),
  rds_path
)
write.csv(slide12, csv12, row.names = FALSE)
write.csv(slide14, csv14, row.names = FALSE)

cat("Saved:", rds_path, "\n")
cat("Saved:", csv12, "\n")
cat("Saved:", csv14, "\n")
cat("  slide12 rows:", nrow(slide12), "| slide14 rows:", nrow(slide14), "\n")
