#!/usr/bin/env Rscript
# gsea_cl1_combined.R
# Cluster 1 broader GSEA — one diverging bar, no facets
#   fill:   database (GOBP=green, REACTOME=orange, IMMUNE=purple)
#   right-side labels: hypothesis group in grey (no hypothesis coloring)
#   dedup:  GLOBAL before hypothesis assignment — prevents same/parental terms
#           appearing from multiple databases or hypothesis groups

library(ggplot2)
library(dplyr)

source("gsea_families_shared.R")

OUT_DIR     <- "bm_analysis_out/figures_update_3"
CANCER_CSV  <- file.path(OUT_DIR, "gsea_cancer_results.csv")
BROADER_CSV <- file.path(OUT_DIR, "gsea_broader_results.csv")

# ── 1. Load all 5 databases (HALLMARK + C6 + GOBP + REACTOME + IMMUNE) ───────
cat("Loading cached GSEA results ...\n")
cancer <- read.csv(CANCER_CSV)
cancer$cluster  <- as.character(cancer$cluster)
cancer$database <- ifelse(grepl("^HALLMARK_", cancer$pathway), "HALLMARK", "C6")

broader <- read.csv(BROADER_CSV)
broader$cluster <- as.character(broader$cluster)

gsea <- bind_rows(cancer, broader) %>% filter(!is.na(NES))
cat("  H+C6 rows:", nrow(cancer), " | Broader rows:", nrow(broader),
    " | Combined:", nrow(gsea), "\n")

cl1 <- gsea %>% filter(cluster == "1")
cat("  Cluster 1 total rows:", nrow(cl1), "\n")

# ── 2. Semantic families — aggressive within-concept collapse ─────────────────
# One family per biological concept; first match wins.
# Key principle: all variant/parental/child terms of the same concept → same family.
families <- shared_families

sem_dedup <- function(df, fams) {
  assign_fam <- function(pw) {
    for (nm in names(fams)) {
      if (grepl(fams[[nm]], pw, ignore.case = TRUE, perl = TRUE)) return(nm)
    }
    paste0("UNIQUE__", pw)
  }
  df$sem_family <- vapply(df$pathway, assign_fam, character(1))
  df %>%
    group_by(sem_family) %>%
    slice_min(padj, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-sem_family)
}

# ── 3. Per-hypothesis dedup + top 5 by padj ──────────────────────────────────
# For each hypothesis group:
#   a) Filter cluster-1 terms matching that group's keywords
#   b) Pool top 30 by |NES| (up and down combined)
#   c) Semantic dedup → one representative per concept family
#   d) Take top 5 by padj

hypo_keywords <- list(
  "ECM remodeling"             = "COLLAGEN|EXTRACELLULAR_MATRIX|ECM|MMP|MATRIX_METAL|FIBRONECT|LAMININ|INTEGRIN|ECM_REMODEL|BASEMENT_MEMBRANE",
  "TGF-b / EMT"                = "TGF.?B|TRANSFORMING_GROWTH|EPITHELIAL_MESENCH|EMT|SMAD|TGFB",
  "Wnt / osteogenic niche"     = "WNT|OSTEOBLAST|OSTEOCLAST|BONE_REMODEL|NOTCH|HEDGEHOG|BMP|OSSIF",
  "Immunosuppression / T cell" = "T_CELL|TCELL|CD8|CD4|TREG|EXHAUST|CHECKPOINT|PD.?1|PDL|CTLA|LAG3|TIM3|IMMUNE_EVASI|IMMUNO_SUPPRESS|IL10|MDSC",
  "Macrophage polarization"    = "MACRO|MACROPHAGE|\\bM1\\b|\\bM2\\b|\\bTAM\\b|POLARIZ|MYELOID|MONOCYTE|ALTERNATIVE_ACTIV|CLASSIC_ACTIV|PPARG|WBP7",
  "Cytokine / NFkB / JAK-STAT" = "CYTOKINE|NFKB|NF.KB|JAK.?STAT|STAT3|IL6|IL4|IL13|INTERFERON|CHEMOKINE|TNF"
)
hypo_order <- names(hypo_keywords)

selected <- bind_rows(lapply(hypo_order, function(h) {
  d <- cl1 %>%
    filter(grepl(hypo_keywords[[h]], pathway, ignore.case = TRUE, perl = TRUE))
  if (nrow(d) == 0) return(NULL)
  pool <- d %>% slice_max(abs(NES), n = 30, with_ties = FALSE)
  dd   <- sem_dedup(pool, families)
  dd   %>% slice_min(padj, n = 5, with_ties = FALSE) %>%
    mutate(hypothesis = h)
}))

cat("  Terms per hypothesis after dedup:\n")
print(selected %>% count(hypothesis), row.names = FALSE)

cat("\n  Selected terms per hypothesis:\n")
print(selected %>% count(hypothesis), row.names = FALSE)

# ── 6. Short labels ───────────────────────────────────────────────────────────
shorten <- function(x, n = 50) {
  x <- gsub("^HALLMARK_",  "H: ", x)
  x <- gsub("^GOBP_",      "",    x)
  x <- gsub("^REACTOME_",  "R: ", x)
  x <- gsub("^GSE[0-9]+_", "",    x)
  x <- gsub("_", " ", x)
  x <- gsub("  +", " ", trimws(x))
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

selected <- selected %>%
  arrange(hypothesis, NES) %>%
  mutate(
    path_short = shorten(pathway),
    path_short = factor(path_short, levels = unique(path_short)),
    sig_label  = ifelse(padj < 0.05, "*", ifelse(padj < 0.1, ".", ""))
  )

cat("\nFinal terms:\n")
print(selected %>% select(hypothesis, path_short, NES, padj, database), row.names = FALSE)

# ── 7. Hypothesis group boundary info for annotations ────────────────────────
group_info <- selected %>%
  group_by(hypothesis) %>%
  summarise(
    y_min = min(as.integer(path_short)),
    y_max = max(as.integer(path_short)),
    .groups = "drop"
  ) %>%
  mutate(y_mid = (y_min + y_max) / 2)

n_terms <- nlevels(selected$path_short)
sep_y   <- (group_info$y_max[-nrow(group_info)] + group_info$y_min[-1]) / 2

# ── 8. Plot ───────────────────────────────────────────────────────────────────
DB_COLORS <- c(
  "HALLMARK" = "#1f77b4",   # blue
  "C6"       = "#aec7e8",   # light blue
  "GOBP"     = "#2ca02c",   # green
  "REACTOME" = "#ff7f0e",   # orange
  "IMMUNE"   = "#9467bd"    # purple
)

x_range <- max(abs(selected$NES), na.rm = TRUE) * 1.12

fig <- ggplot(selected, aes(x = NES, y = path_short, fill = database)) +
  geom_bar(stat = "identity", width = 0.72, color = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "grey30") +
  geom_text(aes(label = sig_label,
                x = NES + ifelse(NES >= 0, 0.05, -0.05)),
            hjust = ifelse(selected$NES >= 0, 0, 1),
            size = 3.5, color = "grey20") +
  scale_fill_manual(values = DB_COLORS, name = "Database") +
  scale_x_continuous(
    limits = c(-x_range, x_range * 1.15),
    expand = expansion(0)
  ) +
  labs(x = "NES", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y     = element_text(size = 9),
    axis.text.x     = element_text(size = 9),
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 8.5, color = "grey40"),
    legend.position = "bottom",
    legend.title    = element_text(face = "bold", size = 10),
    legend.text     = element_text(face = "bold", size = 9),
    plot.margin     = margin(10, 5, 10, 5)
  )

fig_h <- max(8, 0.3 * n_terms + 3.5)
out_pdf <- file.path(OUT_DIR, "05f_cl1_combined_nofacet.pdf")
out_png <- file.path(OUT_DIR, "05f_cl1_combined_nofacet.png")
ggsave(out_pdf, fig, width = 14, height = fig_h, limitsize = FALSE)
ggsave(out_png, fig, width = 14, height = fig_h, dpi = 150, limitsize = FALSE)
cat("\nSaved:", out_png, " (", fig_h, "in tall)\n")
cat("===== gsea_cl1_combined.R complete =====\n")
