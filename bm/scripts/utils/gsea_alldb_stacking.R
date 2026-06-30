#!/usr/bin/env Rscript
# gsea_alldb_stacking.R
# One stacking diverging-bar figure — ALL clusters × ALL 5 databases
# Databases: HALLMARK | C6 | GOBP | REACTOME | IMMUNE
# Deduplication: semantic families (same approach as fix_dedup_gsea_figs.R)
# Selection: top 5 up + top 5 down per cluster (smallest padj after dedup)
# Color: by database
# Reads only cached CSVs — no Seurat object needed

library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)

source("gsea_families_shared.R")

OUT_DIR     <- "bm/figures/update3"
CANCER_CSV  <- file.path(OUT_DIR, "gsea_cancer_results.csv")
BROADER_CSV <- file.path(OUT_DIR, "gsea_broader_results.csv")

# ── 1. Load & label databases ─────────────────────────────────────────────────
cat("Loading cached GSEA results ...\n")

cancer <- read.csv(CANCER_CSV)
cancer$cluster <- as.character(cancer$cluster)
cancer$database <- ifelse(grepl("^HALLMARK_", cancer$pathway), "HALLMARK", "C6")
cat("  H+C6 rows:", nrow(cancer), "\n")

broader <- read.csv(BROADER_CSV)
broader$cluster <- as.character(broader$cluster)
cat("  Broader rows:", nrow(broader), "\n")

gsea_all <- bind_rows(cancer, broader) %>% filter(!is.na(NES))
cat("  Combined rows:", nrow(gsea_all), "\n")
cat("  Databases:", paste(sort(unique(gsea_all$database)), collapse=" | "), "\n\n")

# ── 2. Semantic dedup families ────────────────────────────────────────────────
# First-match wins; singletons (no match) kept as-is
families <- shared_families

sem_dedup <- function(df, fams) {
  assign_fam <- function(pw) {
    for (nm in names(fams)) {
      if (grepl(fams[[nm]], pw, ignore.case = TRUE, perl = TRUE))
        return(nm)
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

# ── 3. Shorten pathway labels ─────────────────────────────────────────────────
shorten <- function(x, n = 42) {
  x <- gsub("^HALLMARK_", "H: ", x)
  x <- gsub("^GOBP_",     "",    x)
  x <- gsub("^REACTOME_", "R: ", x)
  x <- gsub("^GSE[0-9]+_","",    x)
  x <- gsub("_", " ", x)
  x <- gsub("  +", " ", trimws(x))
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

# ── 4. Per-cluster: dedup + top 5 up + top 5 down ────────────────────────────
clusters <- as.character(sort(as.integer(unique(gsea_all$cluster))))
cat("Clusters:", paste(clusters, collapse=" "), "\n")

top_df <- bind_rows(lapply(clusters, function(cl) {
  d <- gsea_all %>% filter(cluster == cl, !is.na(NES))

  up_pool <- d %>% filter(NES > 0) %>% slice_max(abs(NES), n = 20, with_ties = FALSE)
  dn_pool <- d %>% filter(NES < 0) %>% slice_max(abs(NES), n = 20, with_ties = FALSE)

  up_dd <- sem_dedup(up_pool, families)
  dn_dd <- sem_dedup(dn_pool, families)

  up5 <- up_dd %>% slice_min(padj, n = 5, with_ties = FALSE)
  dn5 <- dn_dd %>% slice_min(padj, n = 5, with_ties = FALSE)

  bind_rows(dn5, up5) %>% mutate(cluster = cl)
}))

top_df <- top_df %>%
  mutate(
    path_short  = shorten(pathway),
    direction   = ifelse(NES > 0, "Up", "Down"),
    cluster_lbl = paste0("Cluster ", cluster)
  )

cat("\nTerms per cluster:\n")
print(top_df %>% count(cluster, direction) %>% arrange(as.integer(cluster)), row.names=FALSE)

# ── 5. Database color palette ─────────────────────────────────────────────────
DB_COLORS <- c(
  "HALLMARK" = "#1f77b4",   # blue
  "C6"       = "#aec7e8",   # light blue
  "GOBP"     = "#2ca02c",   # green
  "REACTOME" = "#ff7f0e",   # orange
  "IMMUNE"   = "#9467bd"    # purple
)

# ── 6. Build per-cluster panels ───────────────────────────────────────────────
cl_ord <- as.character(sort(as.integer(unique(top_df$cluster))))

plot_list <- lapply(cl_ord, function(cl) {
  df <- top_df %>%
    filter(cluster == cl) %>%
    arrange(NES) %>%
    mutate(path_short = factor(path_short, levels = unique(path_short)))

  ggplot(df, aes(x = NES, y = path_short, fill = database)) +
    geom_bar(stat = "identity", width = 0.72) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
    scale_fill_manual(values = DB_COLORS, name = "Database",
                      drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(title = paste0("Cluster ", cl), x = "NES", y = NULL) +
    theme_classic(base_size = 7) +
    theme(
      axis.text.y     = element_text(size = 6),
      axis.text.x     = element_text(size = 6),
      axis.title.x    = element_text(size = 6.5),
      plot.title      = element_text(face = "bold", size = 8.5, hjust = 0),
      legend.position = "none"
    )
})

# Shared legend from a dummy plot
leg_dummy <- ggplot(top_df, aes(x = NES, y = path_short, fill = database)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = DB_COLORS, name = "Database",
                    guide = guide_legend(title.position = "top",
                                         override.aes = list(size = 5))) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title    = element_text(face = "bold", size = 10),
        legend.text     = element_text(face = "bold", size = 9))
shared_legend <- get_legend(leg_dummy)

# Pad to 3 × 4 grid (12 slots for 11 clusters)
n_empty <- 12 - length(plot_list)
for (i in seq_len(n_empty)) plot_list[[length(plot_list) + 1]] <- plot_spacer()

grid <- wrap_plots(plot_list, nrow = 3, ncol = 4) +
  plot_annotation(
    title    = "GSEA — Top 5 Up/Down per Cluster | All Databases | Deduplicated",
    subtitle = "HALLMARK (blue) · C6 (lt.blue) · GO:BP (green) · Reactome (orange) · C7-IMMUNESIGDB (purple)  |  semantic dedup · top 5 by padj",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9)
    )
  )

fig_final <- plot_grid(
  grid,
  shared_legend,
  ncol   = 1,
  rel_heights = c(20, 1)
)

out_pdf <- file.path(OUT_DIR, "05e_gsea_alldb_stacking.pdf")
out_png <- file.path(OUT_DIR, "05e_gsea_alldb_stacking.png")
ggsave(out_pdf, fig_final, width = 20, height = 15)
ggsave(out_png, fig_final, width = 20, height = 15, dpi = 150)
cat("\nSaved:", out_pdf, "\n")
cat("Saved:", out_png, "\n")
cat("\n===== gsea_alldb_stacking.R complete =====\n")
