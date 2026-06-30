#!/usr/bin/env Rscript
# Sanity check — macrophage cell counts & % per cancer type and per sample
# Loads mac_update3_best.rds (mac subset only)
# Output: bm_analysis_out/figures_update_3/00_sanity_*

library(Seurat)
library(ggplot2)
library(dplyr)

RDS_IN  <- "bm_analysis_out/figures_update_3/mac_update3_best.rds"
OUT_DIR <- "bm_analysis_out/figures_update_3"

cat("Loading ...\n")
s <- readRDS(RDS_IN)
meta <- s@meta.data
cat("  Total macrophages:", nrow(meta), "\n")

BOLD_LEG <- theme(
  legend.title = element_text(face = "bold", size = 11),
  legend.text  = element_text(face = "bold", size = 10)
)

CANCER_COLORS <- c(
  "BC_ER+"  = "#d62728", "BC_ER-"  = "#ff9896", "BC_ER?"  = "#fa9fb5",
  "BC"      = "#fb6a4a",
  "KC"      = "#8c564b", "LC"      = "#17becf",
  "CC"      = "#2ca02c", "BDC"     = "#ff7f0e",
  "EC"      = "#9467bd", "TC"      = "#7f7f7f",
  "PC"      = "#bcbd22", "ctrl"    = "#aec7e8"
)

# Use cancer_type_er if present, else cancer
ct_col <- if ("cancer_type_er" %in% colnames(meta)) "cancer_type_er" else "cancer"
meta$ct <- meta[[ct_col]]

# ── Fig A: cells per cancer type ─────────────────────────────────────────────
ct_sum <- meta %>%
  group_by(ct) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pct   = 100 * n / sum(n),
         label = sprintf("%d\n(%.1f%%)", n, pct)) %>%
  arrange(desc(n))

ct_sum$ct <- factor(ct_sum$ct, levels = ct_sum$ct)

figA <- ggplot(ct_sum, aes(x = ct, y = n, fill = ct)) +
  geom_bar(stat = "identity", width = 0.75) +
  geom_text(aes(label = label), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = CANCER_COLORS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Cancer type", y = "Macrophage cell count",
       title = "Macrophage cells per cancer type",
       subtitle = sprintf("Total = %d macrophages across %d cancer types",
                          nrow(meta), n_distinct(meta$ct))) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0)) + BOLD_LEG

ggsave(file.path(OUT_DIR, "00_sanity_cells_per_cancer.pdf"), figA, width = 9, height = 5)
ggsave(file.path(OUT_DIR, "00_sanity_cells_per_cancer.png"), figA, width = 9, height = 5, dpi = 150)
cat("  Saved: 00_sanity_cells_per_cancer\n")

# ── Fig B: cells per sample (cancer.id), colored by cancer type ───────────────
# Order samples: by cancer type, then descending n within type
samp_sum <- meta %>%
  group_by(cancer.id, ct) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(ct, desc(n))

# Fix factor order for x axis
samp_sum$cancer.id <- factor(samp_sum$cancer.id,
                              levels = unique(samp_sum$cancer.id))

n_samples <- n_distinct(meta$cancer.id)
bar_w <- max(8, 0.25 * n_samples + 2)

figB <- ggplot(samp_sum, aes(x = cancer.id, y = n, fill = ct)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = CANCER_COLORS, name = "Cancer type") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Sample (cancer.id)", y = "Macrophage cell count",
       title = "Macrophage cells per sample",
       subtitle = sprintf("%d samples | ordered by cancer type then descending count", n_samples)) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) + BOLD_LEG

ggsave(file.path(OUT_DIR, "00_sanity_cells_per_sample.pdf"), figB, width = bar_w, height = 5)
ggsave(file.path(OUT_DIR, "00_sanity_cells_per_sample.png"), figB, width = bar_w, height = 5, dpi = 150)
cat("  Saved: 00_sanity_cells_per_sample\n")

# ── Fig C: combined panel (A | B) ─────────────────────────────────────────────
library(patchwork)
figC <- figA + figB + plot_layout(widths = c(1, 2))
ggsave(file.path(OUT_DIR, "00_sanity_combined.pdf"), figC, width = bar_w + 4, height = 5)
ggsave(file.path(OUT_DIR, "00_sanity_combined.png"), figC, width = bar_w + 4, height = 5, dpi = 150)
cat("  Saved: 00_sanity_combined\n")

# ── Print summary table ───────────────────────────────────────────────────────
cat("\nMacrophage cells per cancer type:\n")
print(ct_sum %>% select(ct, n, pct) %>% mutate(pct = round(pct, 1)), row.names = FALSE)

cat("\nMacrophage cells per sample:\n")
print(samp_sum %>% arrange(ct, cancer.id), row.names = FALSE)

cat("\nDone.\n")
