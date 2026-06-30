#!/usr/bin/env Rscript
# DEG + GSEA + metadata breakdown for macrophage sub-clusters (update 2)
# Loads cached Seurat object from mac_subcluster_update2.R

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(pheatmap)
library(fgsea)
library(msigdbr)
library(RColorBrewer)
library(readxl)

options(future.globals.maxSize = 2 * 1024^3)

RDS_IN      <- "bm_analysis_out/figures_update_2/mac_update2.rds"
PATIENT_XLS <- "bm_analysis_out/raw_geo/Table S1 patient infor_corrected_2025-08.xlsx"
OUT_DIR     <- "bm_analysis_out/figures_update_2"

if (!file.exists(RDS_IN)) stop("Run mac_subcluster_update2.R first.")

cat("Loading cached Seurat object ...\n")
s <- readRDS(RDS_IN)
Idents(s) <- "mac_subcluster"
cat("  Clusters:", length(levels(s)), "| Cells:", ncol(s), "\n")

# ── Join treatment response from patient table ─────────────────────────────────
cat("Joining treatment response metadata ...\n")
pt <- read_excel(PATIENT_XLS, skip = 1)
colnames(pt) <- make.names(colnames(pt))

tx_map <- pt %>%
  select(cancer.id,
         neo_resp = Neoadjuvant.treatment.response,
         met_resp = Metastatic...treatment.response) %>%
  distinct() %>%
  mutate(
    tx_response = case_when(
      !tolower(neo_resp) %in% c("na","") ~ tolower(neo_resp),
      !tolower(met_resp) %in% c("na","") ~ tolower(met_resp),
      TRUE ~ "na"
    )
  )

tx_lookup <- setNames(tx_map$tx_response, tx_map$cancer.id)
s$tx_response <- unname(tx_lookup[s$cancer.id])
s$tx_response[is.na(s$tx_response)] <- "na"
cat("  tx_response distribution:\n")
print(table(s$tx_response))

# Standardise gender
s$gender_clean <- recode(s$gender, "f" = "Female", "m" = "Male", "na" = "Unknown")

# ── 1. DEG — FindAllMarkers (cached) ──────────────────────────────────────────
deg_csv <- file.path(OUT_DIR, "deg_all_markers.csv")
if (file.exists(deg_csv)) {
  cat("\nLoading cached DEG results ...\n")
  markers <- read.csv(deg_csv)
  markers$cluster <- as.character(markers$cluster)
} else {
  cat("\nRunning FindAllMarkers ...\n")
  markers <- FindAllMarkers(s, only.pos = TRUE, min.pct = 0.25,
                            logfc.threshold = 0.25, verbose = FALSE)
  markers <- markers %>% filter(p_val_adj < 0.05) %>%
    arrange(cluster, desc(avg_log2FC))
  write.csv(markers, deg_csv, row.names = FALSE)
}
cat("  Total markers:", nrow(markers), "\n")

## Fig 4a: DEG heatmap — top 5 per cluster
if (!file.exists(file.path(OUT_DIR, "04a_deg_heatmap.png"))) {
  top5 <- markers %>% group_by(cluster) %>%
    slice_max(avg_log2FC, n = 5) %>% ungroup()
  hm <- DoHeatmap(s, features = unique(top5$gene), group.by = "mac_subcluster",
                  angle = 0, size = 3, draw.lines = TRUE) +
    scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu")), na.value = "white") +
    theme(axis.text.y = element_text(size = 7))
  ggsave(file.path(OUT_DIR, "04a_deg_heatmap.pdf"), hm, width = 14, height = 8)
  ggsave(file.path(OUT_DIR, "04a_deg_heatmap.png"), hm, width = 14, height = 8, dpi = 150)
  cat("  Saved: 04a_deg_heatmap\n")
} else { cat("  04a_deg_heatmap exists, skipping.\n") }

## Fig 4b: DEG dot plot — top 3 per cluster
if (!file.exists(file.path(OUT_DIR, "04b_deg_dotplot.png"))) {
  top3 <- markers %>% group_by(cluster) %>%
    slice_max(avg_log2FC, n = 3) %>% ungroup()
  dp <- DotPlot(s, features = unique(top3$gene), group.by = "mac_subcluster") +
    coord_flip() +
    labs(title = "Top 3 DEGs per macrophage sub-cluster") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 0), axis.text.y = element_text(size = 8))
  ggsave(file.path(OUT_DIR, "04b_deg_dotplot.pdf"), dp, width = 12, height = 8)
  ggsave(file.path(OUT_DIR, "04b_deg_dotplot.png"), dp, width = 12, height = 8, dpi = 150)
  cat("  Saved: 04b_deg_dotplot\n")
} else { cat("  04b_deg_dotplot exists, skipping.\n") }

# ── 2. GSEA — cancer-focused gene sets (Hallmark H + Oncogenic C6) ────────────
cat("\nLoading cancer-focused MSigDB gene sets (Hallmark + Oncogenic C6) ...\n")
msig_h  <- msigdbr(species = "Homo sapiens", collection = "H")
msig_c6 <- msigdbr(species = "Homo sapiens", collection = "C6")
msig    <- bind_rows(msig_h, msig_c6)
pathways <- split(msig$gene_symbol, msig$gs_name)
cat("  Total gene sets:", length(pathways), "\n")

cat("Running fgsea per cluster (HVGs only) ...\n")
clusters <- levels(Idents(s))
gsea_results <- list()

for (cl in clusters) {
  cat("  Cluster", cl, "...\n")
  cl_markers <- FindMarkers(s, ident.1 = cl,
                            features        = VariableFeatures(s),
                            only.pos        = FALSE,
                            min.pct         = 0.1,
                            logfc.threshold = 0,
                            verbose         = FALSE)
  ranks <- setNames(cl_markers$avg_log2FC, rownames(cl_markers))
  ranks <- sort(ranks, decreasing = TRUE)

  res <- fgsea(pathways = pathways, stats = ranks,
               minSize = 15, maxSize = 500, nPermSimple = 1000, nproc = 1)
  res$cluster <- cl
  gsea_results[[cl]] <- res
}

gsea_all <- bind_rows(gsea_results)
write.csv(gsea_all %>% select(-leadingEdge),
          file.path(OUT_DIR, "gsea_cancer_results.csv"), row.names = FALSE)

## Fig 5a: GSEA NES heatmap
cat("  Drawing GSEA NES heatmap ...\n")
sig_paths <- gsea_all %>%
  filter(padj < 0.1) %>%
  group_by(pathway) %>% filter(n() >= 2) %>%
  pull(pathway) %>% unique()

if (length(sig_paths) < 5) {
  sig_paths <- gsea_all %>%
    group_by(pathway) %>%
    summarise(mean_abs_nes = mean(abs(NES), na.rm = TRUE)) %>%
    slice_max(mean_abs_nes, n = 25) %>% pull(pathway)
}

nes_mat <- gsea_all %>%
  filter(pathway %in% sig_paths) %>%
  select(pathway, cluster, NES) %>%
  pivot_wider(names_from = cluster, values_from = NES, values_fill = 0) %>%
  as.data.frame() %>% { rownames(.) <- .$pathway; .[ , -1] }

rownames(nes_mat) <- gsub("^HALLMARK_", "H: ", rownames(nes_mat))
rownames(nes_mat) <- gsub("_", " ", rownames(nes_mat))

pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
hm_h <- max(6, 0.28 * nrow(nes_mat) + 2)

pdf(file.path(OUT_DIR, "05a_gsea_nes_heatmap.pdf"), width = 14, height = hm_h)
pheatmap(nes_mat, color = pal, breaks = seq(-3, 3, length.out = 101),
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize_row = 8, fontsize_col = 10, border_color = NA,
         main = "GSEA NES - Hallmark + Oncogenic (C6) per macrophage cluster",
         na_col = "grey90")
dev.off()
png(file.path(OUT_DIR, "05a_gsea_nes_heatmap.png"), width = 14, height = hm_h,
    units = "in", res = 150)
pheatmap(nes_mat, color = pal, breaks = seq(-3, 3, length.out = 101),
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize_row = 8, fontsize_col = 10, border_color = NA,
         main = "GSEA NES - Hallmark + Oncogenic (C6) per macrophage cluster",
         na_col = "grey90")
dev.off()
cat("  Saved: 05a_gsea_nes_heatmap\n")

## Fig 5b: GSEA bubble plot — top pathways per cluster
top_gsea <- gsea_all %>%
  filter(padj < 0.25) %>%
  group_by(cluster) %>%
  slice_max(abs(NES), n = 3) %>%
  ungroup() %>%
  mutate(
    pathway_short = gsub("^HALLMARK_", "H: ", pathway),
    pathway_short = gsub("_", " ", pathway_short)
  )

if (nrow(top_gsea) > 0) {
  bp <- ggplot(top_gsea, aes(x = cluster, y = pathway_short,
                              size = -log10(pmax(padj, 1e-4)), color = NES)) +
    geom_point(alpha = 0.85) +
    scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu")),
                          limits = c(-3, 3), oob = scales::squish, name = "NES") +
    scale_size_continuous(name = "-log10(padj)", range = c(2, 8)) +
    labs(x = "Macrophage sub-cluster", y = NULL,
         title = "Top GSEA pathways per cluster (Hallmark + Oncogenic C6, padj<0.25)") +
    theme_classic(base_size = 12) +
    theme(axis.text.y = element_text(size = 9))
  bub_h <- max(5, 0.3 * length(unique(top_gsea$pathway_short)) + 2)
  ggsave(file.path(OUT_DIR, "05b_gsea_bubble.pdf"), bp, width = 12, height = bub_h)
  ggsave(file.path(OUT_DIR, "05b_gsea_bubble.png"), bp, width = 12, height = bub_h, dpi = 150)
  cat("  Saved: 05b_gsea_bubble\n")
}

# ── 3. Metadata breakdown per cluster ─────────────────────────────────────────
cat("\nGenerating metadata breakdown figures ...\n")

meta_bar <- function(df, fill_var, title, subtitle, fill_colors = NULL, filename) {
  dat <- df %>%
    group_by(mac_subcluster, .data[[fill_var]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(mac_subcluster) %>%
    mutate(frac = n / sum(n))

  p <- ggplot(dat, aes(x = mac_subcluster, y = frac, fill = .data[[fill_var]])) +
    geom_bar(stat = "identity", width = 0.8) +
    labs(x = "Macrophage sub-cluster", y = "Fraction",
         title = title, subtitle = subtitle, fill = fill_var) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 0))

  if (!is.null(fill_colors))
    p <- p + scale_fill_manual(values = fill_colors)
  else
    p <- p + scale_fill_brewer(palette = "Set2")

  ggsave(file.path(OUT_DIR, paste0(filename, ".pdf")), p, width = 10, height = 5)
  ggsave(file.path(OUT_DIR, paste0(filename, ".png")), p, width = 10, height = 5, dpi = 150)
  cat("  Saved:", filename, "\n")
}

## Fig 6: Gender per cluster
meta_bar(s@meta.data, "gender_clean",
         title    = "Gender distribution per macrophage sub-cluster",
         subtitle = "Female / Male / Unknown",
         fill_colors = c("Female" = "#e377c2", "Male" = "#1f77b4", "Unknown" = "#aec7e8"),
         filename = "06_gender_per_cluster")

## Fig 7: Tissue origin per cluster
meta_bar(s@meta.data, "tissue.origin",
         title    = "Tissue / bone site per macrophage sub-cluster",
         subtitle = "Site of bone metastasis resection",
         filename = "07_tissue_origin_per_cluster")

## Fig 8: Treatment response per cluster
meta_bar(s@meta.data, "tx_response",
         title    = "Treatment response per macrophage sub-cluster",
         subtitle = "Neoadjuvant or metastatic treatment response (Table S1)",
         filename = "08_treatment_per_cluster")

cat("\nAll done.\n")
cat("Figures in:", OUT_DIR, "\n")
