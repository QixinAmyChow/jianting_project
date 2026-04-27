#!/usr/bin/env Rscript
# Version A — FindSubCluster on existing SNN graph
# Uses the authors' pre-computed graph; no recomputation of PCA/UMAP/neighbors.
# Only rebuilds neighbors after subset() since subsetting breaks graph edges.

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

ZIP     <- "bm_analysis_out/raw_geo/integrated_Seurat_objects.zip"
RDS_IN  <- "integrated_Seurat_objects/47.integrated_object_subset_by_major_celltypes/Myeloid.rds"
OUT_DIR <- "bm_analysis_out/figures_R_v1"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load & extract macrophages ─────────────────────────────────────────────
cat("Extracting Myeloid.rds from ZIP ...\n")
tmp_rds <- tempfile(fileext = ".rds")
unzip(ZIP, files = RDS_IN, exdir = tempdir())
rds_path <- file.path(tempdir(), RDS_IN)

cat("Loading Myeloid.rds ...\n")
s <- readRDS(rds_path)
cat("  Full myeloid object:", ncol(s), "cells\n")

# Keep only macrophages per authors' celltype_C annotation
s <- subset(s, celltype_C == "M\u03c6")
cat("  After Mφ subset:", ncol(s), "cells\n")

DefaultAssay(s) <- if ("SCT" %in% names(s@assays)) "SCT" else "RNA"
cat("  Assay:", DefaultAssay(s), "\n")

# ── 2. Rebuild neighbors on existing PCA (subset breaks the graph) ─────────────
# Reuse dims already stored — just recompute the SNN on the macrophage subset
cat("Rebuilding SNN on macrophage subset (reusing stored PCA) ...\n")
reduction <- if ("pca" %in% names(s@reductions)) "pca" else names(s@reductions)[1]
s <- FindNeighbors(s, reduction = reduction, dims = 1:20,
                   graph.name = c("mac_nn", "mac_snn"))

# ── 3. FindSubCluster on the rebuilt graph ─────────────────────────────────────
cat("Running FindSubCluster ...\n")
# Treat all cells as one cluster to sub-cluster the whole macrophage population
s$dummy_cluster <- "mac"
Idents(s) <- "dummy_cluster"
s <- FindSubCluster(s,
                    cluster        = "mac",
                    graph.name     = "mac_snn",
                    subcluster.name = "mac_subcluster",
                    resolution     = 0.4,
                    algorithm      = 1)
cat("  Sub-clusters:", length(unique(s$mac_subcluster)), "\n")
print(table(s$mac_subcluster))

# ── 4. BC ER+/ER- split ────────────────────────────────────────────────────────
cat("Splitting BC by ESR1 expression ...\n")
esr1 <- FetchData(s, vars = "ESR1")[, 1]
bc_idx <- which(s$cancer == "BC")
bc_esr1 <- esr1[bc_idx]
thresh <- median(bc_esr1[bc_esr1 > 0])
cat("  BC cells:", length(bc_idx), " | ESR1 threshold:", round(thresh, 4), "\n")

s$cancer_type_er <- as.character(s$cancer)
s$cancer_type_er[bc_idx[bc_esr1 >= thresh]] <- "BC_ER+"
s$cancer_type_er[bc_idx[bc_esr1 <  thresh]] <- "BC_ER-"
print(table(s$cancer_type_er))

Idents(s) <- "mac_subcluster"

# ── 5. Figures ─────────────────────────────────────────────────────────────────
CANCER_COLORS <- c(
  "BC_ER+" = "#d62728", "BC_ER-" = "#ff9896",
  "KC"     = "#8c564b", "LC"     = "#17becf",
  "CC"     = "#2ca02c", "BDC"    = "#ff7f0e",
  "EC"     = "#9467bd", "TC"     = "#7f7f7f",
  "PC"     = "#bcbd22", "ctrl"   = "#aec7e8"
)

umap_key <- if (!is.null(s@reductions$umap)) "umap" else names(s@reductions)[grep("umap", names(s@reductions), ignore.case=TRUE)[1]]

## Fig 1: UMAP sub-clusters + cancer type
p1 <- DimPlot(s, reduction = umap_key, group.by = "mac_subcluster",
              label = TRUE, repel = TRUE, pt.size = 0.3) +
      ggtitle("Macrophage sub-clusters (FindSubCluster)") +
      theme(legend.position = "right")

p2 <- DimPlot(s, reduction = umap_key, group.by = "cancer_type_er",
              cols = CANCER_COLORS, pt.size = 0.3) +
      ggtitle("Cancer type (BC split by ER status)") +
      theme(legend.position = "right")

fig1 <- p1 | p2
ggsave(file.path(OUT_DIR, "01_umap_overview.pdf"), fig1, width = 14, height = 6)
ggsave(file.path(OUT_DIR, "01_umap_overview.png"), fig1, width = 14, height = 6, dpi = 150)
cat("  Saved: 01_umap_overview\n")

## Fig 2: Stacked bar — cancer type per cluster
meta <- s@meta.data %>%
  group_by(mac_subcluster, cancer_type_er) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(mac_subcluster) %>%
  mutate(frac = n / sum(n))

col_order <- c("BC_ER+", "BC_ER-",
               sort(setdiff(unique(meta$cancer_type_er), c("BC_ER+","BC_ER-"))))
meta$cancer_type_er <- factor(meta$cancer_type_er, levels = col_order)

fig2 <- ggplot(meta, aes(x = mac_subcluster, y = frac, fill = cancer_type_er)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = CANCER_COLORS, name = "Cancer type") +
  labs(x = "Macrophage sub-cluster", y = "Fraction",
       title = "Cancer type composition per macrophage cluster",
       subtitle = "BC split by ESR1 expression (threshold = median of positive expressors)") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0))

ggsave(file.path(OUT_DIR, "02_cancer_type_bar_bc_er.pdf"), fig2, width = 10, height = 5)
ggsave(file.path(OUT_DIR, "02_cancer_type_bar_bc_er.png"), fig2, width = 10, height = 5, dpi = 150)
cat("  Saved: 02_cancer_type_bar_bc_er\n")

## Fig 3: ESR1 + SP1 feature plots + dot plot
genes <- intersect(c("ESR1","SP1"), rownames(s))

fp <- FeaturePlot(s, features = genes, reduction = umap_key,
                  cols = c("lightgrey","#d62728"), ncol = 2, pt.size = 0.2) &
      theme(legend.position = "right")
ggsave(file.path(OUT_DIR, "03a_feature_plots_ESR1_SP1.pdf"), fp, width = 12, height = 5)
ggsave(file.path(OUT_DIR, "03a_feature_plots_ESR1_SP1.png"), fp, width = 12, height = 5, dpi = 150)

dp <- DotPlot(s, features = genes, group.by = "mac_subcluster") +
      coord_flip() +
      labs(title = "ESR1 & SP1 across macrophage sub-clusters") +
      theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(angle = 0))
ggsave(file.path(OUT_DIR, "03b_dotplot_ESR1_SP1.pdf"), dp, width = 10, height = 4)
ggsave(file.path(OUT_DIR, "03b_dotplot_ESR1_SP1.png"), dp, width = 10, height = 4, dpi = 150)
cat("  Saved: 03a/03b ESR1_SP1\n")

# Report ER+/SP1-high cluster
avg <- AverageExpression(s, features = genes, group.by = "mac_subcluster")[[DefaultAssay(s)]]
scores <- colMeans(avg)
cat("  ER+/SP1-high cluster:", names(which.max(scores)), "\n")

cat("\nAll done. Figures in:", OUT_DIR, "\n")
cat("Method: FindSubCluster (existing PCA reused, SNN rebuilt post-subset)\n")
