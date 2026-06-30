#!/usr/bin/env Rscript
# Update 2 — based on R_v2 full recomputation
# Changes from R_v2:
#   1. ER+/ER- split uses patient metadata (Table S1 Excel), not ESR1 expression
#   2. Highlight genes from highlight_genes.txt instead of ESR1/SP1
# Output: bm_analysis_out/figures_update_2

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(readxl)

ZIP         <- "bm_analysis_out/raw_geo/integrated_Seurat_objects.zip"
RDS_IN      <- "integrated_Seurat_objects/47.integrated_object_subset_by_major_celltypes/Myeloid.rds"
PATIENT_XLS <- "bm_analysis_out/raw_geo/Table S1 patient infor_corrected_2025-08.xlsx"
HIGHLIGHT_F <- "bm_analysis_out/figures/highlight_genes.txt"
OUT_DIR     <- "bm_analysis_out/figures_update_2"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load & extract macrophages ─────────────────────────────────────────────
cat("Extracting Myeloid.rds from ZIP ...\n")
unzip(ZIP, files = RDS_IN, exdir = tempdir())
rds_path <- file.path(tempdir(), RDS_IN)

cat("Loading Myeloid.rds ...\n")
s <- readRDS(rds_path)
cat("  Full myeloid object:", ncol(s), "cells\n")

s <- subset(s, celltype_C == "M\u03c6")
cat("  After Mφ subset:", ncol(s), "cells\n")

# ── 2. Re-normalise + recompute from scratch ───────────────────────────────────
cat("SCTransform (vars.to.regress = percent.mt) ...\n")
s <- SCTransform(s, vars.to.regress = "percent.mt",
                 variable.features.n = 2000,
                 verbose = FALSE)

cat("PCA ...\n")
s <- RunPCA(s, npcs = 50, verbose = FALSE)

cat("UMAP ...\n")
s <- RunUMAP(s, dims = 1:20, verbose = FALSE)

cat("FindNeighbors + FindClusters ...\n")
s <- FindNeighbors(s, dims = 1:20, verbose = FALSE)
s <- FindClusters(s, resolution = 0.4, verbose = FALSE,
                  cluster.name = "mac_subcluster")
cat("  Sub-clusters:", length(unique(s$mac_subcluster)), "\n")
print(table(s$mac_subcluster))

# ── 3. BC ER+/ER- split from patient metadata ─────────────────────────────────
cat("Reading patient metadata for ER status ...\n")
pt <- read_excel(PATIENT_XLS, skip = 1)  # row 1 is a title line
# Column names: MDACC patient ID, cancer.id, cancer, cancer subtype, ...
# The actual header row is row 2 in the file but skip=1 reads row 2 as header
colnames(pt) <- make.names(colnames(pt))

# Identify ER status per cancer.id
# cancer.subtype values: "ER+,PR+,Her2+", "ER-, PR-, HER2+", "na", etc.
er_map <- pt %>%
  filter(cancer == "BC") %>%
  select(cancer.id, cancer.subtype) %>%
  distinct() %>%
  mutate(
    er_status = case_when(
      grepl("^ER\\+",  cancer.subtype, ignore.case = FALSE) ~ "BC_ER+",
      grepl("^ER-|^ER\\s*-", cancer.subtype, ignore.case = FALSE) ~ "BC_ER-",
      TRUE ~ "BC_ER?"
    )
  )

cat("  ER status per cancer.id:\n")
print(er_map[, c("cancer.id", "cancer.subtype", "er_status")])

# Map to cells via cancer.id in Seurat metadata
er_lookup <- setNames(er_map$er_status, er_map$cancer.id)

s$cancer_type_er <- as.character(s$cancer)
bc_cells <- s$cancer == "BC"
s$cancer_type_er[bc_cells] <- er_lookup[s$cancer.id[bc_cells]]
# Cells with unmatched cancer.id (NA) stay as "BC"
s$cancer_type_er[is.na(s$cancer_type_er)] <- "BC"

cat("\n  cancer_type_er distribution:\n")
print(table(s$cancer_type_er))

Idents(s) <- "mac_subcluster"

# ── 4. Read highlight genes ────────────────────────────────────────────────────
cat("\nReading highlight genes from:", HIGHLIGHT_F, "\n")
raw_genes <- readLines(HIGHLIGHT_F)
raw_genes <- unlist(strsplit(paste(raw_genes, collapse = ","), ","))
raw_genes <- trimws(raw_genes)
raw_genes <- raw_genes[nzchar(raw_genes)]
cat("  Requested:", paste(raw_genes, collapse = ", "), "\n")

# Try exact match first, then TOUPPER (human convention)
all_genes  <- rownames(s)
found <- intersect(raw_genes, all_genes)
upper_genes <- toupper(raw_genes)
found_upper <- intersect(upper_genes[!upper_genes %in% found], all_genes)
highlight_genes <- c(found, found_upper)

cat("  Found in object:", paste(highlight_genes, collapse = ", "), "\n")
if (length(highlight_genes) == 0) stop("No highlight genes found in Seurat object.")

gene_label <- paste(highlight_genes, collapse = "/")

# ── 5. Figures ─────────────────────────────────────────────────────────────────
CANCER_COLORS <- c(
  "BC_ER+"  = "#d62728", "BC_ER-"  = "#ff9896", "BC_ER?"  = "#fa9fb5",
  "BC"      = "#fb6a4a",
  "KC"      = "#8c564b", "LC"      = "#17becf",
  "CC"      = "#2ca02c", "BDC"     = "#ff7f0e",
  "EC"      = "#9467bd", "TC"      = "#7f7f7f",
  "PC"      = "#bcbd22", "ctrl"    = "#aec7e8"
)

## Fig 1: UMAP
p1 <- DimPlot(s, reduction = "umap", group.by = "mac_subcluster",
              label = TRUE, repel = TRUE, pt.size = 0.3) +
      ggtitle("Macrophage sub-clusters (recomputed)") +
      theme(legend.position = "right")

p2 <- DimPlot(s, reduction = "umap", group.by = "cancer_type_er",
              cols = CANCER_COLORS, pt.size = 0.3) +
      ggtitle("Cancer type (BC split by ER status from metadata)") +
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

col_order <- c("BC_ER+", "BC_ER-", "BC_ER?", "BC",
               sort(setdiff(unique(meta$cancer_type_er),
                            c("BC_ER+","BC_ER-","BC_ER?","BC"))))
meta$cancer_type_er <- factor(meta$cancer_type_er, levels = col_order)

fig2 <- ggplot(meta, aes(x = mac_subcluster, y = frac, fill = cancer_type_er)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = CANCER_COLORS, name = "Cancer type") +
  labs(x = "Macrophage sub-cluster", y = "Fraction",
       title = "Cancer type composition per macrophage cluster",
       subtitle = "BC split by ER status from patient metadata (Table S1)") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0))

ggsave(file.path(OUT_DIR, "02_cancer_type_bar_bc_er.pdf"), fig2, width = 10, height = 5)
ggsave(file.path(OUT_DIR, "02_cancer_type_bar_bc_er.png"), fig2, width = 10, height = 5, dpi = 150)
cat("  Saved: 02_cancer_type_bar_bc_er\n")

## Fig 3a: Feature plots — highlight genes
fp <- FeaturePlot(s, features = highlight_genes, reduction = "umap",
                  cols = c("lightgrey","#d62728"),
                  ncol = min(length(highlight_genes), 3),
                  pt.size = 0.2) &
      theme(legend.position = "right")

fig3a_w <- 6 * min(length(highlight_genes), 3)
fig3a_h <- 5 * ceiling(length(highlight_genes) / 3)
ggsave(file.path(OUT_DIR, "03a_feature_plots_highlight.pdf"), fp,
       width = fig3a_w, height = fig3a_h)
ggsave(file.path(OUT_DIR, "03a_feature_plots_highlight.png"), fp,
       width = fig3a_w, height = fig3a_h, dpi = 150)

## Fig 3b: Dot plot — highlight genes
dp <- DotPlot(s, features = highlight_genes, group.by = "mac_subcluster") +
      coord_flip() +
      labs(title = paste(gene_label, "across macrophage sub-clusters")) +
      theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(angle = 0))

dp_h <- max(4, 0.35 * length(highlight_genes) + 2)
ggsave(file.path(OUT_DIR, "03b_dotplot_highlight.pdf"), dp, width = 10, height = dp_h)
ggsave(file.path(OUT_DIR, "03b_dotplot_highlight.png"), dp, width = 10, height = dp_h, dpi = 150)
cat("  Saved: 03a/03b highlight genes\n")

# Report which cluster has highest average expression of highlight genes
avg <- AverageExpression(s, features = highlight_genes, group.by = "mac_subcluster")[["SCT"]]
scores <- colMeans(avg)
cat("  Highlight-gene-high cluster:", names(which.max(scores)), "\n")
cat("  Scores per cluster:\n")
print(round(sort(scores, decreasing = TRUE), 4))

# ── Save Seurat object for downstream DEG / GSEA ──────────────────────────────
rds_out <- file.path(OUT_DIR, "mac_update2.rds")
cat("Saving Seurat object to:", rds_out, "...\n")
saveRDS(s, rds_out)

cat("\nAll done. Figures in:", OUT_DIR, "\n")
cat("Method: Full recompute — SCTransform + PCA + UMAP + FindNeighbors + FindClusters\n")
cat("ER split: metadata-based (Table S1), not expression-based\n")
cat("Marker genes:", gene_label, "\n")
