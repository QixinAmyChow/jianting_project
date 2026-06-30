# =============================================================================
# Reyfman 2019 Annotation Pipeline
# Dataset : Reyfman et al. 2019, AJRCCM — GSE122960
#           17 samples: 8 Donor (Control) + 4 IPF + 1 HP + 2 SSc-ILD +
#                       1 Myositis-ILD + 1 Cryobiopsy
#
# Approach: Reyfman QC/normalization params + Harmony batch correction (v5-compatible)
#   1. Per-sample QC with exact paper thresholds
#   2. Normalize + HVG per sample
#   3. Merge all samples → JoinLayers → PCA
#   4. Harmony integration by sample_id
#   5. Cluster (res=0.5), tSNE + UMAP on harmony embedding
#   6. FindAllMarkers → auto-annotation by canonical markers
#   7. Subset to Donor + IPF only → save
#
# Cache   : ipf/Reyfman2019/cache_reyfman_*.rds
# Figures : ipf/figures_Reyfman2019/
# =============================================================================

library(Seurat)
library(harmony)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)

# ─────────────────────────────────────────────────────────────────────────────
# PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────

DATA_DIR      <- "ipf/Reyfman2019"
OUT_DIR       <- "ipf/figures_Reyfman2019"
CACHE_ALL     <- file.path(DATA_DIR, "cache_reyfman_all_annotated.rds")
CACHE_IPF     <- file.path(DATA_DIR, "cache_reyfman_ipf_annotated.rds")
CACHE_PRE_HARMONY <- file.path(DATA_DIR, "cache_reyfman_pre_harmony.rds")
MARKERS_F     <- file.path(DATA_DIR, "reyfman_cluster_markers.csv")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Per-sample QC thresholds (exact from Reyfman R code)
SAMPLE_META <- data.frame(
  file = c(
    "GSM3489182_Donor_01_filtered_gene_bc_matrices_h5.h5",
    "GSM3489185_Donor_02_filtered_gene_bc_matrices_h5.h5",
    "GSM3489187_Donor_03_filtered_gene_bc_matrices_h5.h5",
    "GSM3489189_Donor_04_filtered_gene_bc_matrices_h5.h5",
    "GSM3489191_Donor_05_filtered_gene_bc_matrices_h5.h5",
    "GSM3489193_Donor_06_filtered_gene_bc_matrices_h5.h5",
    "GSM3489195_Donor_07_filtered_gene_bc_matrices_h5.h5",
    "GSM3489197_Donor_08_filtered_gene_bc_matrices_h5.h5",
    "GSM3489183_IPF_01_filtered_gene_bc_matrices_h5.h5",
    "GSM3489184_IPF_02_filtered_gene_bc_matrices_h5.h5",
    "GSM3489188_IPF_03_filtered_gene_bc_matrices_h5.h5",
    "GSM3489190_IPF_04_filtered_gene_bc_matrices_h5.h5",
    "GSM3489192_HP_01_filtered_gene_bc_matrices_h5.h5",
    "GSM3489194_SSc-ILD_01_filtered_gene_bc_matrices_h5.h5",
    "GSM3489196_Myositis-ILD_01_filtered_gene_bc_matrices_h5.h5",
    "GSM3489198_SSc-ILD_02_filtered_gene_bc_matrices_h5.h5",
    "GSM3489186_Cryobiopsy_01_filtered_gene_bc_matrices_h5.h5"
  ),
  sample_id = c(
    "donor01","donor02","donor03","donor04","donor05","donor06","donor07","donor08",
    "ipf01","ipf02","ipf03","ipf04",
    "hp01","ssc_ild01","myositis_ild01","ssc_ild02","cryo"
  ),
  condition = c(
    rep("Donor", 8), rep("IPF", 4),
    "HP","SSc-ILD","Myositis-ILD","SSc-ILD","Cryo"
  ),
  group = c(rep("donor", 8), rep("fibrosis", 9)),
  nGene_max = c(5000,6000,5000,5000,6000,6000,5000,6000,
                6000,6000,5000,2500,
                6000,6000,4000,6000,5000),
  mito_max  = c(0.20,0.15,0.10,0.15,0.10,0.10,0.10,0.20,
                0.20,0.20,0.20,0.20,
                0.15,0.15,0.30,0.30,0.20),
  stringsAsFactors = FALSE
)

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1-3: Load, QC, normalize, merge, PCA
# ─────────────────────────────────────────────────────────────────────────────

if (file.exists(CACHE_PRE_HARMONY)) {
  message("Loading cached pre-harmony object...")
  hl <- readRDS(CACHE_PRE_HARMONY)
} else {
  message("=== Loading and processing samples ===")
  obj_list <- list()

  for (i in seq_len(nrow(SAMPLE_META))) {
    sid   <- SAMPLE_META$sample_id[i]
    fpath <- file.path(DATA_DIR, SAMPLE_META$file[i])
    message(sprintf("[%d/%d] %s", i, nrow(SAMPLE_META), sid))

    counts <- Read10X_h5(fpath)
    obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200,
                              project = sid)
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj <- subset(obj,
      subset = nFeature_RNA >= 200 &
               nFeature_RNA <= SAMPLE_META$nGene_max[i] &
               percent.mt  <= SAMPLE_META$mito_max[i] * 100)

    obj$sample_id <- sid
    obj$condition <- SAMPLE_META$condition[i]
    obj$group     <- SAMPLE_META$group[i]

    obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                         scale.factor = 10000, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = 2000, verbose = FALSE)
    obj_list[[sid]] <- obj
    message(sprintf("  -> %d cells", ncol(obj)))
  }

  message(sprintf("Total cells: %d", sum(sapply(obj_list, ncol))))

  # Merge all samples
  message("=== Merging all samples ===")
  hl <- merge(obj_list[[1]],
              y = obj_list[-1],
              add.cell.ids = SAMPLE_META$sample_id)
  hl <- JoinLayers(hl)

  # HVG on full merged object
  hl <- FindVariableFeatures(hl, nfeatures = 2000, verbose = FALSE)
  hl <- ScaleData(hl, verbose = FALSE)
  hl <- RunPCA(hl, npcs = 40, verbose = FALSE)

  saveRDS(hl, CACHE_PRE_HARMONY)
  message(sprintf("Pre-harmony object cached: %s", CACHE_PRE_HARMONY))
  rm(obj_list); gc()
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: Harmony integration by sample
# ─────────────────────────────────────────────────────────────────────────────

message("=== Harmony integration ===")
hl <- RunHarmony(hl, group.by.vars = "sample_id", reduction.use = "pca",
                 dims.use = 1:40, verbose = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: Cluster + UMAP + tSNE on harmony embedding
# ─────────────────────────────────────────────────────────────────────────────

message("=== Clustering ===")
hl <- FindNeighbors(hl, reduction = "harmony", dims = 1:31, verbose = FALSE)
hl <- FindClusters(hl, resolution = 0.5, verbose = FALSE)
hl <- RunUMAP(hl, reduction = "harmony", dims = 1:31)
hl <- RunTSNE(hl, reduction = "harmony", dims = 1:31)
message(sprintf("Clusters: %d", length(unique(Idents(hl)))))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: FindAllMarkers + annotate
# ─────────────────────────────────────────────────────────────────────────────

message("=== Finding cluster markers ===")
markers <- FindAllMarkers(hl, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
write.csv(markers, MARKERS_F, row.names = FALSE)
message(sprintf("Markers saved: %s", MARKERS_F))

CANONICAL <- list(
  "AT2 Cells"           = c("SFTPC", "SFTPB", "SFTPD"),
  "AT1 Cells"           = c("AGER", "RTKN2", "CLDN18"),
  "Club Cells"          = c("SCGB3A2", "SCGB1A1"),
  "Ciliated Cells"      = c("TPPP3", "FOXJ1", "CAPS"),
  "Basal Cells"         = c("KRT5", "TP63", "KRT17"),
  "Macrophages"         = c("CD68", "MRC1", "FABP4", "MARCO"),
  "Monocytes"           = c("FCN1", "CD14", "S100A8"),
  "Dendritic Cells"     = c("CLEC10A", "HLA-DQA1", "CD1C"),
  "T/NKT Cells"         = c("CD3D", "CD3E", "TRAC"),
  "B Cells"             = c("MS4A1", "CD79A"),
  "Plasma Cells"        = c("IGHG4", "MZB1", "JCHAIN"),
  "Mast Cells"          = c("TPSB2", "TPSAB1", "CPA3"),
  "Endothelial Cells"   = c("VWF", "PECAM1", "CDH5"),
  "Fibroblasts"         = c("DCN", "COL1A1", "LUM", "PDGFRA"),
  "Smooth Muscle Cells" = c("ACTA2", "MYH11", "CNN1")
)

cluster_ids <- levels(Idents(hl))
score_mat <- sapply(CANONICAL, function(mgs) {
  sapply(cluster_ids, function(cl) {
    sub_m <- markers[markers$cluster == cl & markers$gene %in% mgs, ]
    if (nrow(sub_m) == 0) return(0)
    mean(sub_m$avg_log2FC)
  })
})

auto_labels <- apply(score_mat, 1, function(x) {
  if (max(x, na.rm = TRUE) == 0) return("Unassigned")
  names(which.max(x))
})

message("Auto-annotation:")
for (cl in cluster_ids) message(sprintf("  Cluster %s -> %s", cl, auto_labels[cl]))

new_idents <- auto_labels[as.character(Idents(hl))]
names(new_idents) <- colnames(hl)
hl <- AddMetaData(hl, metadata = new_idents, col.name = "cell_type")
Idents(hl) <- "cell_type"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 7: Save all-sample object + plots
# ─────────────────────────────────────────────────────────────────────────────

saveRDS(hl, CACHE_ALL)
message(sprintf("All-sample object saved: %s", CACHE_ALL))

pdf(file.path(OUT_DIR, "all_umap_condition.pdf"), width = 8, height = 6)
print(DimPlot(hl, reduction = "umap", group.by = "condition", pt.size = 0.3) +
      ggtitle("All samples — condition"))
dev.off()

pdf(file.path(OUT_DIR, "all_umap_celltype.pdf"), width = 10, height = 7)
print(DimPlot(hl, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.3) +
      ggtitle("All samples — cell type"))
dev.off()

# ─────────────────────────────────────────────────────────────────────────────
# STEP 8: Subset to Donor + IPF only
# ─────────────────────────────────────────────────────────────────────────────

message("=== Subsetting to Donor + IPF ===")
hl_ipf <- subset(hl, subset = condition %in% c("Donor", "IPF"))
hl_ipf$condition <- factor(hl_ipf$condition, levels = c("Donor", "IPF"))
hl_ipf <- RunUMAP(hl_ipf, reduction = "harmony", dims = 1:31)
hl_ipf <- RunTSNE(hl_ipf, reduction = "harmony", dims = 1:31)

message(sprintf("IPF+Donor subset: %d cells", ncol(hl_ipf)))
print(table(hl_ipf$condition))
print(table(hl_ipf$cell_type))

saveRDS(hl_ipf, CACHE_IPF)
message(sprintf("IPF+Donor object saved: %s", CACHE_IPF))

pdf(file.path(OUT_DIR, "ipf_umap_condition.pdf"), width = 7, height = 6)
print(DimPlot(hl_ipf, reduction = "umap", group.by = "condition",
              cols = c("Donor" = "#2166AC", "IPF" = "#D6604D"), pt.size = 0.4) +
      ggtitle("Donor vs IPF"))
dev.off()

pdf(file.path(OUT_DIR, "ipf_umap_celltype.pdf"), width = 10, height = 7)
print(DimPlot(hl_ipf, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.4) +
      ggtitle("Cell types — Donor + IPF"))
dev.off()

pdf(file.path(OUT_DIR, "ipf_fibro_markers.pdf"), width = 10, height = 5)
print(FeaturePlot(hl_ipf, features = c("DCN", "COL1A1", "ACTA2", "SFTPC"),
                  reduction = "umap", ncol = 4, pt.size = 0.3))
dev.off()

message("=== Done ===")
message(sprintf("Outputs in: %s", OUT_DIR))
