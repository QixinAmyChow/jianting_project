# =============================================================================
# Reyfman 2019 Annotation Pipeline — v2
# Architecture matches paper's S1 code:
#   - Per-sample: QC → normalize → HVG (dispersion) → ScaleData(regress nUMI+mito)
#                 → PCA → cluster → FindAllMarkers → annotate
#   - Merge annotated objects → Harmony (v5 replacement for CCA, visualization only)
#   - Subset to Donor + IPF
#
# Key differences from v1:
#   1. ScaleData regresses nCount_RNA + percent.mt (paper: nUMI + percent.mito)
#   2. HVG via mean.var.plot dispersion (paper: ExpMean/LogVMR cutoffs)
#   3. Per-sample PCA dims and cluster resolution exactly from paper's S1 code
#   4. Annotation done per-sample on individual clusterings, not on joint clusters
#   5. COL3A1 added to Fibroblast markers (paper uses this)
# =============================================================================

library(Seurat)
library(harmony)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)

DATA_DIR  <- "ipf/data/Reyfman2019"
OUT_DIR   <- "ipf/figures/Reyfman2019/v2"
CACHE_LIST <- file.path(DATA_DIR, "cache_v2_per_sample_list.rds")
CACHE_ALL  <- file.path(DATA_DIR, "cache_v2_all_annotated.rds")
CACHE_IPF  <- file.path(DATA_DIR, "cache_v2_ipf_annotated.rds")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Per-sample metadata: QC thresholds + exact PCA dims + cluster resolution from S1 code
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
  # PCA dims per sample (exact from paper's S1 code; cryo not in S1, use 15)
  pca_dims  = c(13,16,9,9,27,20,22,12,
                26,26,36,11,
                32,19,17,25,15),
  # Cluster resolution per sample (exact from paper's S1 code)
  cluster_res = c(0.3,0.3,0.3,0.3,0.5,0.5,0.3,0.5,
                  0.5,0.7,0.5,0.5,
                  0.7,0.5,0.7,0.5,0.5),
  stringsAsFactors = FALSE
)

# ─────────────────────────────────────────────────────────────────────────────
# Canonical markers — matches paper's annotation markers
# COL3A1 added for Fibroblasts (paper uses this, v1 was missing it)
# ─────────────────────────────────────────────────────────────────────────────
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
  "Fibroblasts"         = c("DCN", "COL1A1", "COL3A1", "LUM", "PDGFRA"),
  "Smooth Muscle Cells" = c("ACTA2", "MYH11", "CNN1")
)

annotate_sample <- function(obj, sid) {
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                            logfc.threshold = 0.25, verbose = FALSE)
  cluster_ids <- levels(Idents(obj))
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
  message(sprintf("  [%s] cluster → cell type:", sid))
  for (cl in cluster_ids) message(sprintf("    %s → %s", cl, auto_labels[cl]))
  new_ids <- auto_labels[as.character(Idents(obj))]
  names(new_ids) <- colnames(obj)
  obj <- AddMetaData(obj, metadata = new_ids, col.name = "cell_type")
  obj
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1–4: Per-sample pipeline (paper's S1 architecture)
# ─────────────────────────────────────────────────────────────────────────────

if (file.exists(CACHE_LIST)) {
  message("Loading cached per-sample annotated list...")
  obj_list <- readRDS(CACHE_LIST)
} else {
  message("=== Per-sample processing (matching paper S1) ===")
  obj_list <- list()

  for (i in seq_len(nrow(SAMPLE_META))) {
    sid   <- SAMPLE_META$sample_id[i]
    fpath <- file.path(DATA_DIR, SAMPLE_META$file[i])
    dims  <- seq_len(SAMPLE_META$pca_dims[i])
    res   <- SAMPLE_META$cluster_res[i]
    message(sprintf("\n[%d/17] %s (dims=1:%d, res=%.1f)", i, sid, max(dims), res))

    counts <- Read10X_h5(fpath)
    obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200,
                              project = sid)
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj <- subset(obj,
      subset = nFeature_RNA >= 200 &
               nFeature_RNA <= SAMPLE_META$nGene_max[i] &
               percent.mt  <= SAMPLE_META$mito_max[i] * 100)
    message(sprintf("  After QC: %d cells", ncol(obj)))

    obj$sample_id <- sid
    obj$condition <- SAMPLE_META$condition[i]
    obj$group     <- SAMPLE_META$group[i]

    obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                         scale.factor = 10000, verbose = FALSE)

    # HVG: dispersion-based (paper uses ExpMean/LogVMR cutoffs in Seurat v2)
    # v5 equivalent: mean.var.plot with same cutoffs
    obj <- FindVariableFeatures(obj,
                                selection.method  = "mean.var.plot",
                                mean.cutoff       = c(0.0125, 3),
                                dispersion.cutoff = c(0.5, Inf),
                                verbose = FALSE)
    message(sprintf("  HVG: %d", length(VariableFeatures(obj))))

    # Regress out nCount_RNA (nUMI) and percent.mt (percent.mito) — paper does this
    obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA", "percent.mt"),
                     verbose = FALSE)
    obj <- RunPCA(obj, npcs = max(dims) + 5, verbose = FALSE)

    obj <- FindNeighbors(obj, dims = dims, verbose = FALSE)
    obj <- FindClusters(obj, resolution = res, verbose = FALSE)
    message(sprintf("  Clusters: %d", length(unique(Idents(obj)))))

    # Annotate per sample
    obj <- annotate_sample(obj, sid)

    obj_list[[sid]] <- obj
  }

  saveRDS(obj_list, CACHE_LIST)
  message(sprintf("\nPer-sample list cached: %s", CACHE_LIST))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: Print per-sample cell type summary
# ─────────────────────────────────────────────────────────────────────────────

message("\n=== Per-sample cell type counts ===")
for (sid in names(obj_list)) {
  ct <- table(obj_list[[sid]]$cell_type)
  message(sprintf("\n%s:", sid))
  for (nm in names(ct)) message(sprintf("  %-25s %d", nm, ct[nm]))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: Merge annotated objects → Harmony for combined visualization
# Cell type labels come from per-sample annotation — not re-derived here
# ─────────────────────────────────────────────────────────────────────────────

if (file.exists(CACHE_ALL)) {
  message("\nLoading cached merged+harmony object...")
  hl <- readRDS(CACHE_ALL)
} else {
  message("\n=== Merging annotated objects ===")
  hl <- merge(obj_list[[1]],
              y    = obj_list[-1],
              add.cell.ids = names(obj_list))
  hl <- JoinLayers(hl)

  # Re-run HVG + scale + PCA on merged object for Harmony input
  hl <- FindVariableFeatures(hl, selection.method = "mean.var.plot",
                             mean.cutoff = c(0.0125, 3),
                             dispersion.cutoff = c(0.5, Inf),
                             verbose = FALSE)
  hl <- ScaleData(hl, vars.to.regress = c("nCount_RNA", "percent.mt"),
                  verbose = FALSE)
  hl <- RunPCA(hl, npcs = 40, verbose = FALSE)

  # Harmony for batch correction across samples (visualization only)
  hl <- RunHarmony(hl, group.by.vars = "sample_id", reduction.use = "pca",
                   dims.use = 1:40, verbose = TRUE)

  hl <- RunUMAP(hl, reduction = "harmony", dims = 1:31)
  hl <- RunTSNE(hl, reduction = "harmony", dims = 1:31)

  saveRDS(hl, CACHE_ALL)
  message(sprintf("Merged object saved: %s", CACHE_ALL))
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 7: Summary + figures (all 17 samples)
# ─────────────────────────────────────────────────────────────────────────────

message("\n=== All-sample summary ===")
message(sprintf("Total cells: %d", ncol(hl)))
print(table(hl$condition))
print(sort(table(hl$cell_type), decreasing = TRUE))

Idents(hl) <- "cell_type"

pdf(file.path(OUT_DIR, "all_umap_condition.pdf"), width = 8, height = 6)
print(DimPlot(hl, reduction = "umap", group.by = "condition", pt.size = 0.3) +
      ggtitle("All samples - condition"))
dev.off()

pdf(file.path(OUT_DIR, "all_umap_celltype.pdf"), width = 11, height = 7)
print(DimPlot(hl, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.3) +
      ggtitle("All samples - cell type (per-sample annotation)"))
dev.off()

# ─────────────────────────────────────────────────────────────────────────────
# STEP 8: Subset to Donor + IPF
# ─────────────────────────────────────────────────────────────────────────────

message("\n=== Subsetting to Donor + IPF ===")
hl_ipf <- subset(hl, subset = condition %in% c("Donor", "IPF"))
hl_ipf$condition <- factor(hl_ipf$condition, levels = c("Donor", "IPF"))
hl_ipf <- RunUMAP(hl_ipf, reduction = "harmony", dims = 1:31)
hl_ipf <- RunTSNE(hl_ipf, reduction = "harmony", dims = 1:31)

message(sprintf("IPF+Donor cells: %d", ncol(hl_ipf)))
print(table(hl_ipf$condition))
print(sort(table(hl_ipf$cell_type), decreasing = TRUE))

saveRDS(hl_ipf, CACHE_IPF)
message(sprintf("IPF+Donor object saved: %s", CACHE_IPF))

pdf(file.path(OUT_DIR, "ipf_umap_condition.pdf"), width = 7, height = 6)
print(DimPlot(hl_ipf, reduction = "umap", group.by = "condition",
              cols = c("Donor" = "#2166AC", "IPF" = "#D6604D"), pt.size = 0.4) +
      ggtitle("Donor vs IPF"))
dev.off()

pdf(file.path(OUT_DIR, "ipf_umap_celltype.pdf"), width = 11, height = 7)
print(DimPlot(hl_ipf, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.4) +
      ggtitle("Cell types - Donor + IPF (per-sample annotation)"))
dev.off()

pdf(file.path(OUT_DIR, "ipf_fibro_markers.pdf"), width = 12, height = 5)
print(FeaturePlot(hl_ipf, features = c("DCN", "COL1A1", "COL3A1", "ACTA2", "SFTPC", "SCGB3A2"),
                  reduction = "umap", ncol = 6, pt.size = 0.2))
dev.off()

message("\n=== Done ===")
message(sprintf("Outputs in: %s", OUT_DIR))
