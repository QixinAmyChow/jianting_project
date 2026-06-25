# =============================================================================
# IPF Myofibroblast / Fibroblast Analysis
# Dataset : Habermann et al. 2020, Science Advances — GSE135893
#           114,396 cells | 12 IPF + 10 Control | Vanderbilt + NTI
#
# Output  : ipf/figures_Habermann2020/
# Cache   : ipf/Habermann2020/cache_*.rds
#
# Figures produced (PDF + TIFF 300 dpi):
#   fig1_all_fibro   — Myofibroblasts + Fibroblasts + PLIN2+ Fib + HAS1 Fib
#   fig2_fibroblasts — Fibroblasts only
#   fig3_myofibro    — Myofibroblasts only
#   Each figure: cell-level + patient-level, both genes, p-value annotated
#
# Version log:
#   v1 (2026-06-24) — initial build: load, subset, VlnPlot, DotPlot, DEG
#   v2 (2026-06-24) — MFAP5/TIMP1; Nature aesthetics; cell + patient level;
#                     p-value annotations; PDF + TIFF output; 3 figures
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# LIBRARIES
# ─────────────────────────────────────────────────────────────────────────────

library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(ggsignif)


# ─────────────────────────────────────────────────────────────────────────────
# PARAMETERS  ← edit here
# ─────────────────────────────────────────────────────────────────────────────

GENES_OI   <- c("MFAP5", "TIMP1")

# Cell type groupings for the three figures
CT_ALL <- c("Myofibroblasts", "Fibroblasts", "PLIN2+ Fibroblasts",
            "HAS1 High Fibroblasts")               # fig 1
CT_FIB <- "Fibroblasts"                            # fig 2
CT_MYO <- "Myofibroblasts"                         # fig 3

DX_COLORS  <- c("Control" = "#2166AC", "IPF" = "#D6604D")   # Nature-style blue/red
CONDITIONS <- c("Control", "IPF")

DATA_DIR   <- "ipf/Habermann2020"
OUT_DIR    <- "ipf/figures_Habermann2020"
TIFF_RES   <- 300   # dpi for TIFF


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# load_geo_sparse()
#   Loads a standard GEO-style directory (matrix.mtx.gz, genes.tsv.gz,
#   barcodes.tsv.gz, *metadata*.csv.gz) into a Seurat object.
# -----------------------------------------------------------------------------
load_geo_sparse <- function(dir,
                             matrix_f   = "GSE135893_matrix.mtx.gz",
                             genes_f    = "GSE135893_genes.tsv.gz",
                             barcodes_f = "GSE135893_barcodes.tsv.gz",
                             metadata_f = "GSE135893_IPF_metadata.csv.gz",
                             project    = "GEO_dataset") {

  cat("Reading matrix:", file.path(dir, matrix_f), "\n")
  mat      <- readMM(file.path(dir, matrix_f))
  genes    <- read.table(file.path(dir, genes_f),    header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.table(file.path(dir, barcodes_f), header = FALSE, stringsAsFactors = FALSE)

  rownames(mat) <- make.unique(genes$V2)   # V2 = gene symbol; V1 = Ensembl ID
  colnames(mat) <- barcodes$V1

  cat("Reading metadata:", file.path(dir, metadata_f), "\n")
  meta   <- read.csv(file.path(dir, metadata_f), row.names = 1)
  shared <- intersect(colnames(mat), rownames(meta))
  cat("Cells — matrix:", ncol(mat), "| metadata:", nrow(meta),
      "| shared:", length(shared), "\n")

  mat  <- mat[, shared]
  meta <- meta[shared, ]

  CreateSeuratObject(counts = mat, meta.data = meta, project = project)
}


# -----------------------------------------------------------------------------
# inspect_seurat()
#   Structured summary: dims, metadata cols, cell type counts, group breakdown.
# -----------------------------------------------------------------------------
inspect_seurat <- function(obj,
                            celltype_col = "celltype",
                            group_col    = "Diagnosis") {

  cat("═══════════════════════════════════════\n")
  cat("Dims    :", nrow(obj), "genes ×", ncol(obj), "cells\n")
  cat("Assays  :", paste(Assays(obj), collapse = ", "), "\n")
  cat("Meta    :", paste(colnames(obj@meta.data), collapse = ", "), "\n\n")

  if (celltype_col %in% colnames(obj@meta.data)) {
    cat("── Cell types ──\n")
    print(sort(table(obj@meta.data[[celltype_col]]), decreasing = TRUE))
    cat("\n")
  }
  if (group_col %in% colnames(obj@meta.data)) {
    cat("── Groups ──\n")
    print(table(obj@meta.data[[group_col]]))
    if (celltype_col %in% colnames(obj@meta.data)) {
      cat("\n── Celltype × Group ──\n")
      print(table(obj@meta.data[[celltype_col]], obj@meta.data[[group_col]]))
    }
    cat("\n")
  }
  cat("═══════════════════════════════════════\n")
}


# -----------------------------------------------------------------------------
# subset_by_celltype()
#   Subsets to one or more cell types and specified conditions.
# -----------------------------------------------------------------------------
subset_by_celltype <- function(obj,
                                celltypes,
                                conditions,
                                celltype_col = "celltype",
                                group_col    = "Diagnosis") {

  keep <- obj@meta.data[[celltype_col]] %in% celltypes &
          obj@meta.data[[group_col]]    %in% conditions
  sub  <- subset(obj, cells = rownames(obj@meta.data)[keep])
  sub@meta.data[[group_col]] <- factor(sub@meta.data[[group_col]],
                                        levels = conditions)
  cat("Subset:", paste(celltypes, collapse = " + "), "—",
      paste(paste0(conditions, ": ",
                   table(sub@meta.data[[group_col]])), collapse = " | "), "\n")
  sub
}


# -----------------------------------------------------------------------------
# check_genes()
#   Returns genes present in object; prints missing ones.
# -----------------------------------------------------------------------------
check_genes <- function(obj, genes) {
  found   <- genes[genes %in% rownames(obj)]
  missing <- genes[!genes %in% rownames(obj)]
  if (length(missing)) cat("Missing:", paste(missing, collapse = ", "), "\n")
  cat("Using  :", paste(found,   collapse = ", "), "\n")
  found
}


# -----------------------------------------------------------------------------
# theme_nature()
#   Minimal theme matching Nature journal figure standards.
#   base_size ~7 pt is standard for Nature figure panels.
# -----------------------------------------------------------------------------
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
  theme(
    text             = element_text(family = "sans", size = base_size,
                                    color  = "black"),
    axis.text        = element_text(size = base_size, color = "black"),
    axis.title       = element_text(size = base_size + 1),
    axis.line        = element_line(linewidth = 0.35, color = "black"),
    axis.ticks       = element_line(linewidth = 0.35, color = "black"),
    axis.ticks.length = unit(1.5, "mm"),
    legend.text      = element_text(size = base_size),
    legend.title     = element_text(size = base_size, face = "bold"),
    legend.key.size  = unit(3,   "mm"),
    legend.position  = "right",
    plot.title       = element_text(size = base_size + 1, face = "bold",
                                    hjust = 0.5),
    plot.subtitle    = element_text(size = base_size, hjust = 0.5,
                                    color = "grey40"),
    strip.background = element_blank(),
    strip.text       = element_text(size = base_size, face = "bold"),
    panel.spacing    = unit(2, "mm")
  )
}


# -----------------------------------------------------------------------------
# save_figure()
#   Saves a ggplot/patchwork object as both PDF and TIFF.
#   TIFF at 300 dpi, PDF vector.
# -----------------------------------------------------------------------------
save_figure <- function(plot, prefix, width_mm, height_mm) {

  w_in <- width_mm  / 25.4
  h_in <- height_mm / 25.4

  ggsave(paste0(prefix, ".pdf"),  plot, width = w_in, height = h_in,
         device = "pdf",  units = "in")
  ggsave(paste0(prefix, ".tiff"), plot, width = w_in, height = h_in,
         device = "tiff", units = "in", dpi = TIFF_RES, compression = "lzw")

  cat("Saved:", basename(prefix), ".pdf + .tiff\n")
}


# -----------------------------------------------------------------------------
# make_cell_level_plot()
#   Violin + boxplot for one gene, Control vs IPF at cell level.
#   Wilcoxon p-value annotation via ggpubr.
# -----------------------------------------------------------------------------
make_cell_level_plot <- function(obj, gene, group_col = "Diagnosis",
                                  colors = DX_COLORS, label = "") {

  df <- data.frame(
    expr  = FetchData(obj, vars = gene)[, 1],
    group = obj@meta.data[[group_col]]
  )

  wt  <- wilcox.test(expr ~ group, data = df, exact = FALSE)
  p   <- wt$p.value
  # Bonferroni correction for number of genes tested
  p_adj <- min(p * length(GENES_OI), 1)
  p_lab <- ifelse(p_adj < 0.001,
                  sprintf("p[adj] == '%.2e'", p_adj),
                  sprintf("p[adj] == '%.3f'", p_adj))

  ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7,
                linewidth = 0.3) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA,
                 linewidth = 0.35) +
    geom_signif(
      comparisons  = list(c("Control", "IPF")),
      annotations  = ifelse(p_adj < 0.05,
                            ifelse(p_adj < 0.001, "***",
                                   ifelse(p_adj < 0.01, "**", "*")),
                            "ns"),
      tip_length   = 0.01,
      textsize     = 2.2,
      size         = 0.3
    ) +
    annotate("text", x = 1.5, y = max(df$expr) * 1.35,
             label = p_lab, parse = TRUE, size = 1.8, color = "grey30") +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(x = NULL, y = "log-norm expression",
         title = gene,
         subtitle = paste0("Cell level  n=", nrow(df))) +
    theme_nature() +
    theme(legend.position = "none")
}


# -----------------------------------------------------------------------------
# make_patient_level_plot()
#   Pseudobulk: mean expression per patient per gene.
#   Dots = individual patients, box = IQR, Wilcoxon test.
# -----------------------------------------------------------------------------
make_patient_level_plot <- function(obj, gene, group_col = "Diagnosis",
                                     patient_col = "Sample_Name",
                                     colors = DX_COLORS) {

  pb <- data.frame(
    expr    = FetchData(obj, vars = gene)[, 1],
    group   = obj@meta.data[[group_col]],
    patient = obj@meta.data[[patient_col]]
  ) %>%
    group_by(patient, group) %>%
    summarise(mean_expr = mean(expr), .groups = "drop")

  wt    <- wilcox.test(mean_expr ~ group, data = pb, exact = FALSE)
  p     <- wt$p.value
  p_adj <- min(p * length(GENES_OI), 1)
  p_lab <- ifelse(p_adj < 0.001,
                  sprintf("p[adj] == '%.2e'", p_adj),
                  sprintf("p[adj] == '%.3f'", p_adj))

  n_ctrl <- sum(pb$group == "Control")
  n_ipf  <- sum(pb$group == "IPF")

  ggplot(pb, aes(x = group, y = mean_expr, fill = group, color = group)) +
    geom_boxplot(width = 0.45, alpha = 0.4, outlier.shape = NA,
                 linewidth = 0.35) +
    geom_jitter(width = 0.1, size = 1.2, alpha = 0.85, shape = 21,
                stroke = 0.3, color = "white") +
    geom_signif(
      comparisons  = list(c("Control", "IPF")),
      annotations  = ifelse(p_adj < 0.05,
                            ifelse(p_adj < 0.001, "***",
                                   ifelse(p_adj < 0.01, "**", "*")),
                            "ns"),
      tip_length   = 0.015,
      textsize     = 2.2,
      size         = 0.3,
      color        = "black"
    ) +
    annotate("text", x = 1.5, y = max(pb$mean_expr) * 1.4,
             label = p_lab, parse = TRUE, size = 1.8, color = "grey30") +
    scale_fill_manual(values  = colors) +
    scale_color_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    labs(x = NULL, y = "mean log-norm expression",
         title = gene,
         subtitle = sprintf("Patient level  ctrl=%d  IPF=%d", n_ctrl, n_ipf)) +
    theme_nature() +
    theme(legend.position = "none")
}


# -----------------------------------------------------------------------------
# make_figure()
#   Assembles a 2-row × n-gene patchwork figure for one cell type subset.
#   Row 1 = cell level, Row 2 = patient level.
#   Returns a patchwork object.
# -----------------------------------------------------------------------------
make_figure <- function(obj, genes, title,
                         group_col   = "Diagnosis",
                         patient_col = "Sample_Name",
                         colors      = DX_COLORS) {

  cell_panels    <- lapply(genes, make_cell_level_plot,
                            obj = obj, group_col = group_col, colors = colors)
  patient_panels <- lapply(genes, make_patient_level_plot,
                            obj = obj, group_col = group_col,
                            patient_col = patient_col, colors = colors)

  all_panels <- c(cell_panels, patient_panels)
  n <- length(genes)

  wrap_plots(all_panels, nrow = 2, ncol = n) +
    plot_annotation(
      title   = title,
      theme   = theme(
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5,
                                   family = "sans")
      )
    )
}


# =============================================================================
# SECTION 1 — LOAD & CACHE SUBSETS
# =============================================================================

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cache_all <- file.path(DATA_DIR, "cache_all_fibro.rds")
cache_fib <- file.path(DATA_DIR, "cache_fibroblasts.rds")
cache_myo <- file.path(DATA_DIR, "cache_myofibroblasts.rds")

if (!all(file.exists(cache_all, cache_fib, cache_myo))) {

  cat("\n[1/4] Loading full dataset...\n")
  seu <- load_geo_sparse(DATA_DIR, project = "GSE135893_ILD")
  inspect_seurat(seu)

  cat("\n[2/4] Normalizing...\n")
  seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                       scale.factor = 10000)

  cat("\n[3/4] Creating subsets...\n")
  sub_all <- subset_by_celltype(seu, CT_ALL, CONDITIONS)
  sub_fib <- subset_by_celltype(seu, CT_FIB, CONDITIONS)
  sub_myo <- subset_by_celltype(seu, CT_MYO, CONDITIONS)

  cat("\n[4/4] Caching subsets...\n")
  saveRDS(sub_all, cache_all)
  saveRDS(sub_fib, cache_fib)
  saveRDS(sub_myo, cache_myo)
  rm(seu); gc()

} else {

  cat("Loading cached subsets...\n")
  sub_all <- readRDS(cache_all)
  sub_fib <- readRDS(cache_fib)
  sub_myo <- readRDS(cache_myo)

}

# Verify genes
genes_found <- check_genes(sub_myo, GENES_OI)


# =============================================================================
# SECTION 2 — FIGURES  (cell level + patient level, p-value annotated)
# =============================================================================
# Each figure: 2 rows (cell / patient) × n genes columns
# Saved as PDF (vector) + TIFF 300 dpi

# ── v2: Figure 1 — All fibroblast clusters ───────────────────────────────────

fig1 <- make_figure(
  obj    = sub_all,
  genes  = genes_found,
  title  = "All fibroblast clusters — IPF vs Control (Habermann 2020)"
)
save_figure(fig1,
            prefix     = file.path(OUT_DIR, "v2_fig1_all_fibro_clusters"),
            width_mm   = 85 * length(genes_found),
            height_mm  = 120)


# ── v2: Figure 2 — Fibroblasts only ──────────────────────────────────────────

fig2 <- make_figure(
  obj    = sub_fib,
  genes  = genes_found,
  title  = "Fibroblasts — IPF vs Control (Habermann 2020)"
)
save_figure(fig2,
            prefix     = file.path(OUT_DIR, "v2_fig2_fibroblasts_only"),
            width_mm   = 85 * length(genes_found),
            height_mm  = 120)


# ── v2: Figure 3 — Myofibroblasts only ───────────────────────────────────────

fig3 <- make_figure(
  obj    = sub_myo,
  genes  = genes_found,
  title  = "Myofibroblasts — IPF vs Control (Habermann 2020)"
)
save_figure(fig3,
            prefix     = file.path(OUT_DIR, "v2_fig3_myofibroblasts_only"),
            width_mm   = 85 * length(genes_found),
            height_mm  = 120)


# =============================================================================
# SECTION 3 — DEG  (kept from v1 for reference)
# =============================================================================

# ── v1: Wilcoxon — Myofibroblasts ────────────────────────────────────────────

Idents(sub_myo) <- "Diagnosis"
deg_myo <- FindMarkers(
  sub_myo,
  ident.1         = "IPF",
  ident.2         = "Control",
  test.use        = "wilcox",
  min.pct         = 0.1,
  logfc.threshold = 0.1
) %>% arrange(desc(avg_log2FC))

write.csv(deg_myo, file.path(OUT_DIR, "v1_myofib_IPF_vs_ctrl_DEG.csv"))

cat("\nTop 20 upregulated in IPF myofibroblasts (padj < 0.05):\n")
print(head(filter(deg_myo, p_val_adj < 0.05), 20))

cat("\nDone. All outputs in:", OUT_DIR, "\n")

