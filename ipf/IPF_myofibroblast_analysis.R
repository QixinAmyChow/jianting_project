# =============================================================================
# IPF Fibroblast / Myofibroblast Gene Expression Analysis
# Dataset : Habermann et al. 2020, Science Advances — GSE135893
#           114,396 cells | 12 IPF + 10 Control | Vanderbilt + NTI
#
# Output  : ipf/figures/Habermann2020/
# Cache   : ipf/Habermann2020/cache_*.rds
#
# Figures (PDF only):
#   v2_violin_Myo_Fib.pdf          cell-level violin, Myo + Fib, faceted by gene
#   v2_boxplot_Myo_Fib.pdf         cell-level boxplot, Myo + Fib, faceted by gene
#   v2_violin_allcelltypes.pdf     violin, all cell types with >=20 cells/group
#   v2_boxplot_allcelltypes.pdf    boxplot, same
#
#   X-axis: cell type; fill: Diagnosis (Control / IPF) — dodged pairs
#   Facets: one panel per gene
#   Stats:  Wilcoxon BH-adjusted p-value per gene × cell type
#
# Version log:
#   v1 (2026-06-24) — load, subset, Seurat VlnPlot + DotPlot + DEG
#   v2 (2026-06-24) — MFAP5/TIMP1; Nature aesthetics; cell + patient level;
#                     p-value; PDF + TIFF; 3 figures (superseded)
#   v3 (2026-06-24) — pure ggplot; faceted by gene; dodged Ctrl/IPF pairs;
#                     violin and boxplot separate; all-celltypes figure;
#                     study subdirectory output; no CSV
#   v4 (2026-06-24) — each gene own PDF; violin = violin only (no inner box);
#                     per-gene figure sizing; aesthetic refinements
#   v5 (2026-06-24) — restore v2 aesthetic (violin + inner boxplot); add
#                     sample-count annotation (ctrl / ipf) in upper-left corner
#   v6 (2026-06-24) — add jitter scatter points (violin → dots → inner box);
#   v7 (2026-06-24) — violin: white inner box, no scatter (v2 clean aesthetic);
#                     boxplot: patient-level means + per-sample dots;
#                     p-value label: stars + exact padj numeric value
#   v8 (2026-06-24) — exact v2 geoms: alpha=0.7, inner box width=0.12,
#                     patient box alpha=0.4, jitter shape=21 white border,
#                     geom_signif + parse annotate for padj; correct colors
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# LIBRARIES
# ─────────────────────────────────────────────────────────────────────────────

library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(rstatix)
library(ggpubr)
library(ggsignif)


# ─────────────────────────────────────────────────────────────────────────────
# PARAMETERS  ← edit here
# ─────────────────────────────────────────────────────────────────────────────

GENES_OI   <- c("MFAP5", "TIMP1")

# Cell types for the specific-pair figure
CT_SPECIFIC <- c("Myofibroblasts", "Fibroblasts")

# Minimum cells per group (Control or IPF) for a cell type to appear
MIN_CELLS  <- 20

CONDITIONS <- c("Control", "IPF")
DX_COLORS  <- c("Control" = "#2166AC", "IPF" = "#D6604D")
DODGE_W    <- 0.8

DATA_DIR   <- "ipf/Habermann2020"
OUT_DIR    <- "ipf/figures/Habermann2020"


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# load_geo_sparse()
#   Loads a GEO-style directory (matrix.mtx.gz, genes.tsv.gz, barcodes.tsv.gz,
#   metadata.csv.gz) into a Seurat object.
# -----------------------------------------------------------------------------
load_geo_sparse <- function(dir,
                             matrix_f   = "GSE135893_matrix.mtx.gz",
                             genes_f    = "GSE135893_genes.tsv.gz",
                             barcodes_f = "GSE135893_barcodes.tsv.gz",
                             metadata_f = "GSE135893_IPF_metadata.csv.gz",
                             project    = "GEO_dataset") {

  cat("Reading matrix:", file.path(dir, matrix_f), "\n")
  mat      <- readMM(file.path(dir, matrix_f))
  genes    <- read.table(file.path(dir, genes_f),    header = FALSE,
                         stringsAsFactors = FALSE)
  barcodes <- read.table(file.path(dir, barcodes_f), header = FALSE,
                         stringsAsFactors = FALSE)

  # Use V2 (symbol) if two columns exist, otherwise V1
  gene_col      <- if (ncol(genes) >= 2) genes$V2 else genes$V1
  rownames(mat) <- make.unique(as.character(gene_col))
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
#   Structured summary: dims, metadata, cell type counts, group breakdown.
# -----------------------------------------------------------------------------
inspect_seurat <- function(obj,
                            celltype_col = "celltype",
                            group_col    = "Diagnosis") {
  cat("═══════════════════════════════════════\n")
  cat("Dims    :", nrow(obj), "genes x", ncol(obj), "cells\n")
  cat("Assays  :", paste(Assays(obj), collapse = ", "), "\n")
  cat("Meta    :", paste(colnames(obj@meta.data), collapse = ", "), "\n\n")
  if (celltype_col %in% colnames(obj@meta.data)) {
    cat("-- Cell types --\n")
    print(sort(table(obj@meta.data[[celltype_col]]), decreasing = TRUE))
    cat("\n")
  }
  if (group_col %in% colnames(obj@meta.data)) {
    cat("-- Groups --\n")
    print(table(obj@meta.data[[group_col]]))
    if (celltype_col %in% colnames(obj@meta.data)) {
      cat("\n-- Celltype x Group --\n")
      print(table(obj@meta.data[[celltype_col]],
                  obj@meta.data[[group_col]]))
    }
    cat("\n")
  }
  cat("═══════════════════════════════════════\n")
}


# -----------------------------------------------------------------------------
# subset_by_celltype()
#   Subsets to specified cell types and conditions; re-levels condition factor.
# -----------------------------------------------------------------------------
subset_by_celltype <- function(obj, celltypes, conditions,
                                celltype_col = "celltype",
                                group_col    = "Diagnosis") {
  keep <- obj@meta.data[[celltype_col]] %in% celltypes &
          obj@meta.data[[group_col]]    %in% conditions
  sub  <- subset(obj, cells = rownames(obj@meta.data)[keep])
  sub@meta.data[[group_col]] <- factor(sub@meta.data[[group_col]],
                                        levels = conditions)
  cat("Subset:", paste(celltypes, collapse = " + "), "--",
      paste(paste0(conditions, ":",
                   table(sub@meta.data[[group_col]])), collapse = "  "), "\n")
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
  cat("Using  :", paste(found, collapse = ", "), "\n")
  found
}


# -----------------------------------------------------------------------------
# theme_nature()
#   Minimal theme following Nature journal figure conventions.
# -----------------------------------------------------------------------------
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
  theme(
    text              = element_text(family = "sans", size = base_size,
                                     color = "black"),
    axis.text         = element_text(size = base_size, color = "black"),
    axis.text.x       = element_text(angle = 35, hjust = 1),
    axis.title        = element_text(size = base_size + 1),
    axis.line         = element_line(linewidth = 0.35, color = "black"),
    axis.ticks        = element_line(linewidth = 0.35, color = "black"),
    axis.ticks.length = unit(1.5, "mm"),
    legend.text       = element_text(size = base_size),
    legend.title      = element_text(size = base_size, face = "bold"),
    legend.key.size   = unit(3, "mm"),
    legend.position   = "right",
    plot.title        = element_text(size = base_size + 1, face = "bold",
                                     hjust = 0.5),
    strip.background  = element_blank(),
    strip.text        = element_text(size = base_size + 1, face = "bold.italic"),
    panel.spacing     = unit(3, "mm")
  )
}


# -----------------------------------------------------------------------------
# get_sample_label()
#   Returns "Control: N    IPF: M" string from a Seurat object's metadata.
#   Tries common sample-ID column names; falls back to NULL if none found.
# -----------------------------------------------------------------------------
get_sample_label <- function(obj, group_col = "Diagnosis") {
  sample_cols <- c("Sample_Name", "sample", "orig.ident",
                   "Patient", "donor_id", "subject")
  scol <- intersect(sample_cols, colnames(obj@meta.data))[1]
  if (is.na(scol)) return(NULL)
  dx  <- obj@meta.data[[group_col]]
  ids <- obj@meta.data[[scol]]
  n_ctrl <- length(unique(ids[dx == "Control"]))
  n_ipf  <- length(unique(ids[dx == "IPF"]))
  paste0("n(Control)=", n_ctrl, "   n(IPF)=", n_ipf)
}


# -----------------------------------------------------------------------------
# save_pdf()
#   Saves a ggplot/patchwork as PDF at specified mm dimensions.
# -----------------------------------------------------------------------------
save_pdf <- function(plot, path, width_mm, height_mm) {
  ggsave(path, plot,
         width  = width_mm / 25.4,
         height = height_mm / 25.4,
         device = "pdf", units = "in")
  cat("Saved:", basename(path), "\n")
}


# -----------------------------------------------------------------------------
# build_long_df()
#   Extracts gene expression + metadata from a Seurat object into a long-
#   format data frame suitable for ggplot.
# -----------------------------------------------------------------------------
build_long_df <- function(obj, genes,
                           celltype_col = "celltype",
                           group_col    = "Diagnosis") {
  expr <- as.data.frame(t(as.matrix(
    GetAssayData(obj, layer = "data")[genes, , drop = FALSE]
  )))
  expr[[celltype_col]] <- obj@meta.data[[celltype_col]]
  expr[[group_col]]    <- obj@meta.data[[group_col]]

  expr %>%
    pivot_longer(cols = all_of(genes),
                 names_to  = "gene",
                 values_to = "expr") %>%
    mutate(
      gene      = factor(gene, levels = genes),
      celltype  = factor(.data[[celltype_col]]),
      Diagnosis = factor(.data[[group_col]], levels = CONDITIONS)
    )
}


# -----------------------------------------------------------------------------
# build_patient_df()
#   Averages cell-level expression per patient × gene × celltype.
#   Returns same long format as build_long_df() but one row per patient group.
# -----------------------------------------------------------------------------
build_patient_df <- function(obj, genes,
                              celltype_col = "celltype",
                              group_col    = "Diagnosis",
                              sample_col   = "Sample_Name") {

  long <- build_long_df(obj, genes, celltype_col, group_col)

  # Cell order in long df after pivot: for n_genes genes, each cell repeats n_genes times
  long$sample <- rep(obj@meta.data[[sample_col]], each = length(genes))

  long %>%
    group_by(gene, celltype, Diagnosis, sample) %>%
    summarise(expr = mean(expr), .groups = "drop")
}


# -----------------------------------------------------------------------------
# compute_stats()
#   Wilcoxon test (Control vs IPF) per gene x cell type, BH adjusted.
#   Returns stat table with xy positions for stat_pvalue_manual().
# -----------------------------------------------------------------------------
compute_stats <- function(df, min_nonzero = 3) {

  df_testable <- df %>%
    group_by(gene, celltype, Diagnosis) %>%
    filter(sum(expr > 0) >= min_nonzero) %>%
    group_by(gene, celltype) %>%
    filter(n_distinct(Diagnosis) == 2) %>%
    ungroup()

  stats <- df_testable %>%
    group_by(gene, celltype) %>%
    wilcox_test(expr ~ Diagnosis, exact = FALSE) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    mutate(label_full = paste0(p.adj.signif, "\n",
                               "padj=", formatC(p.adj, format = "g", digits = 3)))

  y_caps <- df %>%
    group_by(gene, celltype) %>%
    summarise(y_cap = quantile(expr, 0.99), .groups = "drop")

  # x = "Diagnosis": each celltype facet has one bracket (Control → IPF)
  stats %>%
    add_xy_position(x = "Diagnosis") %>%
    left_join(y_caps, by = c("gene", "celltype")) %>%
    mutate(y.position = y_cap * 1.12)
}


# -----------------------------------------------------------------------------
# plot_violin()  /  plot_boxplot()
#   Single-gene figures, faceted by cell type (x = Diagnosis within each facet).
#   Exact v2 geom parameters and p-value style.
#
#   violin  : v2 — geom_violin(alpha=0.7) + white inner geom_boxplot(width=0.12)
#             + geom_signif bracket + annotate padj text (parse=TRUE)
#   boxplot : v2 — geom_boxplot(alpha=0.4) + geom_jitter(shape=21, white outline)
#             + same p-value style
#
#   Per-facet p-values computed from pre-built stats table.
# -----------------------------------------------------------------------------

# Helper: build per-facet annotation df for the padj text label
.padj_labels <- function(stats, df) {
  stats %>%
    left_join(
      df %>%
        group_by(celltype) %>%
        summarise(y_max = max(expr), .groups = "drop"),
      by = "celltype"
    ) %>%
    mutate(
      x     = 1.5,
      y_ann = y_max * 1.38,
      label = ifelse(p.adj < 0.001,
                     sprintf("p[adj] == '%.2e'", p.adj),
                     sprintf("p[adj] == '%.3f'", p.adj))
    )
}


plot_violin <- function(df, stats, title = "", sample_label = NULL) {

  ann <- .padj_labels(stats, df)

  p <- ggplot(df, aes(x = Diagnosis, y = expr, fill = Diagnosis)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7, linewidth = 0.3) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA,
                 linewidth = 0.35) +
    # bracket + stars (stat_pvalue_manual respects facets via celltype column)
    stat_pvalue_manual(stats, label = "p.adj.signif",
                       tip.length = 0.01, bracket.size = 0.3,
                       size = 2.2, hide.ns = FALSE) +
    # exact padj value as parsed text (v2 style)
    geom_text(data = ann, aes(x = x, y = y_ann, label = label),
              inherit.aes = FALSE, parse = TRUE, size = 1.8, color = "grey30") +
    scale_fill_manual(values = DX_COLORS, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    facet_wrap(~ celltype, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = "log-norm expression", title = title,
         subtitle = paste0("Cell level  n=", nrow(df))) +
    theme_nature() +
    theme(legend.position  = "top",
          legend.direction = "horizontal",
          legend.key.size  = unit(2.5, "mm"),
          legend.spacing.x = unit(1,   "mm"),
          strip.text       = element_text(size = 7, face = "bold"))

  if (!is.null(sample_label)) {
    p <- p + labs(tag = sample_label) +
      theme(plot.tag          = element_text(size = 6, color = "grey40",
                                              hjust = 0, vjust = 1),
            plot.tag.position = "topleft")
  }
  p
}


plot_boxplot <- function(df_patient, stats, title = "", sample_label = NULL) {

  ann <- .padj_labels(stats, df_patient)
  n_ctrl <- length(unique(df_patient$sample[df_patient$Diagnosis == "Control"]))
  n_ipf  <- length(unique(df_patient$sample[df_patient$Diagnosis == "IPF"]))

  p <- ggplot(df_patient, aes(x = Diagnosis, y = expr,
                               fill = Diagnosis, color = Diagnosis)) +
    geom_boxplot(width = 0.45, alpha = 0.4, outlier.shape = NA,
                 linewidth = 0.35) +
    # v2: shape=21, white border stroke
    geom_jitter(width = 0.10, size = 1.2, alpha = 0.85,
                shape = 21, stroke = 0.3, color = "white") +
    stat_pvalue_manual(stats, label = "p.adj.signif",
                       tip.length = 0.015, bracket.size = 0.3,
                       size = 2.2, hide.ns = FALSE, color = "black") +
    geom_text(data = ann, aes(x = x, y = y_ann, label = label),
              inherit.aes = FALSE, parse = TRUE, size = 1.8, color = "grey30") +
    scale_fill_manual(values  = DX_COLORS, name = NULL) +
    scale_color_manual(values = DX_COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    facet_wrap(~ celltype, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = "mean log-norm expression", title = title,
         subtitle = sprintf("Patient level  ctrl=%d  IPF=%d", n_ctrl, n_ipf)) +
    theme_nature() +
    theme(legend.position  = "top",
          legend.direction = "horizontal",
          legend.key.size  = unit(2.5, "mm"),
          legend.spacing.x = unit(1,   "mm"),
          strip.text       = element_text(size = 7, face = "bold"))

  if (!is.null(sample_label)) {
    p <- p + labs(tag = sample_label) +
      theme(plot.tag          = element_text(size = 6, color = "grey40",
                                              hjust = 0, vjust = 1),
            plot.tag.position = "topleft")
  }
  p
}


# =============================================================================
# SECTION 1 — LOAD & CACHE SUBSETS
# =============================================================================

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cache_specific  <- file.path(DATA_DIR, "cache_Myo_Fib.rds")
cache_allcells  <- file.path(DATA_DIR, "cache_allcells.rds")

if (!all(file.exists(cache_specific, cache_allcells))) {

  cat("\n[1/4] Loading full dataset...\n")
  seu <- load_geo_sparse(DATA_DIR, project = "GSE135893_ILD")
  inspect_seurat(seu)

  cat("\n[2/4] Normalizing...\n")
  seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                       scale.factor = 10000)

  cat("\n[3/4] Building subsets...\n")

  # Subset 1: Myofibroblasts + Fibroblasts only
  sub_specific <- subset_by_celltype(seu, CT_SPECIFIC, CONDITIONS)
  saveRDS(sub_specific, cache_specific)

  # Subset 2: all cell types that have >= MIN_CELLS in both Control and IPF
  ct_counts <- table(seu@meta.data$celltype, seu@meta.data$Diagnosis)
  keep_ct   <- rownames(ct_counts)[
    ct_counts[, "Control"] >= MIN_CELLS & ct_counts[, "IPF"] >= MIN_CELLS
  ]
  cat("Cell types passing min-cell filter:", length(keep_ct), "\n")
  sub_all <- subset_by_celltype(seu, keep_ct, CONDITIONS)
  saveRDS(sub_all, cache_allcells)

  cat("\n[4/4] Caches saved.\n")
  rm(seu); gc()

} else {

  cat("Loading cached subsets...\n")
  sub_specific <- readRDS(cache_specific)
  sub_allcells <- readRDS(cache_allcells)

}

# Load if not already in memory (second run path above assigns sub_all not sub_allcells)
if (!exists("sub_allcells")) sub_allcells <- readRDS(cache_allcells)

genes_found <- check_genes(sub_specific, GENES_OI)


# =============================================================================
# SECTION 2 — BUILD DATA FRAMES + STATS
# =============================================================================

# Cell-level (for violin)
df_specific <- build_long_df(sub_specific, genes_found)
df_allcells <- build_long_df(sub_allcells, genes_found)

stats_specific <- compute_stats(df_specific)
stats_allcells <- compute_stats(df_allcells)

# Patient-level (for boxplot) — one mean per sample × celltype × gene
pt_specific <- build_patient_df(sub_specific, genes_found)
pt_allcells <- build_patient_df(sub_allcells, genes_found)

stats_pt_specific <- compute_stats(pt_specific, min_nonzero = 1)
stats_pt_allcells <- compute_stats(pt_allcells, min_nonzero = 1)


# =============================================================================
# SECTION 3 — FIGURES
# =============================================================================
# One PDF per gene per plot type per subset = 8 PDFs total.
# Violin and boxplot are always SEPARATE files.
# Each gene is its own figure (no multi-gene faceting).

n_ct_spec <- length(CT_SPECIFIC)
n_ct_all  <- nlevels(df_allcells$celltype)

# Width: ~48mm per cell-type facet panel + margins
w_spec <- 48 * n_ct_spec + 20
w_all  <- 40 * n_ct_all  + 20

# Sample-count labels (upper-left corner annotation)
lbl_spec <- get_sample_label(sub_specific)
lbl_all  <- get_sample_label(sub_allcells)

# ── v8: loop over genes ───────────────────────────────────────────────────────

for (gene in genes_found) {

  df_g_spec  <- dplyr::filter(df_specific,    gene == !!gene)
  df_g_all   <- dplyr::filter(df_allcells,    gene == !!gene)
  pt_g_spec  <- dplyr::filter(pt_specific,    gene == !!gene)
  pt_g_all   <- dplyr::filter(pt_allcells,    gene == !!gene)
  st_g_spec  <- dplyr::filter(stats_specific, gene == !!gene)
  st_g_all   <- dplyr::filter(stats_allcells, gene == !!gene)
  sp_g_spec  <- dplyr::filter(stats_pt_specific, gene == !!gene)
  sp_g_all   <- dplyr::filter(stats_pt_allcells, gene == !!gene)

  # ── Myo + Fib: violin (cell-level, white inner box) ────────────────────────
  save_pdf(
    plot_violin(df_g_spec, st_g_spec,
                title        = paste0(gene, "  |  Myo & Fib  |  IPF vs Control"),
                sample_label = lbl_spec),
    file.path(OUT_DIR, paste0("v8_violin_Myo_Fib_", gene, ".pdf")),
    width_mm = w_spec, height_mm = 95
  )

  # ── Myo + Fib: boxplot (patient-level, per-sample dots) ────────────────────
  save_pdf(
    plot_boxplot(pt_g_spec, sp_g_spec,
                 title        = paste0(gene, "  |  Myo & Fib  |  Patient level"),
                 sample_label = lbl_spec),
    file.path(OUT_DIR, paste0("v8_boxplot_Myo_Fib_", gene, ".pdf")),
    width_mm = w_spec, height_mm = 95
  )

  # ── All cell types: violin ─────────────────────────────────────────────────
  save_pdf(
    plot_violin(df_g_all, st_g_all,
                title        = paste0(gene, "  |  All cell types  |  IPF vs Control"),
                sample_label = lbl_all),
    file.path(OUT_DIR, paste0("v8_violin_allcelltypes_", gene, ".pdf")),
    width_mm = w_all, height_mm = 105
  )

  # ── All cell types: boxplot (patient-level) ────────────────────────────────
  save_pdf(
    plot_boxplot(pt_g_all, sp_g_all,
                 title        = paste0(gene, "  |  All cell types  |  Patient level"),
                 sample_label = lbl_all),
    file.path(OUT_DIR, paste0("v8_boxplot_allcelltypes_", gene, ".pdf")),
    width_mm = w_all, height_mm = 105
  )
}

cat("\nDone. All figures in:", OUT_DIR, "\n")
