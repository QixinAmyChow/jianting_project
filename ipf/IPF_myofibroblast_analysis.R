# =============================================================================
# IPF Fibroblast / Myofibroblast Gene Expression Analysis
# Dataset : Habermann et al. 2020, Science Advances — GSE135893
#           114,396 cells (all ILD diagnoses) | Vanderbilt + NTI
#           Subset to Control + IPF only for all figures.
#
# Output  : ipf/figures_Habermann2020/
#   v2_fig1_all_fibro_clusters.pdf   — all 4 fibro subtypes, MFAP5 + TIMP1
#   v2_fig2_fibroblasts_only.pdf     — PI16+ Fibroblasts only
#   v2_fig3_myofibroblasts_only.pdf  — Myofibroblasts only
#
# Cache   : ipf/Habermann2020/cache_*.rds  (built on first run, reused after)
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


# ─────────────────────────────────────────────────────────────────────────────
# PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────

GENES_OI   <- c("MFAP5", "TIMP1")
CONDITIONS <- c("Control", "IPF")
DX_COLORS  <- c("Control" = "#2166AC", "IPF" = "#D6604D")

DATA_DIR   <- "ipf/Habermann2020"
OUT_DIR    <- "ipf/figures_Habermann2020"

ALL_FIBRO_CT <- c("Myofibroblasts", "Fibroblasts",
                  "PLIN2+ Fibroblasts", "HAS1 High Fibroblasts")


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

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


subset_ctrl_ipf <- function(obj, celltypes,
                             celltype_col = "celltype",
                             group_col    = "Diagnosis") {
  keep <- obj@meta.data[[celltype_col]] %in% celltypes &
          obj@meta.data[[group_col]]    %in% CONDITIONS
  sub  <- subset(obj, cells = rownames(obj@meta.data)[keep])
  sub@meta.data[[group_col]] <- factor(sub@meta.data[[group_col]],
                                        levels = CONDITIONS)
  cat("Subset:", paste(celltypes, collapse = " + "), "--",
      paste(paste0(CONDITIONS, ":",
                   table(sub@meta.data[[group_col]])), collapse = "  "), "\n")
  sub
}


check_genes <- function(obj, genes) {
  found   <- genes[genes %in% rownames(obj)]
  missing <- genes[!genes %in% rownames(obj)]
  if (length(missing)) cat("Missing:", paste(missing, collapse = ", "), "\n")
  cat("Using  :", paste(found, collapse = ", "), "\n")
  found
}


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
    strip.background  = element_blank(),
    strip.text        = element_text(size = base_size + 1, face = "bold"),
    panel.spacing     = unit(3, "mm")
  )
}


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


compute_stats <- function(df, min_nonzero = 3) {
  df_testable <- df %>%
    group_by(gene, celltype, Diagnosis) %>%
    filter(sum(expr > 0) >= min_nonzero) %>%
    group_by(gene, celltype) %>%
    filter(n_distinct(Diagnosis) == 2) %>%
    ungroup()

  if (nrow(df_testable) == 0) return(tibble())

  stats <- df_testable %>%
    group_by(gene, celltype) %>%
    wilcox_test(expr ~ Diagnosis, exact = FALSE) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")

  y_caps <- df %>%
    group_by(gene, celltype) %>%
    summarise(y_cap = quantile(expr, 0.99), .groups = "drop")

  stats %>%
    add_xy_position(x = "Diagnosis") %>%
    left_join(y_caps, by = c("gene", "celltype")) %>%
    mutate(y.position = y_cap * 1.12)
}


.padj_labels <- function(stats, df) {
  if (nrow(stats) == 0) return(tibble())
  stats %>%
    left_join(
      df %>% group_by(celltype) %>%
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


plot_violin <- function(df, stats, title = "") {
  ann <- .padj_labels(stats, df)

  p <- ggplot(df, aes(x = Diagnosis, y = expr, fill = Diagnosis)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7, linewidth = 0.3) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA,
                 linewidth = 0.35) +
    scale_fill_manual(values = DX_COLORS, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    facet_wrap(~ celltype, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = "log-norm expression", title = title,
         subtitle = paste0("n = ", nrow(df), " cells")) +
    theme_nature() +
    theme(legend.position  = "top",
          legend.direction = "horizontal",
          legend.key.size  = unit(2.5, "mm"),
          legend.spacing.x = unit(1, "mm"),
          strip.text       = element_text(size = 7, face = "bold"))

  if (nrow(stats) > 0)
    p <- p + stat_pvalue_manual(stats, label = "p.adj.signif",
                                tip.length = 0.01, bracket.size = 0.3,
                                size = 2.2, hide.ns = FALSE)

  if (nrow(ann) > 0)
    p <- p + geom_text(data = ann, aes(x = x, y = y_ann, label = label),
                       inherit.aes = FALSE, parse = TRUE,
                       size = 1.8, color = "grey30")
  p
}


make_figure <- function(obj, genes, title, out_pdf, n_ct) {
  df    <- build_long_df(obj, genes)
  stats <- compute_stats(df)

  panels <- lapply(genes, function(g) {
    st <- if (nrow(stats) > 0) filter(stats, gene == g) else tibble()
    plot_violin(filter(df, gene == g), st,
                title = if (g == genes[1]) title else "")
  })

  p <- wrap_plots(panels, ncol = 1)

  w_mm <- 48 * n_ct + 24
  h_mm <- 90 * length(genes)

  ggsave(out_pdf, p,
         width  = w_mm / 25.4,
         height = h_mm / 25.4,
         device = "pdf", units = "in")
  cat("Saved:", basename(out_pdf), "\n")
}


# =============================================================================
# SECTION 1 — LOAD DATA & BUILD CACHES
# =============================================================================

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cache_all_fibro   <- file.path(DATA_DIR, "cache_all_fibro.rds")
cache_fibroblasts <- file.path(DATA_DIR, "cache_fibroblasts.rds")
cache_myo         <- file.path(DATA_DIR, "cache_myofibroblasts.rds")

if (!all(file.exists(cache_all_fibro, cache_fibroblasts, cache_myo))) {

  cat("\n[1/4] Loading full dataset...\n")
  seu <- load_geo_sparse(DATA_DIR, project = "GSE135893_ILD")

  cat("\n[2/4] Normalizing...\n")
  seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                       scale.factor = 10000)

  cat("\n[3/4] Creating subsets...\n")
  sub_all_fibro   <- subset_ctrl_ipf(seu, ALL_FIBRO_CT)
  sub_fibroblasts <- subset_ctrl_ipf(seu, "Fibroblasts")
  sub_myo         <- subset_ctrl_ipf(seu, "Myofibroblasts")

  cat("\n[4/4] Caching subsets...\n")
  saveRDS(sub_all_fibro,   cache_all_fibro)
  saveRDS(sub_fibroblasts, cache_fibroblasts)
  saveRDS(sub_myo,         cache_myo)
  print(gc())
  rm(seu); gc()

} else {
  cat("Loading cached subsets...\n")
  sub_all_fibro   <- readRDS(cache_all_fibro)
  sub_fibroblasts <- readRDS(cache_fibroblasts)
  sub_myo         <- readRDS(cache_myo)
}

genes_found <- check_genes(sub_all_fibro, GENES_OI)


# =============================================================================
# SECTION 2 — FIGURES
# =============================================================================

make_figure(
  sub_all_fibro, genes_found,
  title   = "All fibroblast clusters - IPF vs Control (Habermann 2020)",
  out_pdf = file.path(OUT_DIR, "v2_fig1_all_fibro_clusters.pdf"),
  n_ct    = length(ALL_FIBRO_CT)
)

make_figure(
  sub_fibroblasts, genes_found,
  title   = "Fibroblasts - IPF vs Control (Habermann 2020)",
  out_pdf = file.path(OUT_DIR, "v2_fig2_fibroblasts_only.pdf"),
  n_ct    = 1
)

make_figure(
  sub_myo, genes_found,
  title   = "Myofibroblasts - IPF vs Control (Habermann 2020)",
  out_pdf = file.path(OUT_DIR, "v2_fig3_myofibroblasts_only.pdf"),
  n_ct    = 1
)

cat("\nDone. All outputs in:", OUT_DIR, "\n")
