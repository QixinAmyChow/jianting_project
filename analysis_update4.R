#!/usr/bin/env Rscript
# ============================================================
# BM Macrophage Update 4 — New comparative analyses
#
# 1. BC macrophage vs Other cancer macrophages (DEG + GSEA)
# 2. BC ER+ vs BC ER- macrophage comparison (DEG + GSEA)
# 3. 4-group metadata: BC_mac / Other_mac / BC_Tcell / Other_Tcell
# 4. Export BC macrophages for SpaTrack (run spatrack_bc_mac.py next)
#
# Input:  bm_analysis_out/figures_update_3/mac_update3_best.rds  (required)
#         T_cell.rds (optional – section 3 T cell portion added when present)
# Output: bm_analysis_out/figures_update_4/
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(fgsea)
  library(msigdbr)
  library(ggrepel)
  library(Matrix)
})
options(future.globals.maxSize = 4 * 1024^3)

# ── Paths ─────────────────────────────────────────────────────────────────────
CHECKPOINT <- "bm_analysis_out/figures_update_3/mac_update3_best.rds"
TCELL_RDS  <- "bm_analysis_out/raw_geo/integrated_Seurat_objects/47.integrated_object_subset_by_major_celltypes/T_cell.rds"
OUT_DIR    <- "bm_analysis_out/figures_update_4"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(CHECKPOINT))
  stop("Checkpoint not found: ", CHECKPOINT)

cat("Loading mac_update3_best.rds ...\n")
mac <- readRDS(CHECKPOINT)
cat("  Cells:", ncol(mac),
    "| Clusters:", nlevels(factor(mac$mac_subcluster)), "\n")


shorten_label <- function(x, n = 55) {
  x <- gsub("^HALLMARK_", "H: ",  x)
  x <- gsub("^GOBP_",     "GO: ", x)
  x <- gsub("^REACTOME_", "R: ",  x)
  x <- gsub("^GSE[0-9]+_", "",    x)
  x <- gsub("_", " ", x)
  x <- gsub("  +", " ", trimws(x))
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

save_fig <- function(p, name, w = 10, h = 6) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  cat("  Saved:", name, "\n")
}

# ── Load MSigDB gene sets once ────────────────────────────────────────────────
cat("Loading MSigDB gene sets ...\n")
get_gs <- function(cat, subcat = NULL) {
  args <- list(species = "Homo sapiens", category = cat)
  if (!is.null(subcat)) args$subcategory <- subcat
  gs <- do.call(msigdbr, args)
  split(gs$gene_symbol, gs$gs_name)
}
all_gs <- c(
  get_gs("H"),
  get_gs("C6"),
  get_gs("C5", "GO:BP"),
  get_gs("C2", "CP:REACTOME"),
  get_gs("C7", "IMMUNESIGDB")
)
cat("  Gene sets:", length(all_gs), "\n")

# ── Helper: fgsea from FindMarkers result ─────────────────────────────────────
run_fgsea <- function(markers, gs_list, nperm = 1000) {
  markers <- markers %>% filter(!is.na(avg_log2FC))
  rank_vec <- setNames(
    sign(markers$avg_log2FC) * (-log10(markers$p_val + 1e-300)),
    rownames(markers)
  )
  fgsea(pathways = gs_list, stats = sort(rank_vec, decreasing = TRUE),
        nPermSimple = nperm, eps = 0)
}

# ── Helper: volcano plot ──────────────────────────────────────────────────────
volcano_plot <- function(markers, title, top_n = 25, fc_cut = 0.5, padj_cut = 0.05) {
  df <- markers %>%
    mutate(
      direction  = case_when(
        avg_log2FC > fc_cut  & p_val_adj < padj_cut ~ "Up",
        avg_log2FC < -fc_cut & p_val_adj < padj_cut ~ "Down",
        TRUE ~ "NS"
      ),
      neg_log10p = -log10(p_val_adj + 1e-300),
      gene       = rownames(.)
    )
  top_genes <- df %>%
    filter(direction != "NS") %>%
    arrange(desc(abs(avg_log2FC))) %>%
    head(top_n)
  ggplot(df, aes(avg_log2FC, neg_log10p, color = direction)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 3, max.overlaps = 20) +
    scale_color_manual(values = c("Up" = "#d62728", "Down" = "#1f77b4",
                                  "NS" = "grey75")) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", alpha = 0.4) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", alpha = 0.4) +
    labs(title = title, x = "avg log2 FC (ident.1 / ident.2)",
         y = "-log10(padj)", color = "Direction") +
    theme_classic(base_size = 13)
}

# ── Helper: GSEA bubble plot ──────────────────────────────────────────────────
gsea_bubble <- function(gsea_df, title, top_n = 25, padj_cut = 0.05) {
  sig <- gsea_df %>%
    filter(padj < padj_cut) %>%
    arrange(padj) %>%
    head(top_n) %>%
    mutate(label     = shorten_label(pathway),
           direction = ifelse(NES > 0, "Up", "Down"))
  if (nrow(sig) == 0) {
    message("  No GSEA terms pass padj < ", padj_cut, " for: ", title)
    return(NULL)
  }
  sig$label <- factor(sig$label, levels = sig$label[order(sig$NES)])
  ggplot(sig, aes(NES, label, size = abs(NES), color = direction)) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = c("Up" = "#d62728", "Down" = "#1f77b4")) +
    scale_size_continuous(range = c(3, 10)) +
    labs(title = title, x = "NES", y = NULL, size = "|NES|") +
    theme_classic(base_size = 12) +
    theme(axis.text.y = element_text(size = 9))
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: BC macrophage vs Other macrophage
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Section 2: BC mac vs Other mac (DEG + GSEA) ===\n")

mac$mac_vs_other <- dplyr::case_when(
  mac$cancer == "BC"   ~ "BC_mac",
  mac$cancer == "ctrl" ~ NA_character_,
  TRUE                 ~ "Other_mac"
)
Idents(mac) <- "mac_vs_other"

n_bc  <- sum(mac$mac_vs_other == "BC_mac",    na.rm = TRUE)
n_oth <- sum(mac$mac_vs_other == "Other_mac", na.rm = TRUE)
cat("  BC_mac:", n_bc, "| Other_mac:", n_oth, "\n")

DEG_BC_F  <- file.path(OUT_DIR, "deg_bc_vs_other.csv")
GSEA_BC_F <- file.path(OUT_DIR, "gsea_bc_vs_other.csv")

if (file.exists(DEG_BC_F)) {
  cat("  Loading cached DEG ...\n")
  deg_bc <- read.csv(DEG_BC_F, row.names = 1)
} else {
  cat("  FindMarkers (BC_mac vs Other_mac) ...\n")
  deg_bc <- FindMarkers(mac, ident.1 = "BC_mac", ident.2 = "Other_mac",
                        min.pct = 0.1, logfc.threshold = 0.1, test.use = "wilcox")
  write.csv(deg_bc, DEG_BC_F)
}
cat("  DEG rows:", nrow(deg_bc), "\n")

if (file.exists(GSEA_BC_F)) {
  cat("  Loading cached GSEA ...\n")
  gsea_bc <- read.csv(GSEA_BC_F)
} else {
  cat("  Running GSEA ...\n")
  gsea_bc_full <- run_fgsea(deg_bc, all_gs)
  gsea_bc <- gsea_bc_full %>% select(-leadingEdge)
  write.csv(gsea_bc, GSEA_BC_F, row.names = FALSE)
}
cat("  GSEA rows:", nrow(gsea_bc), "| sig (padj<0.05):",
    sum(gsea_bc$padj < 0.05, na.rm = TRUE), "\n")

# Volcano
p_vol_bc <- volcano_plot(deg_bc, "BC macrophages vs Other cancer macrophages")
save_fig(p_vol_bc, "10_bc_vs_other_volcano")

# GSEA bubble
p_gsea_bc <- gsea_bubble(gsea_bc, "BC vs Other Mac — Top Pathways (GSEA)")
if (!is.null(p_gsea_bc)) save_fig(p_gsea_bc, "11_bc_vs_other_gsea", w = 11, h = 8)

# Top DEG dotplot (top 30 by |log2FC|, padj < 0.05)
top_bc_genes <- deg_bc %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(30) %>%
  rownames()

if (length(top_bc_genes) >= 5) {
  mac_sub <- subset(mac, !is.na(mac_vs_other))
  Idents(mac_sub) <- "mac_vs_other"
  p_dot_bc <- DotPlot(mac_sub, features = top_bc_genes, group.by = "mac_vs_other") +
    coord_flip() +
    labs(title = "Top DEGs: BC mac vs Other mac") +
    scale_color_gradient2(low = "#1f77b4", mid = "white", high = "#d62728",
                          midpoint = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p_dot_bc, "12_bc_vs_other_dotplot", w = 7, h = 10)
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: BC ER+ vs BC ER- macrophage comparison
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Section 3: BC ER+ vs BC ER- macrophages (DEG + GSEA) ===\n")

mac_bc <- subset(mac, cancer_type_er %in% c("BC_ER+", "BC_ER-"))
Idents(mac_bc) <- "cancer_type_er"
n_erp <- sum(mac_bc$cancer_type_er == "BC_ER+")
n_ern <- sum(mac_bc$cancer_type_er == "BC_ER-")
cat("  BC_ER+:", n_erp, "| BC_ER-:", n_ern, "\n")

DEG_ER_F  <- file.path(OUT_DIR, "deg_erpos_vs_erneg.csv")
GSEA_ER_F <- file.path(OUT_DIR, "gsea_erpos_vs_erneg.csv")

if (file.exists(DEG_ER_F)) {
  cat("  Loading cached DEG ...\n")
  deg_er <- read.csv(DEG_ER_F, row.names = 1)
} else {
  cat("  FindMarkers (BC_ER+ vs BC_ER-) ...\n")
  deg_er <- FindMarkers(mac_bc, ident.1 = "BC_ER+", ident.2 = "BC_ER-",
                        min.pct = 0.1, logfc.threshold = 0.1, test.use = "wilcox")
  write.csv(deg_er, DEG_ER_F)
}
cat("  DEG rows:", nrow(deg_er), "\n")

if (file.exists(GSEA_ER_F)) {
  cat("  Loading cached GSEA ...\n")
  gsea_er <- read.csv(GSEA_ER_F)
} else {
  cat("  Running GSEA ...\n")
  gsea_er_full <- run_fgsea(deg_er, all_gs)
  gsea_er <- gsea_er_full %>% select(-leadingEdge)
  write.csv(gsea_er, GSEA_ER_F, row.names = FALSE)
}
cat("  GSEA rows:", nrow(gsea_er), "| sig (padj<0.05):",
    sum(gsea_er$padj < 0.05, na.rm = TRUE), "\n")

# Volcano
p_vol_er <- volcano_plot(deg_er, "BC ER+ vs BC ER- macrophages")
save_fig(p_vol_er, "13_erpos_vs_erneg_volcano")

# GSEA bubble
p_gsea_er <- gsea_bubble(gsea_er, "BC ER+ vs ER- Mac — Top Pathways (GSEA)")
if (!is.null(p_gsea_er)) save_fig(p_gsea_er, "14_erpos_vs_erneg_gsea", w = 11, h = 8)

# Top DEG dotplot
top_er_genes <- deg_er %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(30) %>%
  rownames()

if (length(top_er_genes) >= 5) {
  p_dot_er <- DotPlot(mac_bc, features = top_er_genes, group.by = "cancer_type_er") +
    coord_flip() +
    labs(title = "Top DEGs: BC ER+ vs ER- macrophages") +
    scale_color_gradient2(low = "#1f77b4", mid = "white", high = "#d62728",
                          midpoint = 0) +
    scale_fill_manual(values = c("BC_ER+" = "#d62728", "BC_ER-" = "#1f77b4")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p_dot_er, "15_erpos_vs_erneg_dotplot", w = 7, h = 10)
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: 4-group metadata (mac portion now; T cell added when available)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Section 4: 4-group metadata ===\n")

mac$group4 <- dplyr::case_when(
  mac$cancer == "BC"   ~ "BC_mac",
  mac$cancer == "ctrl" ~ NA_character_,
  TRUE                 ~ "Other_mac"
)

cat("  Mac group4 distribution:\n")
print(table(mac$group4, useNA = "ifany"))

# Bar chart: mac group across clusters
grp_dat <- mac@meta.data %>%
  filter(!is.na(group4)) %>%
  group_by(mac_subcluster, group4) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(mac_subcluster) %>%
  mutate(frac = n / sum(n))

p_grp4 <- ggplot(grp_dat, aes(x = mac_subcluster, y = frac, fill = group4)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("BC_mac" = "#d62728", "Other_mac" = "#4878cf")) +
  labs(title = "BC_mac vs Other_mac per macrophage sub-cluster",
       x = "mac_subcluster", y = "Fraction", fill = "Group") +
  theme_classic(base_size = 13) +
  theme(legend.title = element_text(face = "bold"))
save_fig(p_grp4, "16_group4_mac_per_cluster", w = 10, h = 5)

# T cell portion (if T_cell.rds is already present)
if (file.exists(TCELL_RDS)) {
  cat("  Loading T_cell.rds ...\n")
  tcell <- readRDS(TCELL_RDS)
  cat("  T cells:", ncol(tcell), "\n")

  if ("cancer" %in% colnames(tcell@meta.data)) {
    tcell$group4 <- dplyr::case_when(
      tcell$cancer == "BC"   ~ "BC_Tcell",
      tcell$cancer == "ctrl" ~ NA_character_,
      TRUE                   ~ "Other_Tcell"
    )
    cat("  T cell group4:\n"); print(table(tcell$group4, useNA = "ifany"))

    summary_df <- bind_rows(
      as.data.frame(table(mac$group4))   %>% rename(group = Var1) %>% mutate(cell_class = "Macrophage"),
      as.data.frame(table(tcell$group4)) %>% rename(group = Var1) %>% mutate(cell_class = "T cell")
    ) %>% filter(!is.na(group))

    p_grp4_all <- ggplot(summary_df, aes(x = group, y = Freq, fill = group)) +
      geom_bar(stat = "identity", width = 0.7) +
      facet_wrap(~cell_class, scales = "free") +
      scale_fill_manual(values = c(
        "BC_mac"    = "#d62728", "Other_mac"   = "#4878cf",
        "BC_Tcell"  = "#ff7f0e", "Other_Tcell" = "#2ca02c"
      )) +
      labs(title = "4-Group Cell Counts: BC / Other × Mac / T cell",
           x = NULL, y = "Cell count", fill = "Group") +
      theme_classic(base_size = 13) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    save_fig(p_grp4_all, "16b_group4_all_counts", w = 9, h = 5)

    saveRDS(list(mac_group4   = table(mac$group4),
                 tcell_group4 = table(tcell$group4)),
            file.path(OUT_DIR, "group4_summary.rds"))
    cat("  Saved: group4_summary.rds\n")
    rm(tcell); gc()
  } else {
    cat("  WARNING: T_cell.rds missing 'cancer' column — T cell groups skipped\n")
  }
} else {
  cat("  T_cell.rds not yet available — T cell groups (BC_Tcell/Other_Tcell) will be\n")
  cat("  added on next run after download_tcell.sh completes.\n")
}

# Save mac with group4
saveRDS(mac, file.path(OUT_DIR, "mac_update4_group4.rds"))
cat("  Saved: mac_update4_group4.rds\n")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5: Export BC macrophage data for SpaTrack (Python)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Section 5: Export BC macrophages for SpaTrack ===\n")

{
  suppressPackageStartupMessages(library(Matrix))

  SPATRACK_DIR <- "bm_analysis_out/spatrack"
  dir.create(SPATRACK_DIR, showWarnings = FALSE, recursive = TRUE)

  mac_bc_exp <- subset(mac, cancer == "BC")
  n_bc_cells <- ncol(mac_bc_exp)
  cat("  BC macrophages:", n_bc_cells, "\n")
  cat("  Cluster distribution:\n")
  print(sort(table(as.character(mac_bc_exp$mac_subcluster)), decreasing = TRUE))

  # Identify start cluster (most Mono-like archetype)
  start_clust <- mac_bc_exp@meta.data %>%
    group_by(mac_subcluster) %>%
    summarise(mono_frac = mean(archetype == "Mono", na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    filter(n >= 5) %>%
    arrange(desc(mono_frac)) %>%
    slice(1) %>%
    pull(mac_subcluster) %>%
    as.character()
  cat("  Start cluster (most Mono-like):", start_clust, "\n")
  writeLines(start_clust, file.path(SPATRACK_DIR, "start_cluster.txt"))

  # UMAP embeddings
  umap_df <- as.data.frame(Embeddings(mac_bc_exp, "umap"))
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  write.csv(umap_df, file.path(SPATRACK_DIR, "umap.csv"))

  # PCA embeddings (first 30 PCs)
  pca_df <- as.data.frame(Embeddings(mac_bc_exp, "pca")[, 1:30])
  write.csv(pca_df, file.path(SPATRACK_DIR, "pca.csv"))

  # Cell metadata
  meta_cols <- intersect(
    c("mac_subcluster", "cancer_type_er", "archetype", "patient.id", "mac_vs_other"),
    colnames(mac_bc_exp@meta.data)
  )
  write.csv(mac_bc_exp@meta.data[, meta_cols, drop = FALSE],
            file.path(SPATRACK_DIR, "metadata.csv"))

  # Normalized SCT expression for top variable genes (sparse, MTX format)
  vg <- head(VariableFeatures(mac_bc_exp), 2000)
  expr_mat <- tryCatch(
    GetAssayData(mac_bc_exp, assay = "SCT", layer = "data")[vg, ],
    error = function(e)
      GetAssayData(mac_bc_exp, assay = "SCT", slot = "data")[vg, ]
  )
  # writeMM expects cells × genes (transpose)
  Matrix::writeMM(t(expr_mat), file.path(SPATRACK_DIR, "counts.mtx"))
  writeLines(rownames(expr_mat), file.path(SPATRACK_DIR, "genes.txt"))
  writeLines(colnames(expr_mat), file.path(SPATRACK_DIR, "barcodes.txt"))

  cat("  Exported to:", SPATRACK_DIR, "\n")
  cat("  Files: umap.csv, pca.csv, metadata.csv, counts.mtx, genes.txt,",
      "barcodes.txt, start_cluster.txt\n")
  cat("  Next: conda run -n scanpy_stable python spatrack_bc_mac.py\n")
}

cat("\n===== analysis_update4.R complete =====\n")
cat("All figures in:", OUT_DIR, "\n")
if (!file.exists(TCELL_RDS)) {
  cat("NOTE: T_cell.rds not yet present — re-run after download_tcell.sh\n")
  cat("      completes to add BC_Tcell/Other_Tcell groups (Section 4).\n")
}
