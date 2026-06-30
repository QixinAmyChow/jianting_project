#!/usr/bin/env Rscript
# Update 3 — Parameter sweep to separate BC / BC_ER cells
# Loads from mac_update2.rds (SCTransform + PCA already computed)
# Sweeps: UMAP n.neighbors × min.dist × clustering resolution
# Scores each combo by BC/BC_ER separation; saves best result + full panel
# Output: bm_analysis_out/figures_update_3/
# Previous outputs preserved: figures_R_v2/, figures_update_2/

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

RDS_IN  <- "bm_analysis_out/figures_update_2/mac_update2.rds"
OUT_DIR <- "bm_analysis_out/figures_update_3"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load pre-computed object ───────────────────────────────────────────────
cat("Loading mac_update2.rds (SCTransform + PCA already done) ...\n")
s <- readRDS(RDS_IN)
cat("  Cells:", ncol(s), "\n")
cat("  cancer_type_er distribution:\n")
print(table(s$cancer_type_er))

# ── 2. Parameter grid ─────────────────────────────────────────────────────────
# UMAP params: n.neighbors controls global structure; min.dist controls
# tightness of local clusters. Lower min.dist = tighter clusters.
nn_vals   <- c(10, 20, 30, 50)
md_vals   <- c(0.05, 0.1, 0.3)
res_vals  <- c(0.3, 0.5, 0.8, 1.0, 1.2, 1.5)
n_pcs     <- 1:20   # consistent with update2

grid <- expand.grid(n_neighbors = nn_vals,
                    min_dist    = md_vals,
                    resolution  = res_vals,
                    stringsAsFactors = FALSE)
cat("\nGrid size:", nrow(grid), "combinations\n")
cat("  n.neighbors:", paste(nn_vals, collapse=", "), "\n")
cat("  min.dist:   ", paste(md_vals, collapse=", "), "\n")
cat("  resolution: ", paste(res_vals, collapse=", "), "\n\n")

# ── 3. Scoring function ───────────────────────────────────────────────────────
# For each clustering, find the cluster with highest BC/BC_ER enrichment.
# Score = F1-like: harmonic mean of
#   precision = fraction of that cluster that is BC/BC_ER
#   recall    = fraction of all BC/BC_ER cells captured in that cluster
# Computed separately for BC_ER+ and for (BC_ER+ + BC_ER- + BC) combined.
score_separation <- function(cancer_type_er, clusters) {
  df <- data.frame(ct = cancer_type_er, cl = clusters, stringsAsFactors = FALSE)

  # target groups
  targets <- list(
    BC_ER_pos = "BC_ER+",
    BC_family = c("BC_ER+", "BC_ER-", "BC_ER?", "BC")
  )

  scores <- sapply(targets, function(tgt) {
    n_tgt   <- sum(df$ct %in% tgt)
    if (n_tgt == 0) return(0)
    # per cluster: precision and recall
    cl_stats <- df %>%
      group_by(cl) %>%
      summarise(n_cl  = n(),
                n_hit = sum(ct %in% tgt), .groups = "drop") %>%
      mutate(prec = n_hit / n_cl,
             rec  = n_hit / n_tgt)
    # F1 per cluster, take max
    f1_per <- with(cl_stats, ifelse(prec + rec > 0,
                                    2 * prec * rec / (prec + rec), 0))
    max(f1_per)
  })
  scores
}

# ── 4. Run grid search ────────────────────────────────────────────────────────
results <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  nn  <- grid$n_neighbors[i]
  md  <- grid$min_dist[i]
  res <- grid$resolution[i]

  tag <- sprintf("nn%02d_md%.2f_res%.1f", nn, md, res)

  # Run UMAP with these params (uses PCA already in s)
  s_tmp <- RunUMAP(s, dims = n_pcs, n.neighbors = nn, min.dist = md,
                   reduction.name = "umap_tmp",
                   reduction.key  = "UMAPtmp_",
                   verbose = FALSE)

  # Cluster
  s_tmp <- FindNeighbors(s_tmp, dims = n_pcs, verbose = FALSE)
  s_tmp <- FindClusters(s_tmp, resolution = res, verbose = FALSE,
                        cluster.name = "cl_tmp")

  ncl <- length(unique(s_tmp$cl_tmp))
  sc  <- score_separation(s_tmp$cancer_type_er, s_tmp$cl_tmp)

  results[[i]] <- list(
    tag        = tag,
    n_neighbors = nn,
    min_dist   = md,
    resolution = res,
    n_clusters = ncl,
    f1_BC_ER_pos  = sc["BC_ER_pos"],
    f1_BC_family  = sc["BC_family"],
    umap_coords   = Embeddings(s_tmp, "umap_tmp"),
    clusters      = s_tmp$cl_tmp
  )

  if (i %% 10 == 0 || i == nrow(grid))
    cat(sprintf("  [%d/%d] %s  ncl=%d  F1_ER+=%.3f  F1_BC=%.3f\n",
                i, nrow(grid), tag, ncl,
                sc["BC_ER_pos"], sc["BC_family"]))
}

# ── 5. Score table ────────────────────────────────────────────────────────────
score_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(tag        = r$tag,
             n_neighbors = r$n_neighbors,
             min_dist   = r$min_dist,
             resolution = r$resolution,
             n_clusters = r$n_clusters,
             f1_BC_ER_pos  = r$f1_BC_ER_pos,
             f1_BC_family  = r$f1_BC_family,
             stringsAsFactors = FALSE)
}))

score_df <- score_df %>% arrange(desc(f1_BC_ER_pos))
cat("\nTop 15 by BC_ER+ F1-score:\n")
print(head(score_df, 15), row.names = FALSE)

write.csv(score_df, file.path(OUT_DIR, "00_grid_scores.csv"), row.names = FALSE)
cat("  Saved: 00_grid_scores.csv\n")

# ── 6. Panel of top-12 UMAPs ──────────────────────────────────────────────────
CANCER_COLORS <- c(
  "BC_ER+"  = "#d62728", "BC_ER-"  = "#ff9896", "BC_ER?"  = "#fa9fb5",
  "BC"      = "#fb6a4a",
  "KC"      = "#8c564b", "LC"      = "#17becf",
  "CC"      = "#2ca02c", "BDC"     = "#ff7f0e",
  "EC"      = "#9467bd", "TC"      = "#7f7f7f",
  "PC"      = "#bcbd22", "ctrl"    = "#aec7e8"
)

top12_idx <- match(head(score_df$tag, 12), sapply(results, `[[`, "tag"))

panel_plots <- lapply(top12_idx, function(i) {
  r   <- results[[i]]
  coords <- as.data.frame(r$umap_coords)
  colnames(coords) <- c("UMAP1","UMAP2")
  coords$cancer_type_er <- s$cancer_type_er
  coords$cluster        <- r$clusters

  ggplot(coords, aes(UMAP1, UMAP2, color = cancer_type_er)) +
    geom_point(size = 0.2, alpha = 0.6) +
    scale_color_manual(values = CANCER_COLORS) +
    ggtitle(sprintf("%s\nncl=%d  F1_ER+=%.3f",
                    r$tag, r$n_clusters, r$f1_BC_ER_pos)) +
    theme_classic(base_size = 7) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))
})

panel <- wrap_plots(panel_plots, ncol = 4)
ggsave(file.path(OUT_DIR, "00_top12_panel.pdf"), panel, width = 20, height = 12)
ggsave(file.path(OUT_DIR, "00_top12_panel.png"), panel, width = 20, height = 12, dpi = 120)
cat("  Saved: 00_top12_panel\n")

# ── 7. Best config — full figures ─────────────────────────────────────────────
best <- results[[top12_idx[1]]]
cat("\nBest config:", best$tag,
    "  n_neighbors=", best$n_neighbors,
    "  min_dist=", best$min_dist,
    "  resolution=", best$resolution,
    "  n_clusters=", best$n_clusters,
    "  F1_ER+=", round(best$f1_BC_ER_pos, 4),
    "\n")

# Re-run on s to get a clean Seurat object with the best params
cat("Re-running best config on full object ...\n")
s_best <- RunUMAP(s, dims = n_pcs,
                  n.neighbors = best$n_neighbors,
                  min.dist    = best$min_dist,
                  verbose     = FALSE)
s_best <- FindNeighbors(s_best, dims = n_pcs, verbose = FALSE)
s_best <- FindClusters(s_best, resolution = best$resolution,
                       verbose = FALSE, cluster.name = "mac_subcluster")
Idents(s_best) <- "mac_subcluster"

cat("  Clusters:", length(unique(s_best$mac_subcluster)), "\n")
print(table(s_best$mac_subcluster))
cat("  Cancer type distribution per cluster:\n")
print(table(s_best$mac_subcluster, s_best$cancer_type_er))

## Fig 1: UMAP overview
p1 <- DimPlot(s_best, reduction = "umap", group.by = "mac_subcluster",
              label = TRUE, repel = TRUE, pt.size = 0.3) +
      ggtitle(sprintf("Macrophage sub-clusters\n(nn=%d, md=%.2f, res=%.1f)",
                      best$n_neighbors, best$min_dist, best$resolution)) +
      theme(legend.position = "right")

p2 <- DimPlot(s_best, reduction = "umap", group.by = "cancer_type_er",
              cols = CANCER_COLORS, pt.size = 0.3) +
      ggtitle("Cancer type (BC split by ER status)") +
      theme(legend.position = "right")

fig1 <- p1 | p2
ggsave(file.path(OUT_DIR, "01_umap_overview.pdf"), fig1, width = 14, height = 6)
ggsave(file.path(OUT_DIR, "01_umap_overview.png"), fig1, width = 14, height = 6, dpi = 150)
cat("  Saved: 01_umap_overview\n")

## Fig 2: Stacked bar — cancer type per cluster
meta <- s_best@meta.data %>%
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
       subtitle = sprintf("nn=%d  min.dist=%.2f  res=%.1f  (Update 3 best param)",
                          best$n_neighbors, best$min_dist, best$resolution)) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0))

ggsave(file.path(OUT_DIR, "02_cancer_type_bar_bc_er.pdf"), fig2, width = 10, height = 5)
ggsave(file.path(OUT_DIR, "02_cancer_type_bar_bc_er.png"), fig2, width = 10, height = 5, dpi = 150)
cat("  Saved: 02_cancer_type_bar_bc_er\n")

## Fig 3: Highlight genes (reuse list from update2 if available)
HIGHLIGHT_F <- "bm_analysis_out/figures/highlight_genes.txt"
if (file.exists(HIGHLIGHT_F)) {
  raw_genes <- readLines(HIGHLIGHT_F)
  raw_genes <- trimws(unlist(strsplit(paste(raw_genes, collapse=","), ",")))
  raw_genes <- raw_genes[nzchar(raw_genes)]
  all_genes  <- rownames(s_best)
  found <- intersect(raw_genes, all_genes)
  found <- c(found, intersect(toupper(raw_genes[!raw_genes %in% found]), all_genes))
  if (length(found) > 0) {
    fp <- FeaturePlot(s_best, features = found, reduction = "umap",
                      cols = c("lightgrey","#d62728"),
                      ncol = min(length(found), 3), pt.size = 0.2) &
          theme(legend.position = "right")
    fw <- 6 * min(length(found), 3)
    fh <- 5 * ceiling(length(found) / 3)
    ggsave(file.path(OUT_DIR, "03a_feature_plots_highlight.pdf"), fp, width=fw, height=fh)
    ggsave(file.path(OUT_DIR, "03a_feature_plots_highlight.png"), fp, width=fw, height=fh, dpi=150)

    dp <- DotPlot(s_best, features = found, group.by = "mac_subcluster") +
          coord_flip() +
          labs(title = paste(paste(found, collapse="/"), "across macrophage sub-clusters")) +
          theme_classic(base_size = 12) +
          theme(axis.text.x = element_text(angle = 0))
    dp_h <- max(4, 0.35 * length(found) + 2)
    ggsave(file.path(OUT_DIR, "03b_dotplot_highlight.pdf"), dp, width=10, height=dp_h)
    ggsave(file.path(OUT_DIR, "03b_dotplot_highlight.png"), dp, width=10, height=dp_h, dpi=150)
    cat("  Saved: 03a/03b highlight genes\n")
  }
}

## Fig 4: BC_ER+ highlighted alone on UMAP (to see cluster separation clearly)
coords_best <- as.data.frame(Embeddings(s_best, "umap"))
colnames(coords_best) <- c("UMAP1","UMAP2")
coords_best$cancer_type_er <- s_best$cancer_type_er
coords_best$highlight <- ifelse(coords_best$cancer_type_er %in% c("BC_ER+","BC_ER-","BC_ER?","BC"),
                                 coords_best$cancer_type_er, "other")
coords_best$highlight <- factor(coords_best$highlight,
                                 levels = c("BC_ER+","BC_ER-","BC_ER?","BC","other"))
hl_colors <- c("BC_ER+" = "#d62728", "BC_ER-" = "#ff9896",
                "BC_ER?" = "#fa9fb5", "BC" = "#fb6a4a", "other" = "grey85")

fig4 <- ggplot(coords_best %>% arrange(highlight == "other"),
               aes(UMAP1, UMAP2, color = highlight)) +
  geom_point(size = 0.3, alpha = 0.7) +
  scale_color_manual(values = hl_colors, name = "BC status") +
  ggtitle("BC / BC_ER cells highlighted") +
  theme_classic(base_size = 13) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(file.path(OUT_DIR, "04_bc_highlight_umap.pdf"), fig4, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "04_bc_highlight_umap.png"), fig4, width = 8, height = 6, dpi = 150)
cat("  Saved: 04_bc_highlight_umap\n")

# ── 8. Save best Seurat object ────────────────────────────────────────────────
rds_out <- file.path(OUT_DIR, "mac_update3_best.rds")
saveRDS(s_best, rds_out)
cat("  Saved RDS:", rds_out, "\n")

cat("\n===== Update 3 complete =====\n")
cat("Best params: n.neighbors=", best$n_neighbors,
    "  min.dist=", best$min_dist,
    "  resolution=", best$resolution, "\n")
cat("n_clusters:", best$n_clusters,
    "  F1(BC_ER+)=", round(best$f1_BC_ER_pos, 4),
    "  F1(BC_family)=", round(best$f1_BC_family, 4), "\n")
cat("Figures in:", OUT_DIR, "\n")
cat("Previous outputs preserved: figures_R_v2/, figures_update_2/\n")
