#!/usr/bin/env Rscript
# Update 3 — DEG + GSEA + metadata breakdown
# Loads best-param Seurat object from mac_subcluster_update3.R
# Gene set: highlight_genes.txt UNION ESR1/SP1
# Output: bm_analysis_out/figures_update_3/

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

BOLD_LEG <- theme(
  legend.title = element_text(face = "bold", size = 11),
  legend.text  = element_text(face = "bold", size = 10)
)

RDS_IN      <- "bm_analysis_out/figures_update_3/mac_update3_best.rds"
HIGHLIGHT_F <- "bm_analysis_out/figures/highlight_genes.txt"
PATIENT_XLS <- "bm_analysis_out/raw_geo/Table S1 patient infor_corrected_2025-08.xlsx"
OUT_DIR     <- "bm_analysis_out/figures_update_3"

if (!file.exists(RDS_IN)) stop("Run mac_subcluster_update3.R first.")

cat("Loading mac_update3_best.rds ...\n")
s <- readRDS(RDS_IN)
Idents(s) <- "mac_subcluster"
cat("  Clusters:", length(levels(s)), "| Cells:", ncol(s), "\n")
cat("  cancer_type_er:\n"); print(table(s$cancer_type_er))

# ── Gene list: highlight_genes.txt UNION ESR1/SP1 ─────────────────────────────
cat("\nBuilding combined gene list ...\n")
all_genes <- rownames(s)

raw_hl <- if (file.exists(HIGHLIGHT_F)) {
  raw <- readLines(HIGHLIGHT_F)
  trimws(unlist(strsplit(paste(raw, collapse=","), ",")))
} else character(0)
raw_hl <- raw_hl[nzchar(raw_hl)]

# resolve case: exact first, then TOUPPER
resolve_genes <- function(requested, available) {
  found <- intersect(requested, available)
  upper_miss <- toupper(requested[!requested %in% found])
  c(found, intersect(upper_miss, available))
}

hl_genes   <- resolve_genes(raw_hl, all_genes)
esr1_sp1   <- resolve_genes(c("ESR1","SP1"), all_genes)
comb_genes <- unique(c(hl_genes, esr1_sp1))

cat("  highlight_genes.txt →", paste(hl_genes,   collapse=", "), "\n")
cat("  ESR1/SP1           →", paste(esr1_sp1,   collapse=", "), "\n")
cat("  Combined           →", paste(comb_genes, collapse=", "), "\n")
if (length(comb_genes) == 0) stop("No genes found in Seurat object.")

gene_label    <- paste(comb_genes, collapse=" / ")
hl_label      <- paste(hl_genes, collapse=" / ")
esr1_sp1_label <- paste(esr1_sp1, collapse=" / ")

# ── Join treatment response ────────────────────────────────────────────────────
cat("\nJoining treatment response metadata ...\n")
pt <- read_excel(PATIENT_XLS, skip = 1)
colnames(pt) <- make.names(colnames(pt))

tx_map <- pt %>%
  select(cancer.id,
         neo_resp = Neoadjuvant.treatment.response,
         met_resp = Metastatic...treatment.response) %>%
  distinct() %>%
  mutate(tx_response = case_when(
    !tolower(neo_resp) %in% c("na","") ~ tolower(neo_resp),
    !tolower(met_resp) %in% c("na","") ~ tolower(met_resp),
    TRUE ~ "na"
  ))

tx_lookup <- setNames(tx_map$tx_response, tx_map$cancer.id)
s$tx_response <- unname(tx_lookup[s$cancer.id])
s$tx_response[is.na(s$tx_response)] <- "na"

s$gender_clean <- recode(s$gender, "f"="Female", "m"="Male", "na"="Unknown")

# ── 1. Feature plots ──────────────────────────────────────────────────────────
cat("\nGenerating feature plots ...\n")

# 3a: highlight genes only
if (length(hl_genes) > 0) {
  fp_hl <- FeaturePlot(s, features = hl_genes, reduction = "umap",
                       cols = c("lightgrey","#d62728"),
                       ncol = min(length(hl_genes), 3), pt.size = 0.2) &
           theme(legend.position = "right") & BOLD_LEG
  fw <- 6 * min(length(hl_genes), 3)
  fh <- 5 * ceiling(length(hl_genes) / 3)
  ggsave(file.path(OUT_DIR, "03a_feature_plots_highlight.pdf"), fp_hl, width=fw, height=fh)
  ggsave(file.path(OUT_DIR, "03a_feature_plots_highlight.png"), fp_hl, width=fw, height=fh, dpi=150)
  cat("  Saved: 03a_feature_plots_highlight\n")
}

# 3b: ESR1/SP1 only
if (length(esr1_sp1) > 0) {
  fp_es <- FeaturePlot(s, features = esr1_sp1, reduction = "umap",
                       cols = c("lightgrey","#2171b5"),
                       ncol = min(length(esr1_sp1), 3), pt.size = 0.2) &
           theme(legend.position = "right") & BOLD_LEG
  fw2 <- 6 * min(length(esr1_sp1), 3); fh2 <- 5 * ceiling(length(esr1_sp1) / 3)
  ggsave(file.path(OUT_DIR, "03b_feature_plots_esr1_sp1.pdf"), fp_es, width=fw2, height=fh2)
  ggsave(file.path(OUT_DIR, "03b_feature_plots_esr1_sp1.png"), fp_es, width=fw2, height=fh2, dpi=150)
  cat("  Saved: 03b_feature_plots_esr1_sp1\n")
}

# 3c: combined dot plot
dp_comb <- DotPlot(s, features = comb_genes, group.by = "mac_subcluster") +
  coord_flip() +
  labs(title = paste("Gene expression per macrophage sub-cluster\n", gene_label)) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0)) + BOLD_LEG
dp_h <- max(4, 0.35 * length(comb_genes) + 2)
ggsave(file.path(OUT_DIR, "03c_dotplot_combined.pdf"), dp_comb, width=10, height=dp_h)
ggsave(file.path(OUT_DIR, "03c_dotplot_combined.png"), dp_comb, width=10, height=dp_h, dpi=150)
cat("  Saved: 03c_dotplot_combined\n")

# Average expression summary
avg <- AverageExpression(s, features = comb_genes, group.by = "mac_subcluster")[["SCT"]]
scores <- colMeans(avg)
cat("  Gene-high cluster:", names(which.max(scores)), "\n")
cat("  Scores per cluster:\n"); print(round(sort(scores, decreasing=TRUE), 4))

# ── 2. DEG — FindAllMarkers ────────────────────────────────────────────────────
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
top5 <- markers %>% group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>% ungroup()
hm <- DoHeatmap(s, features = unique(top5$gene), group.by = "mac_subcluster",
                angle = 0, size = 3, draw.lines = TRUE) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu")), na.value = "white") +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(OUT_DIR, "04a_deg_heatmap.pdf"), hm, width=14, height=8)
ggsave(file.path(OUT_DIR, "04a_deg_heatmap.png"), hm, width=14, height=8, dpi=150)
cat("  Saved: 04a_deg_heatmap\n")

## Fig 4b: DEG dot plot — top 3 per cluster
top3 <- markers %>% group_by(cluster) %>%
  slice_max(avg_log2FC, n = 3) %>% ungroup()
dp_deg <- DotPlot(s, features = unique(top3$gene), group.by = "mac_subcluster") +
  coord_flip() +
  labs(title = "Top 3 DEGs per macrophage sub-cluster") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0), axis.text.y = element_text(size = 8)) + BOLD_LEG
ggsave(file.path(OUT_DIR, "04b_deg_dotplot.pdf"), dp_deg, width=12, height=8)
ggsave(file.path(OUT_DIR, "04b_deg_dotplot.png"), dp_deg, width=12, height=8, dpi=150)
cat("  Saved: 04b_deg_dotplot\n")

## Fig 4c: DEG dot plot — combined genes across clusters (for context)
dp_comb_cl1 <- DotPlot(s, features = comb_genes, group.by = "mac_subcluster",
                        idents = c("0","1","2","3","4","5")) +
  coord_flip() +
  labs(title = paste("Highlight + ESR1/SP1 in top 6 clusters\n", gene_label)) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0)) + BOLD_LEG
ggsave(file.path(OUT_DIR, "04c_dotplot_comb_top_clusters.pdf"), dp_comb_cl1, width=10, height=dp_h)
ggsave(file.path(OUT_DIR, "04c_dotplot_comb_top_clusters.png"), dp_comb_cl1, width=10, height=dp_h, dpi=150)
cat("  Saved: 04c_dotplot_comb_top_clusters\n")

# ── 3. GSEA ───────────────────────────────────────────────────────────────────
cat("\nLoading MSigDB gene sets (Hallmark + Oncogenic C6) ...\n")
msig_h  <- msigdbr(species = "Homo sapiens", collection = "H")
msig_c6 <- msigdbr(species = "Homo sapiens", collection = "C6")
pathways <- split(bind_rows(msig_h, msig_c6)$gene_symbol,
                  bind_rows(msig_h, msig_c6)$gs_name)
cat("  Gene sets:", length(pathways), "\n")

gsea_csv <- file.path(OUT_DIR, "gsea_cancer_results.csv")
if (file.exists(gsea_csv)) {
  cat("Loading cached GSEA results ...\n")
  gsea_all <- read.csv(gsea_csv)
  gsea_all$cluster <- as.character(gsea_all$cluster)
} else {
  cat("Running fgsea per cluster (HVGs only) ...\n")
  gsea_results <- lapply(levels(Idents(s)), function(cl) {
    cat("  Cluster", cl, "...\n")
    cl_markers <- FindMarkers(s, ident.1 = cl,
                              features        = VariableFeatures(s),
                              only.pos        = FALSE,
                              min.pct         = 0.1,
                              logfc.threshold = 0,
                              verbose         = FALSE)
    ranks <- sort(setNames(cl_markers$avg_log2FC, rownames(cl_markers)), decreasing=TRUE)
    res <- fgsea(pathways=pathways, stats=ranks,
                 minSize=15, maxSize=500, nPermSimple=1000, nproc=1)
    res$cluster <- cl; res
  })
  gsea_all <- bind_rows(gsea_results)
  write.csv(gsea_all %>% select(-leadingEdge), gsea_csv, row.names=FALSE)
}
cat("  GSEA rows:", nrow(gsea_all), "\n")

## Fig 5b: bubble plot
top_gsea <- gsea_all %>% filter(padj < 0.25) %>%
  group_by(cluster) %>% slice_max(abs(NES), n=3) %>% ungroup() %>%
  mutate(pathway_short = gsub("_"," ", gsub("^HALLMARK_","H: ", pathway)))

if (nrow(top_gsea) > 0) {
  bp <- ggplot(top_gsea, aes(x=cluster, y=pathway_short,
                              size=-log10(pmax(padj,1e-4)), color=NES)) +
    geom_point(alpha=0.85) +
    scale_color_gradientn(colors=rev(brewer.pal(11,"RdBu")),
                          limits=c(-3,3), oob=scales::squish, name="NES") +
    scale_size_continuous(name="-log10(padj)", range=c(2,8)) +
    labs(x="Macrophage sub-cluster", y=NULL,
         title="Top GSEA pathways per cluster (Hallmark + Oncogenic C6, padj<0.25)") +
    theme_classic(base_size=12) + theme(axis.text.y=element_text(size=9)) + BOLD_LEG
  bub_h <- max(5, 0.3 * length(unique(top_gsea$pathway_short)) + 2)
  ggsave(file.path(OUT_DIR,"05b_gsea_bubble.pdf"), bp, width=12, height=bub_h)
  ggsave(file.path(OUT_DIR,"05b_gsea_bubble.png"), bp, width=12, height=bub_h, dpi=150)
  cat("  Saved: 05b_gsea_bubble\n")
}

# ── 4. Metadata breakdown ─────────────────────────────────────────────────────
cat("\nGenerating metadata breakdown figures ...\n")

meta_bar <- function(df, fill_var, title, subtitle, fill_colors=NULL, filename) {
  dat <- df %>%
    group_by(mac_subcluster, .data[[fill_var]]) %>%
    summarise(n=n(), .groups="drop") %>%
    group_by(mac_subcluster) %>% mutate(frac=n/sum(n))
  p <- ggplot(dat, aes(x=mac_subcluster, y=frac, fill=.data[[fill_var]])) +
    geom_bar(stat="identity", width=0.8) +
    labs(x="Macrophage sub-cluster", y="Fraction",
         title=title, subtitle=subtitle, fill=fill_var) +
    theme_classic(base_size=13) + theme(axis.text.x=element_text(angle=0)) + BOLD_LEG
  if (!is.null(fill_colors)) p <- p + scale_fill_manual(values=fill_colors)
  else p <- p + scale_fill_brewer(palette="Set2")
  ggsave(file.path(OUT_DIR, paste0(filename,".pdf")), p, width=10, height=5)
  ggsave(file.path(OUT_DIR, paste0(filename,".png")), p, width=10, height=5, dpi=150)
  cat("  Saved:", filename, "\n")
}

meta_bar(s@meta.data, "gender_clean",
         title="Gender distribution per macrophage sub-cluster",
         subtitle="Female / Male / Unknown",
         fill_colors=c("Female"="#e377c2","Male"="#1f77b4","Unknown"="#aec7e8"),
         filename="06_gender_per_cluster")

meta_bar(s@meta.data, "tissue.origin",
         title="Tissue / bone site per macrophage sub-cluster",
         subtitle="Site of bone metastasis resection",
         filename="07_tissue_origin_per_cluster")

meta_bar(s@meta.data, "tx_response",
         title="Treatment response per macrophage sub-cluster",
         subtitle="Neoadjuvant or metastatic treatment response (Table S1)",
         filename="08_treatment_per_cluster")

cat("\n===== deg_gsea_update3.R complete =====\n")
cat("Genes (highlight):", hl_label, "\n")
cat("Genes (ESR1/SP1):", esr1_sp1_label, "\n")
cat("Figures in:", OUT_DIR, "\n")
