#!/usr/bin/env Rscript
# Regenerate clustering figures with bold legends
# Loads mac_update3_best.rds — no re-clustering needed
# Replaces: 01, 02, 03a, 03b, 03c, 04_bc_highlight in figures_update_3/

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

RDS_IN      <- "bm_analysis_out/figures_update_3/mac_update3_best.rds"
HIGHLIGHT_F <- "bm_analysis_out/figures/highlight_genes.txt"
OUT_DIR     <- "bm_analysis_out/figures_update_3"

BEST_NN  <- 10; BEST_MD <- 0.05; BEST_RES <- 0.3

CANCER_COLORS <- c(
  "BC_ER+"  = "#d62728", "BC_ER-"  = "#ff9896", "BC_ER?"  = "#fa9fb5",
  "BC"      = "#fb6a4a",
  "KC"      = "#8c564b", "LC"      = "#17becf",
  "CC"      = "#2ca02c", "BDC"     = "#ff7f0e",
  "EC"      = "#9467bd", "TC"      = "#7f7f7f",
  "PC"      = "#bcbd22", "ctrl"    = "#aec7e8"
)

BOLD_LEG <- theme(
  legend.title = element_text(face = "bold", size = 11),
  legend.text  = element_text(face = "bold", size = 10)
)

cat("Loading RDS ...\n")
s <- readRDS(RDS_IN)
Idents(s) <- "mac_subcluster"

# resolve genes
all_genes <- rownames(s)
raw_hl <- trimws(unlist(strsplit(paste(readLines(HIGHLIGHT_F), collapse=","), ",")))
raw_hl <- raw_hl[nzchar(raw_hl)]
resolve <- function(g, av) c(intersect(g, av), intersect(toupper(g[!g %in% av]), av))
hl_genes   <- resolve(raw_hl, all_genes)
esr1_sp1   <- resolve(c("ESR1","SP1"), all_genes)
comb_genes <- unique(c(hl_genes, esr1_sp1))

cat("  Genes:", paste(comb_genes, collapse=", "), "\n")

# ── Fig 01: UMAP overview ─────────────────────────────────────────────────────
p1 <- DimPlot(s, reduction = "umap", group.by = "mac_subcluster",
              label = TRUE, repel = TRUE, pt.size = 0.3) +
      ggtitle(sprintf("Macrophage sub-clusters\n(nn=%d, md=%.2f, res=%.1f)",
                      BEST_NN, BEST_MD, BEST_RES)) +
      theme(legend.position = "right") + BOLD_LEG

p2 <- DimPlot(s, reduction = "umap", group.by = "cancer_type_er",
              cols = CANCER_COLORS, pt.size = 0.3) +
      ggtitle("Cancer type (BC split by ER status)") +
      theme(legend.position = "right") + BOLD_LEG

ggsave(file.path(OUT_DIR, "01_umap_overview.pdf"), p1|p2, width=14, height=6)
ggsave(file.path(OUT_DIR, "01_umap_overview.png"), p1|p2, width=14, height=6, dpi=150)
cat("  Saved: 01_umap_overview\n")

# ── Fig 02: Cancer type bar ───────────────────────────────────────────────────
meta <- s@meta.data %>%
  group_by(mac_subcluster, cancer_type_er) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(mac_subcluster) %>% mutate(frac=n/sum(n))

col_order <- c("BC_ER+","BC_ER-","BC_ER?","BC",
               sort(setdiff(unique(meta$cancer_type_er), c("BC_ER+","BC_ER-","BC_ER?","BC"))))
meta$cancer_type_er <- factor(meta$cancer_type_er, levels=col_order)

fig2 <- ggplot(meta, aes(x=mac_subcluster, y=frac, fill=cancer_type_er)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values=CANCER_COLORS, name="Cancer type") +
  labs(x="Macrophage sub-cluster", y="Fraction",
       title="Cancer type composition per macrophage cluster",
       subtitle=sprintf("nn=%d  min.dist=%.2f  res=%.1f", BEST_NN, BEST_MD, BEST_RES)) +
  theme_classic(base_size=13) +
  theme(axis.text.x=element_text(angle=0)) + BOLD_LEG

ggsave(file.path(OUT_DIR, "02_cancer_type_bar_bc_er.pdf"), fig2, width=10, height=5)
ggsave(file.path(OUT_DIR, "02_cancer_type_bar_bc_er.png"), fig2, width=10, height=5, dpi=150)
cat("  Saved: 02_cancer_type_bar_bc_er\n")

# ── Fig 03a: Feature plots — highlight genes ──────────────────────────────────
fp_hl <- FeaturePlot(s, features=hl_genes, reduction="umap",
                     cols=c("lightgrey","#d62728"),
                     ncol=min(length(hl_genes),3), pt.size=0.2) &
         theme(legend.position="right") & BOLD_LEG
fw <- 6*min(length(hl_genes),3); fh <- 5*ceiling(length(hl_genes)/3)
ggsave(file.path(OUT_DIR,"03a_feature_plots_highlight.pdf"), fp_hl, width=fw, height=fh)
ggsave(file.path(OUT_DIR,"03a_feature_plots_highlight.png"), fp_hl, width=fw, height=fh, dpi=150)
cat("  Saved: 03a_feature_plots_highlight\n")

# ── Fig 03b: Feature plots — ESR1/SP1 ────────────────────────────────────────
fp_es <- FeaturePlot(s, features=esr1_sp1, reduction="umap",
                     cols=c("lightgrey","#2171b5"),
                     ncol=min(length(esr1_sp1),3), pt.size=0.2) &
         theme(legend.position="right") & BOLD_LEG
fw2 <- 6*min(length(esr1_sp1),3); fh2 <- 5*ceiling(length(esr1_sp1)/3)
ggsave(file.path(OUT_DIR,"03b_feature_plots_esr1_sp1.pdf"), fp_es, width=fw2, height=fh2)
ggsave(file.path(OUT_DIR,"03b_feature_plots_esr1_sp1.png"), fp_es, width=fw2, height=fh2, dpi=150)
cat("  Saved: 03b_feature_plots_esr1_sp1\n")

# ── Fig 03c: Combined dot plot ────────────────────────────────────────────────
dp_comb <- DotPlot(s, features=comb_genes, group.by="mac_subcluster") +
  coord_flip() +
  labs(title=paste("Gene expression per macrophage sub-cluster\n",
                   paste(comb_genes, collapse=" / "))) +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=0)) + BOLD_LEG
dp_h <- max(4, 0.35*length(comb_genes)+2)
ggsave(file.path(OUT_DIR,"03c_dotplot_combined.pdf"), dp_comb, width=10, height=dp_h)
ggsave(file.path(OUT_DIR,"03c_dotplot_combined.png"), dp_comb, width=10, height=dp_h, dpi=150)
cat("  Saved: 03c_dotplot_combined\n")

# ── Fig 04: BC highlight UMAP ─────────────────────────────────────────────────
coords <- as.data.frame(Embeddings(s,"umap"))
colnames(coords) <- c("UMAP1","UMAP2")
coords$cancer_type_er <- s$cancer_type_er
coords$highlight <- ifelse(coords$cancer_type_er %in% c("BC_ER+","BC_ER-","BC_ER?","BC"),
                            coords$cancer_type_er, "other")
coords$highlight <- factor(coords$highlight,
                            levels=c("BC_ER+","BC_ER-","BC_ER?","BC","other"))
hl_colors <- c("BC_ER+"="#d62728","BC_ER-"="#ff9896",
                "BC_ER?"="#fa9fb5","BC"="#fb6a4a","other"="grey85")

fig4 <- ggplot(coords %>% arrange(highlight=="other"),
               aes(UMAP1, UMAP2, color=highlight)) +
  geom_point(size=0.3, alpha=0.7) +
  scale_color_manual(values=hl_colors, name="BC status") +
  ggtitle("BC / BC_ER cells highlighted") +
  theme_classic(base_size=13) +
  guides(color=guide_legend(override.aes=list(size=3))) + BOLD_LEG

ggsave(file.path(OUT_DIR,"04_bc_highlight_umap.pdf"), fig4, width=8, height=6)
ggsave(file.path(OUT_DIR,"04_bc_highlight_umap.png"), fig4, width=8, height=6, dpi=150)
cat("  Saved: 04_bc_highlight_umap\n")

cat("\nDone — clustering figures regenerated with bold legends.\n")
