#!/usr/bin/env Rscript
# Broader GSEA — 3 databases relevant to BM manuscript hypotheses
#
# Manuscript hypotheses:
#   1. ER+ macrophage → CD8+ T cell suppression (S2C2 crosstalk)
#   2. ECM remodeling: Timp2/Mmp14/Pdgfa/Fn1/Vim in niche transition
#   3. TGF-β / Wnt / PD-1 immunosuppression (W3/metastatic niche)
#   4. MSC–DTC co-migration (perivascular → osteogenic)
#   5. SP1 + ESR1 co-expression as macrophage effector
#
# Databases:
#   A. GO:BP       — ECM organization, immune process, cytokine signaling
#   B. Reactome    — TGF-β, Wnt, ECM, detailed immune pathways
#   C. C7 IMMUNESIGDB — macrophage M1/M2, T cell exhaustion, cancer immune

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(fgsea)
library(msigdbr)
library(RColorBrewer)

options(future.globals.maxSize = 2 * 1024^3)

RDS_IN  <- "bm/data/rds/mac_update3_best.rds"
OUT_DIR <- "bm/figures/update3"
CACHE   <- file.path(OUT_DIR, "gsea_broader_results.csv")

BOLD_LEG <- theme(
  legend.title = element_text(face = "bold", size = 11),
  legend.text  = element_text(face = "bold", size = 10)
)

shorten <- function(x, n = 55) {
  x <- gsub("^GOBP_|^REACTOME_|^GSE[0-9]+_", "", x)
  x <- gsub("_", " ", x)
  x <- gsub("  +", " ", trimws(x))
  ifelse(nchar(x) > n, paste0(substr(x, 1, n-3), "..."), x)
}

# ── 1. Load Seurat + build per-cluster rank lists ─────────────────────────────
cat("Loading Seurat object ...\n")
s <- readRDS(RDS_IN)
Idents(s) <- "mac_subcluster"
clusters  <- as.character(sort(as.integer(levels(Idents(s)))))

rank_lists <- lapply(clusters, function(cl) {
  m <- FindMarkers(s, ident.1 = cl,
                   features        = VariableFeatures(s),
                   only.pos        = FALSE,
                   min.pct         = 0.1,
                   logfc.threshold = 0,
                   verbose         = FALSE)
  sort(setNames(m$avg_log2FC, rownames(m)), decreasing = TRUE)
})
names(rank_lists) <- clusters
cat("  Rank lists built for clusters:", paste(clusters, collapse=" "), "\n")

# ── 2. Load gene sets ─────────────────────────────────────────────────────────
cat("\nLoading gene sets ...\n")

gobp     <- msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP")
reactome <- msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:REACTOME")
# C7 IMMUNESIGDB — filter to macrophage/monocyte/DC/T cell keywords for relevance + speed
immune   <- msigdbr(species="Homo sapiens", collection="C7", subcollection="IMMUNESIGDB")
immune   <- immune %>%
  filter(grepl("MACRO|MONOCYTE|MYELOID|MDSC|DC_|DENDR|CD8|CD4|TREG|EXHAUST|POLARIZ|
                M1_|M2_|TAM|TUMOR_INFILT|NK_|IFN|TGFB|IL10|IL6|PD1|PDL|CHECKPOINT|
                BREAST|BONE|METAST",
               gs_name, ignore.case = TRUE, perl = TRUE))

pathways <- list(
  GOBP     = split(gobp$gene_symbol,     gobp$gs_name),
  REACTOME = split(reactome$gene_symbol, reactome$gs_name),
  IMMUNE   = split(immune$gene_symbol,   immune$gs_name)
)
cat("  GO:BP:", length(pathways$GOBP),
    "| Reactome:", length(pathways$REACTOME),
    "| C7-immune:", length(pathways$IMMUNE), "\n\n")

# ── 3. Run GSEA (or load cache) ───────────────────────────────────────────────
if (file.exists(CACHE)) {
  cat("Loading cached broader GSEA results ...\n")
  gsea_all <- read.csv(CACHE)
  gsea_all$cluster <- as.character(gsea_all$cluster)
  cat("  Rows:", nrow(gsea_all), "\n")
} else {
  cat("Running fgsea across 3 databases × 11 clusters ...\n")
  results <- list()
  for (db in names(pathways)) {
    pw <- pathways[[db]]
    cat("  Database:", db, "(", length(pw), "gene sets)\n")
    for (cl in clusters) {
      cat("    Cluster", cl, "...\n")
      res <- tryCatch(
        fgsea(pathways=pw, stats=rank_lists[[cl]],
              minSize=15, maxSize=500, nPermSimple=1000, nproc=1),
        error = function(e) { cat("    ERROR:", conditionMessage(e), "\n"); NULL }
      )
      if (!is.null(res)) {
        res$cluster  <- cl
        res$database <- db
        results[[paste(db, cl, sep="_")]] <- res
      }
    }
  }
  gsea_all <- bind_rows(results) %>% select(-leadingEdge)
  write.csv(gsea_all, CACHE, row.names = FALSE)
  cat("  Saved cache:", CACHE, "\n")
}

# ── 4. BM-hypothesis keyword filters ─────────────────────────────────────────
# Group pathways by manuscript hypothesis
hypo_keywords <- list(
  "ECM remodeling"          = "COLLAGEN|EXTRACELLULAR_MATRIX|ECM|MMP|MATRIX_METAL|
                               FIBRONECT|LAMININ|INTEGRIN|ECM_REMODEL|BASEMENT_MEMBRANE",
  "TGF-b / EMT"             = "TGF.?B|TRANSFORMING_GROWTH|EPITHELIAL_MESENCH|EMT|
                               SMAD|TGFB",
  "Wnt / osteogenic niche"  = "WNT|OSTEOBLAST|OSTEOCLAST|BONE_REMODEL|
                               NOTCH|HEDGEHOG|BMP",
  "Immunosuppression / T cell crosstalk" = "T_CELL|CD8|CD4|TREG|EXHAUST|
                               CHECKPOINT|PD.?1|PDL|CTLA|LAG3|TIM3|
                               IMMUNE_EVASI|IMMUNO_SUPPRESS|IL10|TGFB.*IMMUNE",
  "Macrophage polarization"  = "MACRO|M1|M2|TAM|POLARIZ|MYELOID|MONOCYTE|
                               ALTERNATIVE_ACTIV|CLASSIC_ACTIV|INFLAMMATORY_MAC|
                               ANTI.INFLAM",
  "Cytokine / NFkB / JAK-STAT" = "CYTOKINE|NFKB|NF.KB|JAK.?STAT|STAT3|IL6|
                               IL4|IL13|INTERFERON|CHEMOKINE|TNF"
)

gsea_hypo <- bind_rows(lapply(names(hypo_keywords), function(h) {
  kw  <- gsub("\\s+","",hypo_keywords[[h]])  # collapse whitespace in pattern
  gsea_all %>%
    filter(grepl(kw, pathway, ignore.case=TRUE, perl=TRUE),
           !is.na(NES)) %>%
    mutate(hypothesis = h)
}))

cat("\nPathways per hypothesis group:\n")
print(gsea_hypo %>% group_by(hypothesis) %>% summarise(n=n_distinct(pathway)), row.names=FALSE)

# ── 5. Figure A: Cluster 1 (BC_ER+) spotlight ────────────────────────────────
cat("\n--- Cluster 1 (BC_ER+) hypothesis-relevant hits ---\n")
cl1 <- gsea_hypo %>%
  filter(cluster == "1") %>%
  mutate(path_short = shorten(pathway)) %>%
  group_by(hypothesis) %>%
  slice_max(abs(NES), n = 5) %>%
  ungroup() %>%
  arrange(hypothesis, NES) %>%
  mutate(path_short = factor(path_short, levels = unique(path_short)),
         dir = ifelse(NES > 0, "Up-reg", "Down-reg"))

print(cl1 %>% select(hypothesis, path_short, NES, padj, database) %>%
        arrange(hypothesis, NES), row.names=FALSE)

fig_cl1 <- ggplot(cl1, aes(x=NES, y=path_short, fill=dir)) +
  geom_bar(stat="identity", width=0.75) +
  geom_vline(xintercept=0, linewidth=0.4, color="grey40") +
  scale_fill_manual(values=c("Up-reg"="#d62728","Down-reg"="#2171b5"), name=NULL) +
  facet_wrap(~hypothesis, scales="free_y", ncol=2) +
  labs(title="Cluster 1 (94.8% BC_ER+) — Hypothesis-Relevant GSEA Pathways",
       subtitle="GO:BP + Reactome + C7-IMMUNESIGDB | top 5 per hypothesis group | ranked by NES",
       x="NES", y=NULL) +
  theme_bw(base_size=10) +
  theme(strip.text      = element_text(face="bold", size=9),
        axis.text.y     = element_text(size=7.5),
        plot.title      = element_text(face="bold", size=12),
        plot.subtitle   = element_text(size=9),
        legend.position = "bottom") + BOLD_LEG

ggsave(file.path(OUT_DIR,"05c_cluster1_hypo_gsea.pdf"), fig_cl1, width=16, height=12)
ggsave(file.path(OUT_DIR,"05c_cluster1_hypo_gsea.png"), fig_cl1, width=16, height=12, dpi=150)
cat("  Saved: 05c_cluster1_hypo_gsea\n")


cat("\n===== gsea_broader_update3.R complete =====\n")
cat("Figures: 05c_cluster1_hypo_gsea\n")
cat("Cache:   gsea_broader_results.csv\n")
