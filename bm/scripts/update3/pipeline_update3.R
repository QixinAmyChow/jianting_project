#!/usr/bin/env Rscript
# ============================================================
# BM Macrophage Sub-clustering Pipeline — Update 3
# Input:  Myeloid.rds  (GEO: GSE266330 Liu et al. 2025)
# Output: bm/figures/update3/  (all figures)
#         run build_slides.py afterwards for the pptx
#
# Checkpoint: if mac_update3_best.rds already exists and passes
#   the BC_ER+ validation (Cluster 1 = 1791/1889 cells), the
#   expensive clustering steps are skipped so figures are fully
#   reproducible.  Remove the checkpoint to re-cluster from scratch.
# ============================================================

suppressPackageStartupMessages({
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
  library(cowplot)
})
options(future.globals.maxSize = 2 * 1024^3)

# ── Paths ──────────────────────────────────────────────────────────────────────
MYELOID_RDS <- "bm/data/raw/geo/integrated_Seurat_objects/47.integrated_object_subset_by_major_celltypes/Myeloid.rds"
CHECKPOINT  <- "bm/data/rds/mac_update3_best.rds"
FALLBACK    <- "bm/data/mac_update3_final.rds"
HIGHLIGHT_F <- "bm/figures/update1/highlight_genes.txt"
PATIENT_XLS <- "bm/data/raw/geo/Table S1 patient infor_corrected_2025-08.xlsx"
OUT_DIR     <- "bm/figures/update3"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Shared constants ───────────────────────────────────────────────────────────
CANCER_COLORS <- c(
  "BC_ER+"="#d62728","BC_ER-"="#ff9896","BC_ER?"="#fa9fb5","BC"="#fb6a4a",
  "KC"="#8c564b","LC"="#17becf","CC"="#2ca02c","BDC"="#ff7f0e",
  "EC"="#9467bd","TC"="#7f7f7f","PC"="#bcbd22","ctrl"="#aec7e8"
)
BOLD_LEG <- theme(
  legend.title = element_text(face="bold", size=11),
  legend.text  = element_text(face="bold", size=10)
)
BEST_NN <- 10; BEST_MD <- 0.05; BEST_RES <- 0.3

# ── Semantic dedup families (shared across all GSEA scripts) ───────────────────
FAMILIES <- list(
  H_MYC="MYC_TARGET|HALLMARK_MYC", H_MTORC1="MTORC1",
  H_OXPHOS="OXIDATIVE_PHOSPH", H_GLYC="GLYCOLYSIS", H_HYPOXIA="HYPOXIA",
  H_APOP="APOPTOSIS", H_G2M="G2M_CHECKPOINT", H_E2F="E2F_TARGET",
  H_ANDROGEN="ANDROGEN", H_ESTROGEN="ESTROGEN", H_P53="P53|TP53",
  H_PI3K="PI3K|AKT_MTOR", H_COAGUL="COAGULATION", H_COMPLEM="COMPLEMENT",
  H_FATTY="FATTY_ACID", H_CHOLEST="CHOLESTEROL|BILE", H_ANGIO="ANGIOGENESIS",
  H_UNFOLDED="UNFOLDED_PROTEIN", H_DNA_REP="DNA_REPAIR",
  H_MITOTIC="MITOTIC_SPINDLE", H_UV="UV_RESPONSE",
  H_XENO="XENOBIOTIC", H_INFLAM="INFLAMMATORY_RESPONSE",
  C6_KRAS="KRAS", C6_PRC1="BMI1|MEL18|PRC1", C6_PRC2="PRC2|SUZ12|EZH2",
  C6_ESC="ESC_V6|ESC_J", C6_MORF="MORF_", C6_STK33="STK33",
  C6_AKT="^AKT", C6_CYCLIN="CYCLIN", C6_ERBB="ERBB2", C6_EGFR="EGFR",
  COLLAGEN="COLLAGEN",
  ECM_ORG="ECM_ORGAN|EXTRACELLULAR_MATRIX_ORGAN|EXTRACELLULAR_MATRIX_ASSEMB|ECM_ASSEMB|ECM_REMODEL|EXTRACELLULAR_STRUCT",
  MMP="MMP|MATRIX_METALLOPROTEIN|METALLOENDOPEPTI",
  BASEMENT_MEM="BASEMENT_MEMBRANE", INTEGRIN="INTEGRIN",
  FIBRONECTIN="FIBRONECT", LAMININ="LAMININ", ECM_GEN="EXTRACELLULAR_MATRIX",
  TGFB="TGFB|TGF.?B|TRANSFORMING_GROWTH|SMAD",
  EMT="EPITHELIAL_MESENCH|EMT_TRANS|MESENCHYMAL_TRANS|\\bEMT\\b",
  WNT="WNT|BETA_CATENIN",
  OSTEOBLAST="OSTEOBLAST|OSTEOGENIC_DIFF|BONE_FORM",
  OSTEOCLAST="OSTEOCLAST|BONE_RESORPT",
  OSSIFICATION="OSSIF|BONE_MINERAL|BONE_REMODEL|BONE_DEVELOP",
  BMP="\\bBMP\\b|BONE_MORPHOGEN", NOTCH="NOTCH", HEDGEHOG="HEDGEHOG",
  CD8_EXHAUST="CD8.*EXHAUST|EXHAUST.*CD8|T_CELL_EXHAUST",
  TREG="TREG|REGULATORY_T_CELL|T_REGULATORY",
  CHECKPOINT="CHECKPOINT|\\bPD.?1\\b|\\bPDL1\\b|CTLA4|LAG3|TIM3",
  IL10="\\bIL10\\b|\\bIL_10\\b",
  T_CELL_ACTIV="T_CELL_ACTIV|T_CELL_STIMUL|T_CELL_PROLIF|T_CELL_RECEPTOR_SIGNAL",
  CD8_GEN="\\bCD8\\b", CD4_GEN="\\bCD4\\b", T_CELL_GEN="T_CELL|TCELL",
  PPARG="PPARG|PPAR.?GAMMA", WBP7="WBP7",
  MDSC="MDSC|MYELOID_DERIVED_SUPPRESS",
  M2_POLAR="M2_POLARIZ|ALTERNATIVE_ACTIV|ALTERNATIVELY_ACTIV|ANTI.INFLAM.*MAC",
  M1_POLAR="M1_POLARIZ|CLASSIC_ACTIV|CLASSICALLY_ACTIV|INFLAM.*MAC",
  TAM="\\bTAM\\b|TUMOR.*MACRO|MACRO.*TUMOR|TUMOR_ASSOC.*MACRO",
  MYELOID_DIFF="MYELOID_DIFF|MYELOID_CELL_DIFF|MYELOID_LEUK_DIFF",
  MONOCYTE="MONOCYTE", MYELOID_GEN="MYELOID", MACRO_GEN="MACROPHAGE|\\bMACRO\\b",
  TNFA="\\bTNF\\b|TNF_ALPHA|TNFA", NFKB="NFKB|NF.KB|NFKAPPAB",
  JAK_STAT="JAK.?STAT|JAK_STAT", STAT3="STAT3", STAT_GEN="\\bSTAT\\b",
  IFN_GAMMA="INTERFERON_GAMMA|IFN_GAMMA|IFNG|TYPE_II_IFN",
  IFN_ALPHA="INTERFERON_ALPHA|IFN_ALPHA|IFNA|TYPE_I_IFN",
  IFN_GEN="INTERFERON|\\bIFN\\b",
  IL6="\\bIL6\\b|\\bIL_6\\b|INTERLEUKIN_6", IL4_IL13="\\bIL4\\b|\\bIL13\\b|IL_4|IL_13",
  CHEMOKINE="CHEMOKINE",
  CYTOKINE_PROD="CYTOKINE_PRODUCT|CYTOKINE_SECRET|CYTOKINE_BIOSYN",
  CYTOKINE_SIG="CYTOKINE_SIGNAL|CYTOKINE_MEDIAT|CYTOKINE_NETWORK",
  CYTOKINE_GEN="CYTOKINE"
)

# ── Helper functions ───────────────────────────────────────────────────────────
sem_dedup <- function(df, fams) {
  assign_fam <- function(pw) {
    for (nm in names(fams))
      if (grepl(fams[[nm]], pw, ignore.case=TRUE, perl=TRUE)) return(nm)
    paste0("UNIQUE__", pw)
  }
  df$sem_family <- vapply(df$pathway, assign_fam, character(1))
  df %>% group_by(sem_family) %>%
    slice_min(padj, n=1, with_ties=FALSE) %>%
    ungroup() %>% select(-sem_family)
}

shorten_label <- function(x, n=50) {
  x <- gsub("^HALLMARK_","H: ", x); x <- gsub("^GOBP_","", x)
  x <- gsub("^REACTOME_","R: ", x); x <- gsub("^GSE[0-9]+_","", x)
  x <- gsub("_"," ", x); x <- gsub("  +"," ", trimws(x))
  ifelse(nchar(x)>n, paste0(substr(x,1,n-3),"..."), x)
}

resolve_genes <- function(requested, available) {
  found <- intersect(requested, available)
  c(found, intersect(toupper(requested[!requested %in% found]), available))
}

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
  if (!is.null(fill_colors)) {
    p <- p + scale_fill_manual(values=fill_colors)
  } else {
    p <- p + scale_fill_brewer(palette="Set2")
  }
  ggsave(file.path(OUT_DIR, paste0(filename,".pdf")), p, width=10, height=5)
  ggsave(file.path(OUT_DIR, paste0(filename,".png")), p, width=10, height=5, dpi=150)
  cat("  Saved:", filename, "\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Clustering: Myeloid.rds → mac_update3_best.rds
# ═══════════════════════════════════════════════════════════════════════════════

checkpoint_ok <- function(path) {
  if (!file.exists(path)) return(FALSE)
  tryCatch({
    s <- readRDS(path)
    "mac_subcluster" %in% colnames(s@meta.data) &&
      sum(s$mac_subcluster == "1" & s$cancer_type_er == "BC_ER+",
          na.rm=TRUE) == 1791
  }, error=function(e) FALSE)
}

if (checkpoint_ok(CHECKPOINT)) {
  cat("=== Checkpoint found:", CHECKPOINT, "— loading ===\n")
  s <- readRDS(CHECKPOINT)

} else if (checkpoint_ok(FALLBACK)) {
  cat("=== Using fallback:", FALLBACK, "===\n")
  s <- readRDS(FALLBACK)
  saveRDS(s, CHECKPOINT)
  cat("  Saved as:", CHECKPOINT, "\n")

} else {
  cat("=== Section 1: Clustering from Myeloid.rds ===\n")

  # 1a. Load myeloid object and subset to Mφ
  cat("Loading Myeloid.rds ...\n")
  if (!file.exists(MYELOID_RDS)) stop("Myeloid.rds not found at: ", MYELOID_RDS)
  s_mye <- readRDS(MYELOID_RDS)
  cat("  Full myeloid:", ncol(s_mye), "cells\n")
  s <- subset(s_mye, celltype_C == "Mφ")
  rm(s_mye); gc()
  cat("  After Mφ subset:", ncol(s), "cells\n")

  # 1b. ER status from patient metadata
  cat("Reading patient metadata ...\n")
  pt <- read_excel(PATIENT_XLS, skip=1)
  colnames(pt) <- make.names(colnames(pt))
  er_map <- pt %>%
    filter(cancer == "BC") %>%
    select(cancer.id, cancer.subtype) %>% distinct() %>%
    mutate(er_status = case_when(
      grepl("^ER\\+",        cancer.subtype) ~ "BC_ER+",
      grepl("^ER-|^ER\\s*-", cancer.subtype) ~ "BC_ER-",
      TRUE ~ "BC_ER?"
    ))
  er_lookup <- setNames(er_map$er_status, er_map$cancer.id)
  s$cancer_type_er <- as.character(s$cancer)
  bc_cells <- s$cancer == "BC"
  s$cancer_type_er[bc_cells] <- er_lookup[s$cancer.id[bc_cells]]
  s$cancer_type_er[is.na(s$cancer_type_er)] <- "BC"
  cat("  cancer_type_er:\n"); print(table(s$cancer_type_er))

  # 1c. SCTransform + PCA
  cat("SCTransform ...\n")
  s <- SCTransform(s, vars.to.regress="percent.mt",
                   variable.features.n=2000, verbose=FALSE)
  cat("PCA ...\n")
  s <- RunPCA(s, npcs=50, verbose=FALSE)

  # 1d. Parameter sweep: n.neighbors × min.dist × resolution
  cat("Parameter sweep (72 combos) ...\n")
  nn_vals  <- c(10,20,30,50); md_vals <- c(0.05,0.1,0.3)
  res_vals <- c(0.3,0.5,0.8,1.0,1.2,1.5); n_pcs <- 1:20
  grid <- expand.grid(n_neighbors=nn_vals, min_dist=md_vals,
                      resolution=res_vals, stringsAsFactors=FALSE)

  score_sep <- function(ct, cl) {
    df <- data.frame(ct=ct, cl=cl)
    sapply(list(BC_ER_pos="BC_ER+", BC_family=c("BC_ER+","BC_ER-","BC_ER?","BC")),
           function(tgt) {
             n_tgt <- sum(df$ct %in% tgt)
             if (n_tgt==0) return(0)
             cs <- df %>% group_by(cl) %>%
               summarise(n_cl=n(), n_hit=sum(ct %in% tgt), .groups="drop") %>%
               mutate(prec=n_hit/n_cl, rec=n_hit/n_tgt)
             max(ifelse(cs$prec+cs$rec>0, 2*cs$prec*cs$rec/(cs$prec+cs$rec), 0))
           })
  }

  results <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    nn <- grid$n_neighbors[i]; md <- grid$min_dist[i]; res <- grid$resolution[i]
    s_tmp <- RunUMAP(s, dims=n_pcs, n.neighbors=nn, min.dist=md,
                     reduction.name="umap_tmp", reduction.key="UMAPtmp_", verbose=FALSE)
    s_tmp <- FindNeighbors(s_tmp, dims=n_pcs, verbose=FALSE)
    s_tmp <- FindClusters(s_tmp, resolution=res, verbose=FALSE, cluster.name="cl_tmp")
    sc <- score_sep(s_tmp$cancer_type_er, s_tmp$cl_tmp)
    results[[i]] <- list(tag=sprintf("nn%02d_md%.2f_res%.1f",nn,md,res),
                         n_neighbors=nn, min_dist=md, resolution=res,
                         n_clusters=length(unique(s_tmp$cl_tmp)),
                         f1_BC_ER_pos=sc["BC_ER_pos"], f1_BC_family=sc["BC_family"],
                         umap_coords=Embeddings(s_tmp,"umap_tmp"),
                         clusters=s_tmp$cl_tmp)
    if (i%%10==0||i==nrow(grid))
      cat(sprintf("  [%d/%d] %s  ncl=%d  F1_ER+=%.3f\n",i,nrow(grid),
                  results[[i]]$tag, results[[i]]$n_clusters, sc["BC_ER_pos"]))
  }

  score_df <- do.call(rbind, lapply(results, function(r)
    data.frame(tag=r$tag, n_neighbors=r$n_neighbors, min_dist=r$min_dist,
               resolution=r$resolution, n_clusters=r$n_clusters,
               f1_BC_ER_pos=r$f1_BC_ER_pos, f1_BC_family=r$f1_BC_family))) %>%
    arrange(desc(f1_BC_ER_pos))
  write.csv(score_df, file.path(OUT_DIR,"00_grid_scores.csv"), row.names=FALSE)
  cat("Top 5:\n"); print(head(score_df,5), row.names=FALSE)

  # Top-12 UMAP panel
  top12_idx <- match(head(score_df$tag,12), sapply(results,`[[`,"tag"))
  panel_plots <- lapply(top12_idx, function(i) {
    r <- results[[i]]
    coords <- as.data.frame(r$umap_coords); colnames(coords) <- c("UMAP1","UMAP2")
    coords$cancer_type_er <- s$cancer_type_er; coords$cluster <- r$clusters
    ggplot(coords, aes(UMAP1,UMAP2,color=cancer_type_er)) +
      geom_point(size=0.2, alpha=0.6) +
      scale_color_manual(values=CANCER_COLORS) +
      ggtitle(sprintf("%s\nncl=%d  F1=%.3f",r$tag,r$n_clusters,r$f1_BC_ER_pos)) +
      theme_classic(base_size=7) + theme(legend.position="none",plot.title=element_text(size=6))
  })
  ggsave(file.path(OUT_DIR,"00_top12_panel.pdf"), wrap_plots(panel_plots,ncol=4),width=20,height=12)
  ggsave(file.path(OUT_DIR,"00_top12_panel.png"), wrap_plots(panel_plots,ncol=4),width=20,height=12,dpi=120)

  # 1e. Final object with best params
  best <- results[[top12_idx[1]]]
  cat(sprintf("\nBest: %s  ncl=%d  F1_ER+=%.4f\n", best$tag, best$n_clusters, best$f1_BC_ER_pos))
  s <- RunUMAP(s, dims=n_pcs, n.neighbors=best$n_neighbors, min.dist=best$min_dist, verbose=FALSE)
  s <- FindNeighbors(s, dims=n_pcs, verbose=FALSE)
  s <- FindClusters(s, resolution=best$resolution, verbose=FALSE, cluster.name="mac_subcluster")
  Idents(s) <- "mac_subcluster"
  cat("  Clusters:", length(unique(s$mac_subcluster)), "\n")
  print(table(s$mac_subcluster, s$cancer_type_er))

  # Add treatment response
  tx_map <- pt %>%
    select(cancer.id,
           neo_resp=Neoadjuvant.treatment.response,
           met_resp=Metastatic...treatment.response) %>% distinct() %>%
    mutate(tx_response=case_when(
      !tolower(neo_resp) %in% c("na","") ~ tolower(neo_resp),
      !tolower(met_resp) %in% c("na","") ~ tolower(met_resp),
      TRUE ~ "na"))
  tx_lookup <- setNames(tx_map$tx_response, tx_map$cancer.id)
  s$tx_response <- unname(tx_lookup[s$cancer.id])
  s$tx_response[is.na(s$tx_response)] <- "na"
  s$gender_clean <- recode(s$gender, "f"="Female","m"="Male","na"="Unknown")

  saveRDS(s, CHECKPOINT)
  cat("  Saved:", CHECKPOINT, "\n")
}

# Fix cluster IDs: rename "11" → "0" to restore original 0-based labeling
cl_vals <- as.character(s$mac_subcluster)
if ("11" %in% cl_vals && !"0" %in% cl_vals) {
  cat("Renaming cluster 11 → 0 (restoring original labeling)\n")
  cl_vals[cl_vals == "11"] <- "0"
  s$mac_subcluster <- factor(cl_vals, levels=as.character(0:10))
  Idents(s) <- "mac_subcluster"
}
clusters <- as.character(sort(as.integer(levels(Idents(s)))))
cat("Active clusters:", paste(clusters, collapse=" "), "| Cells:", ncol(s), "\n")

# Add tx_response + gender_clean if missing (loaded from checkpoint)
if (!"tx_response" %in% colnames(s@meta.data)) {
  cat("Adding treatment response metadata ...\n")
  pt <- read_excel(PATIENT_XLS, skip=1); colnames(pt) <- make.names(colnames(pt))
  tx_map <- pt %>%
    select(cancer.id,
           neo_resp=Neoadjuvant.treatment.response,
           met_resp=Metastatic...treatment.response) %>% distinct() %>%
    mutate(tx_response=case_when(
      !tolower(neo_resp) %in% c("na","") ~ tolower(neo_resp),
      !tolower(met_resp) %in% c("na","") ~ tolower(met_resp),
      TRUE ~ "na"))
  tx_lookup <- setNames(tx_map$tx_response, tx_map$cancer.id)
  s$tx_response <- unname(tx_lookup[s$cancer.id])
  s$tx_response[is.na(s$tx_response)] <- "na"
}
if (!"gender_clean" %in% colnames(s@meta.data))
  s$gender_clean <- recode(as.character(s$gender), "f"="Female","m"="Male","na"="Unknown")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Sanity check figures
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Section 2: Sanity figures ===\n")
meta <- s@meta.data
ct_col <- if ("cancer_type_er" %in% colnames(meta)) "cancer_type_er" else "cancer"
meta$ct <- meta[[ct_col]]

ct_sum <- meta %>% group_by(ct) %>% summarise(n=n(),.groups="drop") %>%
  mutate(pct=100*n/sum(n), label=sprintf("%d\n(%.1f%%)",n,pct)) %>% arrange(desc(n))
ct_sum$ct <- factor(ct_sum$ct, levels=ct_sum$ct)

figA <- ggplot(ct_sum, aes(x=ct,y=n,fill=ct)) +
  geom_bar(stat="identity",width=0.75) +
  geom_text(aes(label=label),vjust=-0.3,size=3.2) +
  scale_fill_manual(values=CANCER_COLORS,guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(x="Cancer type",y="Cell count",title="Macrophage cells per cancer type",
       subtitle=sprintf("Total = %d macrophages across %d cancer types",nrow(meta),n_distinct(meta$ct))) +
  theme_classic(base_size=13) + BOLD_LEG

samp_sum <- meta %>% group_by(cancer.id, ct) %>% summarise(n=n(),.groups="drop") %>%
  arrange(ct, desc(n))
samp_sum$cancer.id <- factor(samp_sum$cancer.id, levels=unique(samp_sum$cancer.id))
bar_w <- max(8, 0.25*n_distinct(meta$cancer.id)+2)

figB <- ggplot(samp_sum, aes(x=cancer.id,y=n,fill=ct)) +
  geom_bar(stat="identity",width=0.8) +
  scale_fill_manual(values=CANCER_COLORS,name="Cancer type") +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(x="Sample (cancer.id)",y="Cell count",title="Macrophage cells per sample") +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8)) + BOLD_LEG

ggsave(file.path(OUT_DIR,"00_sanity_cells_per_cancer.pdf"),figA,width=9,height=5)
ggsave(file.path(OUT_DIR,"00_sanity_cells_per_cancer.png"),figA,width=9,height=5,dpi=150)
ggsave(file.path(OUT_DIR,"00_sanity_cells_per_sample.pdf"),figB,width=bar_w,height=5)
ggsave(file.path(OUT_DIR,"00_sanity_cells_per_sample.png"),figB,width=bar_w,height=5,dpi=150)
figC <- figA + figB + plot_layout(widths=c(1,2))
ggsave(file.path(OUT_DIR,"00_sanity_combined.pdf"),figC,width=bar_w+4,height=5)
ggsave(file.path(OUT_DIR,"00_sanity_combined.png"),figC,width=bar_w+4,height=5,dpi=150)
cat("  Saved: 00_sanity_*\n")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Clustering figures
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Section 3: Clustering figures ===\n")

# Resolve highlight genes
raw_hl <- if (file.exists(HIGHLIGHT_F)) {
  trimws(unlist(strsplit(paste(readLines(HIGHLIGHT_F), collapse=","), ",")))
} else {
  character(0)
}
raw_hl    <- raw_hl[nzchar(raw_hl)]
hl_genes  <- resolve_genes(raw_hl, rownames(s))
esr1_sp1  <- resolve_genes(c("ESR1","SP1"), rownames(s))
comb_genes <- unique(c(hl_genes, esr1_sp1))
cat("  Genes:", paste(comb_genes,collapse=", "), "\n")

# Fig 01: UMAP overview
p1 <- DimPlot(s, reduction="umap", group.by="mac_subcluster",
              label=TRUE, repel=TRUE, pt.size=0.3) +
      ggtitle(sprintf("Macrophage sub-clusters\n(nn=%d, md=%.2f, res=%.1f)",BEST_NN,BEST_MD,BEST_RES)) +
      theme(legend.position="right") + BOLD_LEG
p2 <- DimPlot(s, reduction="umap", group.by="cancer_type_er",
              cols=CANCER_COLORS, pt.size=0.3) +
      ggtitle("Cancer type (BC split by ER status)") +
      theme(legend.position="right") + BOLD_LEG
ggsave(file.path(OUT_DIR,"01_umap_overview.pdf"),p1|p2,width=14,height=6)
ggsave(file.path(OUT_DIR,"01_umap_overview.png"),p1|p2,width=14,height=6,dpi=150)
cat("  Saved: 01_umap_overview\n")

# Fig 02: Cancer type bar
m2 <- s@meta.data %>%
  group_by(mac_subcluster,cancer_type_er) %>% summarise(n=n(),.groups="drop") %>%
  group_by(mac_subcluster) %>% mutate(frac=n/sum(n))
col_ord <- c("BC_ER+","BC_ER-","BC_ER?","BC",
             sort(setdiff(unique(m2$cancer_type_er),c("BC_ER+","BC_ER-","BC_ER?","BC"))))
m2$cancer_type_er <- factor(m2$cancer_type_er, levels=col_ord)
fig2 <- ggplot(m2, aes(x=mac_subcluster,y=frac,fill=cancer_type_er)) +
  geom_bar(stat="identity",width=0.8) +
  scale_fill_manual(values=CANCER_COLORS,name="Cancer type") +
  labs(x="Macrophage sub-cluster",y="Fraction",
       title="Cancer type composition per macrophage cluster",
       subtitle=sprintf("nn=%d  min.dist=%.2f  res=%.1f",BEST_NN,BEST_MD,BEST_RES)) +
  theme_classic(base_size=13) + theme(axis.text.x=element_text(angle=0)) + BOLD_LEG
ggsave(file.path(OUT_DIR,"02_cancer_type_bar_bc_er.pdf"),fig2,width=10,height=5)
ggsave(file.path(OUT_DIR,"02_cancer_type_bar_bc_er.png"),fig2,width=10,height=5,dpi=150)
cat("  Saved: 02_cancer_type_bar_bc_er\n")

# Fig 03a: Feature plots — highlight genes
if (length(hl_genes)>0) {
  fp_hl <- FeaturePlot(s, features=hl_genes, reduction="umap",
                       cols=c("lightgrey","#d62728"),
                       ncol=min(length(hl_genes),3), pt.size=0.2) &
           theme(legend.position="right") & BOLD_LEG
  fw <- 6*min(length(hl_genes),3); fh <- 5*ceiling(length(hl_genes)/3)
  ggsave(file.path(OUT_DIR,"03a_feature_plots_highlight.pdf"),fp_hl,width=fw,height=fh)
  ggsave(file.path(OUT_DIR,"03a_feature_plots_highlight.png"),fp_hl,width=fw,height=fh,dpi=150)
  cat("  Saved: 03a_feature_plots_highlight\n")
}

# Fig 03b: Feature plots — ESR1/SP1
if (length(esr1_sp1)>0) {
  fp_es <- FeaturePlot(s, features=esr1_sp1, reduction="umap",
                       cols=c("lightgrey","#2171b5"),
                       ncol=min(length(esr1_sp1),3), pt.size=0.2) &
           theme(legend.position="right") & BOLD_LEG
  fw2 <- 6*min(length(esr1_sp1),3); fh2 <- 5*ceiling(length(esr1_sp1)/3)
  ggsave(file.path(OUT_DIR,"03b_feature_plots_esr1_sp1.pdf"),fp_es,width=fw2,height=fh2)
  ggsave(file.path(OUT_DIR,"03b_feature_plots_esr1_sp1.png"),fp_es,width=fw2,height=fh2,dpi=150)
  cat("  Saved: 03b_feature_plots_esr1_sp1\n")
}

# Fig 03c: Combined dot plot
dp_comb <- DotPlot(s, features=comb_genes, group.by="mac_subcluster") +
  coord_flip() +
  labs(title=paste("Gene expression per macrophage sub-cluster\n",paste(comb_genes,collapse=" / "))) +
  theme_classic(base_size=12) + theme(axis.text.x=element_text(angle=0)) + BOLD_LEG
dp_h <- max(4, 0.35*length(comb_genes)+2)
ggsave(file.path(OUT_DIR,"03c_dotplot_combined.pdf"),dp_comb,width=10,height=dp_h)
ggsave(file.path(OUT_DIR,"03c_dotplot_combined.png"),dp_comb,width=10,height=dp_h,dpi=150)
cat("  Saved: 03c_dotplot_combined\n")

# Fig 03b (dot): highlight genes dot plot
if (length(hl_genes)>0) {
  dp_hl <- DotPlot(s, features=hl_genes, group.by="mac_subcluster") +
    coord_flip() +
    labs(title=paste(paste(hl_genes,collapse="/"),"across macrophage sub-clusters")) +
    theme_classic(base_size=12) + theme(axis.text.x=element_text(angle=0)) + BOLD_LEG
  ggsave(file.path(OUT_DIR,"03b_dotplot_highlight.pdf"),dp_hl,width=10,height=dp_h)
  ggsave(file.path(OUT_DIR,"03b_dotplot_highlight.png"),dp_hl,width=10,height=dp_h,dpi=150)
  cat("  Saved: 03b_dotplot_highlight\n")
}

# Fig 04: BC highlight UMAP
coords <- as.data.frame(Embeddings(s,"umap")); colnames(coords) <- c("UMAP1","UMAP2")
coords$cancer_type_er <- s$cancer_type_er
coords$highlight <- ifelse(coords$cancer_type_er %in% c("BC_ER+","BC_ER-","BC_ER?","BC"),
                            coords$cancer_type_er, "other")
coords$highlight <- factor(coords$highlight, levels=c("BC_ER+","BC_ER-","BC_ER?","BC","other"))
hl_colors <- c("BC_ER+"="#d62728","BC_ER-"="#ff9896","BC_ER?"="#fa9fb5","BC"="#fb6a4a","other"="grey85")
fig4 <- ggplot(coords %>% arrange(highlight=="other"), aes(UMAP1,UMAP2,color=highlight)) +
  geom_point(size=0.3,alpha=0.7) +
  scale_color_manual(values=hl_colors,name="BC status") +
  ggtitle("BC / BC_ER cells highlighted") + theme_classic(base_size=13) +
  guides(color=guide_legend(override.aes=list(size=3))) + BOLD_LEG
ggsave(file.path(OUT_DIR,"04_bc_highlight_umap.pdf"),fig4,width=8,height=6)
ggsave(file.path(OUT_DIR,"04_bc_highlight_umap.png"),fig4,width=8,height=6,dpi=150)
cat("  Saved: 04_bc_highlight_umap\n")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — DEG + Hallmark/C6 GSEA + metadata figures
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Section 4: DEG + GSEA (Hallmark+C6) + metadata ===\n")

# DEG (cached)
deg_csv <- file.path(OUT_DIR,"deg_all_markers.csv")
if (file.exists(deg_csv)) {
  cat("Loading cached DEG ...\n")
  markers <- read.csv(deg_csv); markers$cluster <- as.character(markers$cluster)
} else {
  cat("Running FindAllMarkers ...\n")
  markers <- FindAllMarkers(s, only.pos=TRUE, min.pct=0.25,
                            logfc.threshold=0.25, verbose=FALSE) %>%
    filter(p_val_adj<0.05) %>% arrange(cluster, desc(avg_log2FC))
  write.csv(markers, deg_csv, row.names=FALSE)
}
cat("  Markers:", nrow(markers), "\n")

# Fig 04a: DEG heatmap
top5 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC,n=5) %>% ungroup()
hm <- DoHeatmap(s, features=unique(top5$gene), group.by="mac_subcluster",
                angle=0, size=3, draw.lines=TRUE) +
  scale_fill_gradientn(colors=rev(brewer.pal(11,"RdYlBu")),na.value="white") +
  theme(axis.text.y=element_text(size=7))
ggsave(file.path(OUT_DIR,"04a_deg_heatmap.pdf"),hm,width=14,height=8)
ggsave(file.path(OUT_DIR,"04a_deg_heatmap.png"),hm,width=14,height=8,dpi=150)
cat("  Saved: 04a_deg_heatmap\n")

# Fig 04b: DEG dot plot (top 3)
top3 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC,n=3) %>% ungroup()
dp_deg <- DotPlot(s, features=unique(top3$gene), group.by="mac_subcluster") +
  coord_flip() + labs(title="Top 3 DEGs per macrophage sub-cluster") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=0),axis.text.y=element_text(size=8)) + BOLD_LEG
ggsave(file.path(OUT_DIR,"04b_deg_dotplot.pdf"),dp_deg,width=12,height=8)
ggsave(file.path(OUT_DIR,"04b_deg_dotplot.png"),dp_deg,width=12,height=8,dpi=150)
cat("  Saved: 04b_deg_dotplot\n")

# Fig 04c: Combined dot plot — top 6 clusters
top6_ids <- as.character(sort(as.integer(clusters))[1:min(6,length(clusters))])
dp_comb_cl1 <- DotPlot(s, features=comb_genes, group.by="mac_subcluster",
                        idents=top6_ids) +
  coord_flip() +
  labs(title=paste("Highlight + ESR1/SP1 in top 6 clusters\n",paste(comb_genes,collapse=" / "))) +
  theme_classic(base_size=11) + theme(axis.text.x=element_text(angle=0)) + BOLD_LEG
ggsave(file.path(OUT_DIR,"04c_dotplot_comb_top_clusters.pdf"),dp_comb_cl1,width=10,height=dp_h)
ggsave(file.path(OUT_DIR,"04c_dotplot_comb_top_clusters.png"),dp_comb_cl1,width=10,height=dp_h,dpi=150)
cat("  Saved: 04c_dotplot_comb_top_clusters\n")

# GSEA Hallmark + C6 (cached)
gsea_csv <- file.path(OUT_DIR,"gsea_cancer_results.csv")
if (file.exists(gsea_csv)) {
  cat("Loading cached GSEA (Hallmark+C6) ...\n")
  gsea_all <- read.csv(gsea_csv); gsea_all$cluster <- as.character(gsea_all$cluster)
} else {
  cat("Running fgsea (Hallmark+C6) ...\n")
  msig_h  <- msigdbr(species="Homo sapiens", collection="H")
  msig_c6 <- msigdbr(species="Homo sapiens", collection="C6")
  pathways <- split(bind_rows(msig_h,msig_c6)$gene_symbol,
                    bind_rows(msig_h,msig_c6)$gs_name)
  gsea_results <- lapply(clusters, function(cl) {
    cat("  Cluster", cl, "...\n")
    m <- FindMarkers(s, ident.1=cl, features=VariableFeatures(s),
                     only.pos=FALSE, min.pct=0.1, logfc.threshold=0, verbose=FALSE)
    ranks <- sort(setNames(m$avg_log2FC, rownames(m)), decreasing=TRUE)
    res <- fgsea(pathways=pathways, stats=ranks, minSize=15, maxSize=500,
                 nPermSimple=1000, nproc=1)
    res$cluster <- cl; res
  })
  gsea_all <- bind_rows(gsea_results)
  write.csv(gsea_all %>% select(-leadingEdge), gsea_csv, row.names=FALSE)
}
cat("  GSEA rows:", nrow(gsea_all), "\n")

# Fig 05b: Bubble plot
top_gsea <- gsea_all %>% filter(padj<0.25) %>%
  group_by(cluster) %>% slice_max(abs(NES),n=3) %>% ungroup() %>%
  mutate(pathway_short=gsub("_"," ",gsub("^HALLMARK_","H: ",pathway)))
if (nrow(top_gsea)>0) {
  bp <- ggplot(top_gsea, aes(x=cluster,y=pathway_short,
                              size=-log10(pmax(padj,1e-4)),color=NES)) +
    geom_point(alpha=0.85) +
    scale_color_gradientn(colors=rev(brewer.pal(11,"RdBu")),limits=c(-3,3),
                          oob=scales::squish,name="NES") +
    scale_size_continuous(name="-log10(padj)",range=c(2,8)) +
    labs(x="Macrophage sub-cluster",y=NULL,
         title="Top GSEA pathways per cluster (Hallmark + Oncogenic C6, padj<0.25)") +
    theme_classic(base_size=12) + theme(axis.text.y=element_text(size=9)) + BOLD_LEG
  bub_h <- max(5, 0.3*length(unique(top_gsea$pathway_short))+2)
  ggsave(file.path(OUT_DIR,"05b_gsea_bubble.pdf"),bp,width=12,height=bub_h)
  ggsave(file.path(OUT_DIR,"05b_gsea_bubble.png"),bp,width=12,height=bub_h,dpi=150)
  cat("  Saved: 05b_gsea_bubble\n")
}

# Metadata bars
meta_bar(s@meta.data,"gender_clean","Gender distribution per macrophage sub-cluster","Female/Male/Unknown",
         c("Female"="#e377c2","Male"="#1f77b4","Unknown"="#aec7e8"),"06_gender_per_cluster")
meta_bar(s@meta.data,"tissue.origin","Tissue / bone site per macrophage sub-cluster","Site of bone metastasis",NULL,"07_tissue_origin_per_cluster")
meta_bar(s@meta.data,"tx_response","Treatment response per macrophage sub-cluster","Neoadjuvant or metastatic treatment response",NULL,"08_treatment_per_cluster")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — Broader GSEA (GO:BP + Reactome + C7-IMMUNESIGDB)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Section 5: Broader GSEA ===\n")

broader_csv <- file.path(OUT_DIR,"gsea_broader_results.csv")
if (file.exists(broader_csv)) {
  cat("Loading cached broader GSEA ...\n")
  gsea_broad <- read.csv(broader_csv); gsea_broad$cluster <- as.character(gsea_broad$cluster)
} else {
  cat("Running fgsea (GO:BP + Reactome + C7) ...\n")
  rank_lists <- lapply(clusters, function(cl) {
    m <- FindMarkers(s, ident.1=cl, features=VariableFeatures(s),
                     only.pos=FALSE, min.pct=0.1, logfc.threshold=0, verbose=FALSE)
    sort(setNames(m$avg_log2FC, rownames(m)), decreasing=TRUE)
  })
  names(rank_lists) <- clusters

  gobp     <- msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP")
  reactome <- msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:REACTOME")
  immune   <- msigdbr(species="Homo sapiens", collection="C7", subcollection="IMMUNESIGDB") %>%
    filter(grepl("MACRO|MONOCYTE|MYELOID|MDSC|DC_|DENDR|CD8|CD4|TREG|EXHAUST|POLARIZ|M1_|M2_|TAM|TUMOR_INFILT|NK_|IFN|TGFB|IL10|IL6|PD1|PDL|CHECKPOINT|BREAST|BONE|METAST",
                 gs_name, ignore.case=TRUE, perl=TRUE))
  pathways_broad <- list(
    GOBP     = split(gobp$gene_symbol,     gobp$gs_name),
    REACTOME = split(reactome$gene_symbol, reactome$gs_name),
    IMMUNE   = split(immune$gene_symbol,   immune$gs_name)
  )
  cat("  GO:BP:",length(pathways_broad$GOBP),
      "| Reactome:",length(pathways_broad$REACTOME),
      "| C7:",length(pathways_broad$IMMUNE),"\n")

  res_list <- list()
  for (db in names(pathways_broad)) {
    cat("  DB:", db, "\n")
    for (cl in clusters) {
      cat("    Cluster", cl, "...\n")
      res <- tryCatch(
        fgsea(pathways=pathways_broad[[db]], stats=rank_lists[[cl]],
              minSize=15, maxSize=500, nPermSimple=1000, nproc=1),
        error=function(e) NULL)
      if (!is.null(res)) {
        res$cluster <- cl; res$database <- db
        res_list[[paste(db,cl,sep="_")]] <- res
      }
    }
  }
  gsea_broad <- bind_rows(res_list) %>% select(-leadingEdge)
  write.csv(gsea_broad, broader_csv, row.names=FALSE)
}
cat("  Broader GSEA rows:", nrow(gsea_broad), "\n")

# Fig 05c: Cluster 1 hypothesis GSEA
hypo_kw <- list(
  "ECM remodeling"             = "COLLAGEN|EXTRACELLULAR_MATRIX|ECM|MMP|MATRIX_METAL|FIBRONECT|LAMININ|INTEGRIN|ECM_REMODEL|BASEMENT_MEMBRANE",
  "TGF-b / EMT"                = "TGF.?B|TRANSFORMING_GROWTH|EPITHELIAL_MESENCH|EMT|SMAD|TGFB",
  "Wnt / osteogenic niche"     = "WNT|OSTEOBLAST|OSTEOCLAST|BONE_REMODEL|NOTCH|HEDGEHOG|BMP|OSSIF",
  "Immunosuppression / T cell" = "T_CELL|TCELL|CD8|CD4|TREG|EXHAUST|CHECKPOINT|PD.?1|PDL|CTLA|LAG3|TIM3|IMMUNE_EVASI|IMMUNO_SUPPRESS|IL10|MDSC",
  "Macrophage polarization"    = "MACRO|MACROPHAGE|\\bM1\\b|\\bM2\\b|\\bTAM\\b|POLARIZ|MYELOID|MONOCYTE|ALTERNATIVE_ACTIV|CLASSIC_ACTIV|PPARG|WBP7",
  "Cytokine / NFkB / JAK-STAT" = "CYTOKINE|NFKB|NF.KB|JAK.?STAT|STAT3|IL6|IL4|IL13|INTERFERON|CHEMOKINE|TNF"
)

gsea_hypo <- bind_rows(lapply(names(hypo_kw), function(h)
  gsea_broad %>% filter(grepl(hypo_kw[[h]],pathway,ignore.case=TRUE,perl=TRUE),!is.na(NES)) %>%
    mutate(hypothesis=h)))

cl1_hypo <- gsea_hypo %>% filter(cluster=="1") %>%
  mutate(path_short=shorten_label(pathway,n=55)) %>%
  group_by(hypothesis) %>% slice_max(abs(NES),n=5) %>% ungroup() %>%
  arrange(hypothesis,NES) %>%
  mutate(path_short=factor(path_short,levels=unique(path_short)),
         dir=ifelse(NES>0,"Up-reg","Down-reg"))

fig_c1_dir <- ggplot(cl1_hypo, aes(x=NES,y=path_short,fill=dir)) +
  geom_bar(stat="identity",width=0.75) +
  geom_vline(xintercept=0,linewidth=0.4,color="grey40") +
  scale_fill_manual(values=c("Up-reg"="#d62728","Down-reg"="#2171b5"),name=NULL) +
  facet_wrap(~hypothesis,scales="free_y",ncol=2) +
  labs(title="Cluster 1 (94.8% BC_ER+) — Hypothesis-Relevant GSEA Pathways",
       subtitle="GO:BP + Reactome + C7-IMMUNESIGDB | top 5 by |NES| | red=up, blue=down",
       x="NES",y=NULL) +
  theme_bw(base_size=10) +
  theme(strip.text=element_text(face="bold",size=9),
        axis.text.y=element_text(size=7.5),
        plot.title=element_text(face="bold",size=12),
        legend.position="bottom") + BOLD_LEG
ggsave(file.path(OUT_DIR,"05c_cluster1_hypo_gsea.pdf"),fig_c1_dir,width=16,height=12)
ggsave(file.path(OUT_DIR,"05c_cluster1_hypo_gsea.png"),fig_c1_dir,width=16,height=12,dpi=150)

db_colors <- c("GOBP"="#2ca02c","REACTOME"="#ff7f0e","IMMUNE"="#9467bd")

# Semantic-dedup version
broad_broad_fams <- list(
  ECM_ORG="ECM_ORGANIZATION|EXTRACELLULAR_MATRIX_ORGAN|EXTRACELLULAR_MATRIX_ASSEMBLY",
  COLLAGEN_FORM="COLLAGEN_FIBRIL|COLLAGEN_BIOSYN",COLLAGEN_GEN="COLLAGEN",
  MMP="MMP|MATRIX_METALLOPROTEIN",INTEGRIN="INTEGRIN",FIBRONECT="FIBRONECT",
  LAMININ="LAMININ",BASEMENT_MEM="BASEMENT_MEMBRANE",
  TGFB_SIG="TGFB_SIGNALING|TGF_BETA_SIGNAL|SIGNALING_BY_TGF|TGF_BETA_RECEPT",
  TGFB_GEN="TGFB|TGF.?B",EMT_TRANS="EPITHELIAL_MESENCH|EMT_TRANS",EMT_GEN="EMT",SMAD="SMAD",
  WNT_CANON="CANONICAL_WNT|WNT_CANONICAL|BETA_CATENIN|WNT_LIGAND",WNT_GEN="WNT",
  OSTEOBLAST="OSTEOBLAST|OSTEOGENIC_DIFF",OSTEOCLAST="OSTEOCLAST|BONE_RESORPT",
  BONE_REMODEL="BONE_REMODEL|BONE_MINERAL|OSSIFICATION",BMP="BMP",NOTCH="NOTCH",HEDGEHOG="HEDGEHOG",
  CD8_EXHAUST="CD8.*EXHAUST|EXHAUST.*CD8|T_CELL_EXHAUST",CD8_GEN="CD8",
  TREG="TREG|REGULATORY_T_CELL|T_REG",CHECKPOINT="CHECKPOINT|PD.?1|PDL1|CTLA4|LAG3|TIM3",
  IL10="IL10|IL_10",T_CELL_ACTIV="T_CELL_ACTIV|T_CELL_STIMUL|T_CELL_PROLIF",T_CELL_GEN="T_CELL",
  M2_POLARIZ="M2_POLARIZ|ALTERNATIVE_ACTIV|ANTI.INFLAM.*MAC",
  M1_POLARIZ="M1_POLARIZ|CLASSIC_ACTIV|INFLAM.*MAC",
  TAM="TAM|TUMOR.*MACRO|MACRO.*TUMOR",MYELOID_DIFF="MYELOID_DIFF|MYELOID_CELL_DIFF",
  MYELOID_GEN="MYELOID",MONOCYTE="MONOCYTE",
  NFKB="NFKB|NF.KB",JAK_STAT="JAK.?STAT|STAT3_TARGET",STAT_GEN="STAT",
  IFN_GAMMA="INTERFERON_GAMMA|IFN_GAMMA",IFN_ALPHA="INTERFERON_ALPHA|IFN_ALPHA|TYPE_I_IFN",
  IFN_GEN="INTERFERON",CHEMOKINE="CHEMOKINE",TNF="TNF",IL6="IL6|IL_6",
  IL4_IL13="IL4|IL13|IL_4|IL_13",
  CYTOKINE_PROD="CYTOKINE_PRODUCT|CYTOKINE_SECRET",CYTOKINE_GEN="CYTOKINE"
)

cl1_db <- gsea_hypo %>% filter(cluster=="1") %>%
  group_by(hypothesis) %>%
  group_modify(~{
    pool_dd <- sem_dedup(.x %>% slice_max(abs(NES),n=15,with_ties=FALSE), broad_broad_fams)
    pool_dd %>% slice_min(padj,n=5,with_ties=FALSE)
  }) %>% ungroup() %>% arrange(hypothesis,NES) %>%
  mutate(path_short=shorten_label(pathway,n=55),
         path_short=factor(path_short,levels=unique(path_short)),
         dir=ifelse(NES>0,"Up-reg","Down-reg"))

fig_c1_db <- ggplot(cl1_db, aes(x=NES,y=path_short,fill=database)) +
  geom_bar(stat="identity",width=0.75) +
  geom_vline(xintercept=0,linewidth=0.4,color="grey40") +
  scale_fill_manual(values=db_colors,name="Database") +
  facet_wrap(~hypothesis,scales="free_y",ncol=2) +
  labs(title="Cluster 1 (94.8% BC_ER+) — Hypothesis-Relevant GSEA Pathways",
       subtitle="GO:BP + Reactome + C7-IMMUNESIGDB | top 5 by padj (semantically deduplicated)",
       x="NES",y=NULL) +
  theme_bw(base_size=10) +
  theme(strip.text=element_text(face="bold",size=9),
        axis.text.y=element_text(size=7.5),
        plot.title=element_text(face="bold",size=12),
        legend.position="bottom") + BOLD_LEG
ggsave(file.path(OUT_DIR,"05c_cluster1_db_gsea.pdf"),fig_c1_db,width=16,height=12)
ggsave(file.path(OUT_DIR,"05c_cluster1_db_gsea.png"),fig_c1_db,width=16,height=12,dpi=150)
cat("  Saved: 05c_cluster1_hypo_gsea + 05c_cluster1_db_gsea\n")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — All-database stacking figure (05e)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Section 6: All-DB stacking ===\n")

gsea_cancer <- read.csv(gsea_csv)
gsea_cancer$cluster  <- as.character(gsea_cancer$cluster)
gsea_cancer$database <- ifelse(grepl("^HALLMARK_",gsea_cancer$pathway),"HALLMARK","C6")
gsea_combined <- bind_rows(gsea_cancer, gsea_broad) %>% filter(!is.na(NES))

all_clusters <- as.character(sort(as.integer(unique(gsea_combined$cluster))))
top_df <- bind_rows(lapply(all_clusters, function(cl) {
  d <- gsea_combined %>% filter(cluster==cl)
  up5 <- sem_dedup(d %>% filter(NES>0) %>% slice_max(abs(NES),n=20,with_ties=FALSE),FAMILIES) %>%
    slice_min(padj,n=5,with_ties=FALSE)
  dn5 <- sem_dedup(d %>% filter(NES<0) %>% slice_max(abs(NES),n=20,with_ties=FALSE),FAMILIES) %>%
    slice_min(padj,n=5,with_ties=FALSE)
  bind_rows(dn5,up5) %>% mutate(cluster=cl)
})) %>% mutate(path_short=shorten_label(pathway,n=42),
               direction=ifelse(NES>0,"Up","Down"),
               cluster_lbl=paste0("Cluster ",cluster))

DB_COLORS <- c("HALLMARK"="#1f77b4","C6"="#aec7e8","GOBP"="#2ca02c",
               "REACTOME"="#ff7f0e","IMMUNE"="#9467bd")
cl_ord <- as.character(sort(as.integer(unique(top_df$cluster))))

plot_list <- lapply(cl_ord, function(cl) {
  df <- top_df %>% filter(cluster==cl) %>% arrange(NES) %>%
    mutate(path_short=factor(path_short,levels=unique(path_short)))
  ggplot(df, aes(x=NES,y=path_short,fill=database)) +
    geom_bar(stat="identity",width=0.72) +
    geom_vline(xintercept=0,linewidth=0.3,color="grey40") +
    scale_fill_manual(values=DB_COLORS,name="Database",drop=FALSE) +
    scale_x_continuous(expand=expansion(mult=c(0.05,0.05))) +
    labs(title=paste0("Cluster ",cl),x="NES",y=NULL) +
    theme_classic(base_size=7) +
    theme(axis.text.y=element_text(size=6),axis.text.x=element_text(size=6),
          axis.title.x=element_text(size=6.5),
          plot.title=element_text(face="bold",size=8.5,hjust=0),
          legend.position="none")
})

leg_dummy <- ggplot(top_df,aes(x=NES,y=path_short,fill=database)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=DB_COLORS,name="Database",
                    guide=guide_legend(title.position="top",override.aes=list(size=5))) +
  theme_classic(base_size=10) +
  theme(legend.position="bottom",
        legend.title=element_text(face="bold",size=10),
        legend.text=element_text(face="bold",size=9))
shared_legend <- get_legend(leg_dummy)

n_empty <- 12 - length(plot_list)
for (i in seq_len(n_empty)) plot_list[[length(plot_list)+1]] <- plot_spacer()

grid_fig <- wrap_plots(plot_list,nrow=3,ncol=4) +
  plot_annotation(
    title="GSEA — Top 5 Up/Down per Cluster | All Databases | Deduplicated",
    subtitle="HALLMARK (blue) · C6 (lt.blue) · GO:BP (green) · Reactome (orange) · C7-IMMUNESIGDB (purple)  |  semantic dedup · top 5 by padj",
    theme=theme(plot.title=element_text(face="bold",size=13),plot.subtitle=element_text(size=9))
  )
fig_stack <- plot_grid(grid_fig, shared_legend, ncol=1, rel_heights=c(20,1))
ggsave(file.path(OUT_DIR,"05e_gsea_alldb_stacking.pdf"),fig_stack,width=20,height=15)
ggsave(file.path(OUT_DIR,"05e_gsea_alldb_stacking.png"),fig_stack,width=20,height=15,dpi=150)
cat("  Saved: 05e_gsea_alldb_stacking\n")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — Cluster 1 combined GSEA (05f)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Section 7: Cluster 1 combined GSEA ===\n")

cl1_all <- gsea_combined %>% filter(cluster=="1")
selected <- bind_rows(lapply(names(hypo_kw), function(h) {
  d <- cl1_all %>% filter(grepl(hypo_kw[[h]],pathway,ignore.case=TRUE,perl=TRUE))
  if (nrow(d)==0) return(NULL)
  sem_dedup(d %>% slice_max(abs(NES),n=30,with_ties=FALSE),FAMILIES) %>%
    slice_min(padj,n=5,with_ties=FALSE) %>% mutate(hypothesis=h)
})) %>% arrange(hypothesis,NES) %>%
  mutate(path_short=shorten_label(pathway,n=50),
         path_short=factor(path_short,levels=unique(path_short)),
         sig_label=ifelse(padj<0.05,"*",ifelse(padj<0.1,".","")))

cat("  Terms per hypothesis:\n"); print(selected %>% count(hypothesis), row.names=FALSE)

group_info <- selected %>% group_by(hypothesis) %>%
  summarise(y_min=min(as.integer(path_short)),y_max=max(as.integer(path_short)),.groups="drop") %>%
  mutate(y_mid=(y_min+y_max)/2)
x_range <- max(abs(selected$NES),na.rm=TRUE)*1.12

fig_cl1f <- ggplot(selected, aes(x=NES,y=path_short,fill=database)) +
  geom_bar(stat="identity",width=0.72,color="white",linewidth=0.15) +
  geom_vline(xintercept=0,linewidth=0.45,color="grey30") +
  geom_text(aes(label=sig_label,x=NES+ifelse(NES>=0,0.05,-0.05)),
            hjust=ifelse(selected$NES>=0,0,1),size=3.5,color="grey20") +
  scale_fill_manual(values=DB_COLORS,name="Database") +
  scale_x_continuous(limits=c(-x_range,x_range*1.15),expand=expansion(0)) +
  labs(x="NES",y=NULL) +
  theme_classic(base_size=11) +
  theme(axis.text.y=element_text(size=9),axis.text.x=element_text(size=9),
        plot.title=element_text(face="bold",size=13),
        legend.position="bottom",
        legend.title=element_text(face="bold",size=10),
        legend.text=element_text(face="bold",size=9))
fig_h <- max(8, 0.3*nlevels(selected$path_short)+3.5)
ggsave(file.path(OUT_DIR,"05f_cl1_combined_nofacet.pdf"),fig_cl1f,width=14,height=fig_h,limitsize=FALSE)
ggsave(file.path(OUT_DIR,"05f_cl1_combined_nofacet.png"),fig_cl1f,width=14,height=fig_h,dpi=150,limitsize=FALSE)
cat("  Saved: 05f_cl1_combined_nofacet\n")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 8 — Save GSEA selections (slide 12 + slide 14)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Section 8: Save GSEA selections ===\n")

slide12 <- bind_rows(lapply(all_clusters, function(cl) {
  d <- gsea_combined %>% filter(cluster==cl)
  up5 <- sem_dedup(d %>% filter(NES>0) %>% slice_max(abs(NES),n=20,with_ties=FALSE),FAMILIES) %>%
    slice_min(padj,n=5,with_ties=FALSE)
  dn5 <- sem_dedup(d %>% filter(NES<0) %>% slice_max(abs(NES),n=20,with_ties=FALSE),FAMILIES) %>%
    slice_min(padj,n=5,with_ties=FALSE)
  bind_rows(dn5,up5) %>% mutate(cluster=cl)
})) %>% mutate(direction=ifelse(NES>0,"Up","Down"))

slide14 <- bind_rows(lapply(names(hypo_kw), function(h) {
  d <- cl1_all %>% filter(grepl(hypo_kw[[h]],pathway,ignore.case=TRUE,perl=TRUE))
  if (nrow(d)==0) return(NULL)
  sem_dedup(d %>% slice_max(abs(NES),n=30,with_ties=FALSE),FAMILIES) %>%
    slice_min(padj,n=5,with_ties=FALSE) %>% mutate(hypothesis=h)
})) %>% mutate(direction=ifelse(NES>0,"Up","Down"))

saveRDS(list(slide12=slide12,slide14=slide14,hypo_keywords=hypo_kw,families=FAMILIES),
        file.path(OUT_DIR,"gsea_slide12_14_selections.rds"))
write.csv(slide12, file.path(OUT_DIR,"gsea_slide12_selection.csv"), row.names=FALSE)
write.csv(slide14, file.path(OUT_DIR,"gsea_slide14_cl1_selection.csv"), row.names=FALSE)
cat(sprintf("  slide12 rows: %d | slide14 rows: %d\n", nrow(slide12), nrow(slide14)))

cat("\n===== pipeline_update3.R complete =====\n")
cat("All figures in:", OUT_DIR, "\n")
cat("Run build_slides.py to generate the pptx\n")
