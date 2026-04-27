#!/usr/bin/env Rscript
# annotate_integrated_obj.R
# Loads Myeloid.rds (home) + T_NK.rds (/tmp) from the integrated object subsets.
# Adds two metadata columns:
#   col1  = cell_group:    BC macrophage | other cancer macrophage |
#                          BC T cell     | other cancer T cell     | other
#   col2  = bc_er_status:  ER+ | ER- | ER? | NA (non-BC cells)
# Merges and saves integrated_annotated.rds to jianting_project/

library(Seurat)
library(dplyr)

MYELOID_RDS <- "bm_analysis_out/raw_geo/integrated_Seurat_objects/47.integrated_object_subset_by_major_celltypes/Myeloid.rds"
TNK_RDS     <- "/tmp/T_NK.rds"
RDS_OUT     <- "integrated_annotated.rds"

MAC_PATTERN <- "macro|M[φΦ]|Mph|monocyte|myeloid|MDM|TAM"
T_PATTERN   <- "\\bT\\b|T_cell|T.cell|CD4|CD8|Treg|Teff|Tex|Ttrans|Tmem|Tfh|cycling.T|T.cycling|NKT"

classify_cells <- function(meta, source_label) {
  # Find cell type column
  id_col <- NA
  for (col in c("predicted.id", "celltype", "archetype", "cell_type")) {
    if (col %in% colnames(meta)) { id_col <- col; break }
  }

  if (is.na(id_col)) {
    cat("  WARNING: no cell-type column found in", source_label, "— using source label\n")
    meta$cell_type_raw <- source_label
  } else {
    cat("  Using", id_col, "for cell type\n")
    meta$cell_type_raw <- as.character(meta[[id_col]])
  }

  is_mac <- grepl(MAC_PATTERN, meta$cell_type_raw, ignore.case = TRUE, perl = TRUE)
  is_t   <- grepl(T_PATTERN,   meta$cell_type_raw, ignore.case = TRUE, perl = TRUE)
  is_bc  <- meta$cancer == "BC"

  # col1: cell_group
  meta$cell_group <- case_when(
    is_mac &  is_bc  ~ "BC macrophage",
    is_mac & !is_bc  ~ "other cancer macrophage",
    is_t   &  is_bc  ~ "BC T cell",
    is_t   & !is_bc  ~ "other cancer T cell",
    TRUE             ~ "other"
  )

  # col2: bc_er_status — ER status for BC cells only
  meta$bc_er_status <- case_when(
    meta$cancer != "BC" ~ NA_character_,
    "cancer_type_er" %in% colnames(meta) &
      grepl("ER\\+",  meta$cancer_type_er, ignore.case = TRUE) ~ "ER+",
    "cancer_type_er" %in% colnames(meta) &
      grepl("ER-",    meta$cancer_type_er, ignore.case = TRUE) ~ "ER-",
    "cancer_type_er" %in% colnames(meta) &
      grepl("ER\\?",  meta$cancer_type_er, ignore.case = TRUE) ~ "ER?",
    "PAM50" %in% colnames(meta) &
      grepl("LumA|LumB",      meta$PAM50, ignore.case = TRUE)  ~ "ER+",
    "PAM50" %in% colnames(meta) &
      grepl("Basal|TNBC|Her2",meta$PAM50, ignore.case = TRUE)  ~ "ER-",
    "cancer.subtype" %in% colnames(meta) &
      grepl("LumA|LumB|ER\\+",meta$cancer.subtype, ignore.case = TRUE) ~ "ER+",
    "cancer.subtype" %in% colnames(meta) &
      grepl("TNBC|Basal|ER-", meta$cancer.subtype, ignore.case = TRUE) ~ "ER-",
    meta$cancer == "BC" ~ "ER?"   # fallback
  )
  meta
}

# ── Myeloid ──────────────────────────────────────────────────────────────────
cat("Loading Myeloid.rds ...\n")
mye <- readRDS(MYELOID_RDS)
cat("  Cells:", ncol(mye), "\n")
cat("  Metadata cols:", paste(colnames(mye@meta.data), collapse=", "), "\n\n")

if ("predicted.id" %in% colnames(mye@meta.data)) {
  cat("  predicted.id (top 15):\n")
  print(sort(table(mye$predicted.id), decreasing=TRUE)[1:15])
}

mye@meta.data <- classify_cells(mye@meta.data, "Myeloid")
cat("\n  cell_group:\n");    print(table(mye$cell_group))
cat("  bc_er_status:\n"); print(table(mye$bc_er_status, useNA="ifany"))

# ── T_NK ─────────────────────────────────────────────────────────────────────
cat("\nLoading T_NK.rds from /tmp ...\n")
tnk <- readRDS(TNK_RDS)
cat("  Cells:", ncol(tnk), "\n")

if ("predicted.id" %in% colnames(tnk@meta.data)) {
  cat("  predicted.id (top 15):\n")
  print(sort(table(tnk$predicted.id), decreasing=TRUE)[1:15])
}

tnk@meta.data <- classify_cells(tnk@meta.data, "T_NK")
cat("\n  cell_group:\n");    print(table(tnk$cell_group))
cat("  bc_er_status:\n"); print(table(tnk$bc_er_status, useNA="ifany"))

# ── Merge ─────────────────────────────────────────────────────────────────────
cat("\nMerging Myeloid + T_NK ...\n")
merged <- merge(mye, tnk, project = "integrated_annotated")
cat("  Total cells:", ncol(merged), "\n")

cat("\n===== Final metadata summary =====\n")
cat("cell_group:\n");    print(table(merged$cell_group))
cat("\nbc_er_status:\n"); print(table(merged$bc_er_status, useNA="ifany"))

cat("\nSaving to", RDS_OUT, "...\n")
saveRDS(merged, RDS_OUT)
cat("  Done. File size:", round(file.size(RDS_OUT)/1024^2, 1), "MB\n")
cat("\n===== annotate_integrated_obj.R complete =====\n")
