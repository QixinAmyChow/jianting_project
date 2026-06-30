#!/usr/bin/env Rscript
# Load all per-sample Seurat objects, score for macrophages using canonical
# markers, extract macrophage cells, save counts + metadata as CSV/MTX for
# Python/scanpy to load.

library(Seurat)
library(Matrix)

RDS_DIR <- "bm_analysis_out/raw_geo/per_sample_seurat_objects"
OUT_DIR <- "bm_analysis_out/data/macrophages"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MACRO_MARKERS <- c("CD68","CSF1R","MRC1","CD163","FCGR3A","MSR1",
                   "C1QA","C1QB","AIF1","ITGAM")

rds_files <- list.files(RDS_DIR, pattern = "\\.rds$", full.names = TRUE)
cat("Processing", length(rds_files), "samples\n")

all_meta  <- list()
all_counts <- list()

for (f in rds_files) {
    sample_id <- tools::file_path_sans_ext(basename(f))
    cat("  ", sample_id, "...\n")

    s <- tryCatch(readRDS(f), error = function(e) {
        cat("    ERROR loading:", conditionMessage(e), "\n"); NULL
    })
    if (is.null(s)) next

    # Use SCT assay if available, else RNA
    assay <- if ("SCT" %in% names(s@assays)) "SCT" else "RNA"
    DefaultAssay(s) <- assay

    # Score macrophage markers — simple mean of log-normalised expression
    markers_present <- intersect(MACRO_MARKERS, rownames(s))
    if (length(markers_present) == 0) {
        cat("    No macrophage markers found, skipping\n"); next
    }
    expr_mat <- GetAssayData(s, assay = assay, layer = "data")
    macro_score <- colMeans(expr_mat[markers_present, , drop = FALSE])
    s$macro_score1 <- macro_score

    threshold <- quantile(macro_score, 0.75)
    macro_cells <- colnames(s)[macro_score > threshold]
    cat("    Macrophages:", length(macro_cells), "/", ncol(s), "\n")

    if (length(macro_cells) < 10) next

    s_mac <- subset(s, cells = macro_cells)

    # Get normalised counts (SCT or log-normalised RNA)
    mat <- GetAssayData(s_mac, assay = assay, layer = "data")

    meta <- s_mac@meta.data
    meta$sample_id  <- sample_id
    meta$barcode_id <- paste(sample_id, rownames(meta), sep = "__")
    rownames(meta)  <- meta$barcode_id
    colnames(mat)   <- meta$barcode_id

    all_meta[[sample_id]]  <- meta
    all_counts[[sample_id]] <- mat
}

cat("\nMerging", length(all_counts), "samples...\n")

# Intersect genes across all samples
gene_sets <- lapply(all_counts, rownames)
common_genes <- Reduce(intersect, gene_sets)
cat("Common genes:", length(common_genes), "\n")

merged_mat <- do.call(cbind, lapply(all_counts, function(m) m[common_genes, ]))
merged_meta <- do.call(rbind, all_meta)

cat("Total macrophage cells:", ncol(merged_mat), "\n")

# Save
saveRDS(merged_mat, file.path(OUT_DIR, "mac_counts.rds"))
write.csv(merged_meta, file.path(OUT_DIR, "mac_metadata.csv"))
writeMM(merged_mat, file.path(OUT_DIR, "mac_counts.mtx"))
writeLines(rownames(merged_mat), file.path(OUT_DIR, "genes.txt"))
writeLines(colnames(merged_mat), file.path(OUT_DIR, "barcodes.txt"))

cat("Saved to:", OUT_DIR, "\n")
cat("DONE\n")
