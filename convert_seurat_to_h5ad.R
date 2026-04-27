#!/usr/bin/env Rscript
# Convert per-sample Seurat .rds objects to a single merged h5ad
# Retains cell type annotations + key metadata columns

library(Seurat)

rds_dir  <- "bm_analysis_out/raw_geo/per_sample_seurat_objects"
out_dir  <- "bm_analysis_out/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
cat("Found", length(rds_files), "rds files\n")

# Peek at first object to find annotation column name
s0 <- readRDS(rds_files[1])
cat("Metadata columns in first object:\n")
print(colnames(s0@meta.data))

# Collect metadata + write per-sample CSV files for Python to load
meta_list <- list()
for (f in rds_files) {
    sample_name <- tools::file_path_sans_ext(basename(f))
    cat("  Processing:", sample_name, "\n")
    s <- readRDS(f)
    meta <- s@meta.data
    meta$sample_id <- sample_name
    meta$barcode   <- rownames(meta)
    meta_list[[sample_name]] <- meta
}

meta_all <- do.call(rbind, meta_list)
rownames(meta_all) <- NULL
out_csv <- file.path(out_dir, "all_samples_metadata.csv")
write.csv(meta_all, out_csv, row.names = FALSE)
cat("Saved metadata to:", out_csv, "\n")
cat("Total cells:", nrow(meta_all), "\n")
cat("Columns:", paste(colnames(meta_all), collapse=", "), "\n")
