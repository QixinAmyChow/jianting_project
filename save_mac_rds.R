#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(readxl)

RDS_IN   <- "bm_analysis_out/figures_update_3/mac_update3_best.rds"
RDS_OUT  <- "mac_update3_final.rds"
PATIENT_XLS <- "bm_analysis_out/raw_geo/Table S1 patient infor_corrected_2025-08.xlsx"

cat("Loading Seurat object ...\n")
s <- readRDS(RDS_IN)
cat("  Cells:", ncol(s), "| Clusters:", length(levels(s$mac_subcluster)), "\n")

# Add tx_response + gender_clean if not already present
if (!"tx_response" %in% colnames(s@meta.data)) {
  cat("Adding tx_response from patient table ...\n")
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
}

if (!"gender_clean" %in% colnames(s@meta.data)) {
  s$gender_clean <- dplyr::recode(s$gender, "f"="Female", "m"="Male", "na"="Unknown")
}

# Ensure identity is set to mac_subcluster
Idents(s) <- "mac_subcluster"

cat("\nSaving to", RDS_OUT, "...\n")
saveRDS(s, RDS_OUT)
cat("  Done. File size:", round(file.size(RDS_OUT)/1024^2, 1), "MB\n")

cat("\n===== Metadata columns =====\n")
print(colnames(s@meta.data))

cat("\n===== Metadata head (first 6 rows) =====\n")
print(head(s@meta.data))

cat("\n===== Summary of key columns =====\n")
cat("mac_subcluster:\n"); print(table(s$mac_subcluster))
cat("\ncancer_type_er:\n");  print(table(s$cancer_type_er))
cat("\ntx_response:\n");    print(table(s$tx_response))
cat("\ngender_clean:\n");   print(table(s$gender_clean))
cat("\ntissue.origin:\n");  print(table(s$tissue.origin))
